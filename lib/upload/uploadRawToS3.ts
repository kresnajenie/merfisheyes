// Client-side orchestration for the server-side ingestion raw-upload path.
//
// Uploads raw file bytes straight to S3 (no browser parsing) via the
// /api/ingest/* routes, then hands off to the queue (Dataset → QUEUED).
// Uses XMLHttpRequest (not fetch) because only XHR exposes upload byte
// progress (upload.onprogress). Large files (> server threshold) come back
// from initiate as multipart targets; each part is a presigned PUT and its
// ETag is read from the response — this requires the bucket CORS to expose
// the ETag header (ExposeHeaders: ["ETag"]).

import type { ViewerConfig } from "@/lib/utils/viewer-config";

export type RawUploadKind = "single_cell" | "single_molecule";

export interface RawUploadFile {
  /** Relative key under raw/{datasetId}/ (folders keep their path). */
  key: string;
  file: File;
  contentType?: string;
}

export interface RawUploadProgress {
  loaded: number;
  total: number;
  fileIndex: number; // 0-based index of the file currently uploading
  fileCount: number;
}

export type RawUploadPhase = "preparing" | "uploading";

export interface UploadRawParams {
  kind: RawUploadKind;
  title?: string;
  processingParams: Record<string, unknown>;
  /** Owner defaults the new dataset starts with, e.g. { viewMode: "3D" }. */
  viewerConfig?: ViewerConfig;
  files: RawUploadFile[];
  signal?: AbortSignal;
  onProgress?: (p: RawUploadProgress) => void;
  /**
   * Coarse phase signal. "preparing" while the initiate request fetches
   * presigned URLs (no byte progress yet), "uploading" once bytes start
   * flowing. Lets the UI avoid a stalled-looking 0 MB during initiate.
   */
  onPhase?: (phase: RawUploadPhase) => void;
  /**
   * Called with the new dataset id as soon as `initiate` returns it — before
   * any bytes are uploaded. Lets a caller wire up a linked upload / project
   * without waiting for the whole transfer to finish.
   */
  onInitiated?: (datasetId: string) => void;
}

type UploadTarget =
  | { mode: "single"; url: string }
  | {
      mode: "multipart";
      s3UploadId: string;
      partSize: number;
      parts: Array<{ partNumber: number; url: string }>;
    };

async function readError(res: Response): Promise<string> {
  try {
    const j = await res.json();

    return j?.message || j?.error || res.statusText;
  } catch {
    return res.statusText || `HTTP ${res.status}`;
  }
}

/** PUT a body to a presigned URL via XHR, reporting byte progress. */
function xhrPut(
  url: string,
  body: Blob,
  contentType: string | null,
  onProgress: (loaded: number) => void,
  signal?: AbortSignal,
): Promise<{ etag: string | null }> {
  return new Promise((resolve, reject) => {
    if (signal?.aborted) {
      reject(new DOMException("Aborted", "AbortError"));

      return;
    }

    const xhr = new XMLHttpRequest();

    xhr.open("PUT", url);
    // Only set Content-Type for single-PUT (its presigned URL is signed with
    // it). Multipart UploadPart URLs are not signed with Content-Type.
    if (contentType) xhr.setRequestHeader("Content-Type", contentType);

    const onAbort = () => xhr.abort();

    signal?.addEventListener("abort", onAbort);
    const cleanup = () => signal?.removeEventListener("abort", onAbort);

    xhr.upload.onprogress = (e) => {
      if (e.lengthComputable) onProgress(e.loaded);
    };
    xhr.onload = () => {
      cleanup();
      if (xhr.status >= 200 && xhr.status < 300) {
        resolve({ etag: xhr.getResponseHeader("ETag") });
      } else {
        reject(new Error(`Upload failed with status ${xhr.status}`));
      }
    };
    xhr.onerror = () => {
      cleanup();
      reject(new Error("Network error during upload"));
    };
    xhr.onabort = () => {
      cleanup();
      reject(new DOMException("Aborted", "AbortError"));
    };
    xhr.send(body);
  });
}

async function markFileComplete(
  datasetId: string,
  uploadId: string,
  key: string,
  multipart:
    | { s3UploadId: string; parts: Array<{ partNumber: number; etag: string }> }
    | undefined,
  signal?: AbortSignal,
): Promise<void> {
  const res = await fetch(
    `/api/ingest/${datasetId}/files/${encodeURIComponent(key)}/complete`,
    {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(
        multipart
          ? {
              uploadId,
              s3UploadId: multipart.s3UploadId,
              parts: multipart.parts,
            }
          : { uploadId },
      ),
      signal,
    },
  );

  if (!res.ok) {
    throw new Error(`Marking ${key} complete failed: ${await readError(res)}`);
  }
}

/**
 * Parts of one file in flight at once. One connection to S3 carries only a
 * few MB/s from a lab network (3.3 MB/s, uploading Beaumont set1 from the
 * lab), while many at once carry what the link has (26–29 MB/s from the same
 * machine): a 6 GB table went up part after part in about half an hour.
 */
export const PART_CONCURRENCY = 6;

/**
 * Upload a multipart target's parts, `PART_CONCURRENCY` at a time, and give
 * their ETags in part order. `onLoaded` gets the bytes of the parts still in
 * flight, `onPartDone` each finished part's size. A part that fails stops the
 * others: the upload has failed.
 */
export async function uploadParts(
  target: {
    partSize: number;
    parts: Array<{ partNumber: number; url: string }>;
  },
  file: Blob,
  key: string,
  onLoaded: (inFlightBytes: number) => void,
  onPartDone: (bytes: number) => void,
  signal?: AbortSignal,
  put: typeof xhrPut = xhrPut,
): Promise<Array<{ partNumber: number; etag: string }>> {
  const etags: Array<{ partNumber: number; etag: string }> = [];
  const inFlight = new Map<number, number>();
  const loaded = () => {
    let sum = 0;

    inFlight.forEach((bytes) => {
      sum += bytes;
    });
    onLoaded(sum);
  };
  const stop = new AbortController();
  const onAbort = () => stop.abort();

  signal?.addEventListener("abort", onAbort);
  if (signal?.aborted) stop.abort();

  let next = 0;
  const worker = async () => {
    while (next < target.parts.length) {
      const p = next++;
      const start = p * target.partSize;
      const end = Math.min(start + target.partSize, file.size);

      inFlight.set(p, 0);
      const { etag } = await put(
        target.parts[p].url,
        file.slice(start, end),
        null,
        (bytes) => {
          inFlight.set(p, bytes);
          loaded();
        },
        stop.signal,
      );

      if (!etag) {
        throw new Error(
          `Missing ETag for part ${p + 1} of ${key}. Ensure the S3 bucket CORS exposes the ETag header.`,
        );
      }
      etags.push({ partNumber: target.parts[p].partNumber, etag });
      inFlight.delete(p);
      onPartDone(end - start);
      loaded();
    }
  };

  try {
    await Promise.all(
      Array.from(
        { length: Math.min(PART_CONCURRENCY, target.parts.length) },
        () =>
          worker().catch((e) => {
            stop.abort();
            throw e;
          }),
      ),
    );
  } finally {
    signal?.removeEventListener("abort", onAbort);
  }

  return etags.sort((a, b) => a.partNumber - b.partNumber);
}

export interface RawUploadResult {
  datasetId: string;
  /** Null when the processing job could not be submitted — see submitError. */
  batchJobId: string | null;
  /**
   * Set when the raw upload succeeded but the processing job was NOT submitted
   * (e.g. Batch misconfigured). The dataset stays QUEUED and can be retried —
   * the UI must not present this as a clean success.
   */
  submitError: string | null;
}

/**
 * Upload raw files to S3 and queue the dataset for server-side processing.
 * Resolves once the dataset reaches QUEUED.
 */
export async function uploadRawToS3(
  params: UploadRawParams,
): Promise<RawUploadResult> {
  const {
    kind,
    title,
    processingParams,
    viewerConfig,
    files,
    signal,
    onProgress,
    onPhase,
    onInitiated,
  } = params;

  if (files.length === 0) throw new Error("No files to upload");

  const total = files.reduce((sum, f) => sum + f.file.size, 0);
  let completedBytes = 0;

  const report = (currentLoaded: number, fileIndex: number) =>
    onProgress?.({
      loaded: completedBytes + currentLoaded,
      total,
      fileIndex,
      fileCount: files.length,
    });

  // 1. Initiate: create the Dataset (UPLOADING) + presigned upload targets.
  // This round-trip has no byte progress, so flag it as "preparing".
  onPhase?.("preparing");
  const initRes = await fetch("/api/ingest/initiate", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({
      metadata: { title },
      processingParams: { ...processingParams, kind },
      viewerConfig,
      files: files.map((f) => ({
        key: f.key,
        size: f.file.size,
        contentType: f.contentType || "application/octet-stream",
      })),
    }),
    signal,
  });

  if (!initRes.ok) {
    throw new Error(`Could not start upload: ${await readError(initRes)}`);
  }

  const { datasetId, uploadId, uploadUrls } = await initRes.json();

  onInitiated?.(datasetId);

  onPhase?.("uploading");
  report(0, 0);

  // 2. Upload each file (sequentially for simple, correct aggregate progress).
  for (let i = 0; i < files.length; i++) {
    const { key, file } = files[i];
    const contentType = files[i].contentType || "application/octet-stream";
    const target: UploadTarget = uploadUrls[key];

    if (!target) throw new Error(`Server returned no upload URL for ${key}`);

    if (target.mode === "single") {
      await xhrPut(
        target.url,
        file,
        contentType,
        (loaded) => report(loaded, i),
        signal,
      );
      completedBytes += file.size;
      await markFileComplete(datasetId, uploadId, key, undefined, signal);
    } else {
      const etags = await uploadParts(
        target,
        file,
        key,
        (loaded) => report(loaded, i),
        (bytes) => {
          completedBytes += bytes;
        },
        signal,
      );

      await markFileComplete(
        datasetId,
        uploadId,
        key,
        { s3UploadId: target.s3UploadId, parts: etags },
        signal,
      );
    }

    report(0, Math.min(i + 1, files.length - 1));
  }

  // 3. Complete: verify all files uploaded → QUEUED.
  const compRes = await fetch(`/api/ingest/${datasetId}/complete`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ uploadId }),
    signal,
  });

  if (!compRes.ok) {
    throw new Error(`Could not finalize upload: ${await readError(compRes)}`);
  }

  const completed = await compRes.json().catch(() => ({}));

  return {
    datasetId,
    batchJobId: completed?.batchJobId ?? null,
    submitError: completed?.submitError ?? null,
  };
}

export type IngestStatusValue =
  | "UPLOADING"
  | "QUEUED"
  | "PROCESSING"
  | "COMPLETE"
  | "FAILED";

export interface IngestStatus {
  status: IngestStatusValue;
  progress?: { stage?: string; percent?: number } | null;
  errorMessage?: string | null;
  viewerUrl?: string | null;
  numCells?: number;
  numGenes?: number;
}

/**
 * Poll server-side processing status until it reaches a terminal state.
 * Returns a stop function; safe to call after the user dismisses the UI.
 *
 * (Supabase Realtime can replace this later — see design §5.5 — but polling
 * needs no extra infrastructure and works everywhere.)
 */
export function pollIngestStatus(
  datasetId: string,
  onUpdate: (status: IngestStatus) => void,
  opts: { intervalMs?: number } = {},
): () => void {
  const intervalMs = opts.intervalMs ?? 3000;
  let stopped = false;
  let timer: ReturnType<typeof setTimeout> | undefined;

  const tick = async () => {
    if (stopped) return;
    try {
      const res = await fetch(`/api/ingest/${datasetId}/status`);

      if (res.ok) {
        const data: IngestStatus = await res.json();

        if (stopped) return;
        onUpdate(data);
        if (data.status === "COMPLETE" || data.status === "FAILED") {
          stopped = true;

          return;
        }
      }
    } catch {
      // Transient network errors shouldn't kill the poll loop.
    }
    if (!stopped) timer = setTimeout(tick, intervalMs);
  };

  tick();

  return () => {
    stopped = true;
    if (timer) clearTimeout(timer);
  };
}

/**
 * Best-effort cancel: tell the server to abort the multipart uploads, delete
 * any raw objects, and remove the dataset row. Fire-and-forget from the UI.
 */
export async function abortRawUpload(datasetId: string): Promise<void> {
  try {
    await fetch(`/api/ingest/${datasetId}/abort`, { method: "POST" });
  } catch {
    // Cancellation is best-effort; the S3 lifecycle rule is the backstop.
  }
}
