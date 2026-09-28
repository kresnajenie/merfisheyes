/**
 * Upload a chunked single-cell dataset (manifest.json, coords/, expr/, obs/,
 * palettes/, de/) to S3 using ONE presigned-POST policy scoped to
 * `datasets/{datasetId}/`, instead of one presigned PUT URL per file.
 *
 * - Files go up with bounded concurrency.
 * - The policy is re-signed via /upload-policy shortly before it expires (or
 *   on a 403), so uploads longer than an hour keep going.
 * - UploadFile rows are still created at initiate (the viewer lists a
 *   chunked dataset's files from them) and are marked COMPLETE in batches.
 * - manifest.json is uploaded last, once everything it points at is in place.
 *
 * Returns once every file is uploaded and registered; the caller then calls
 * /api/datasets/{id}/complete.
 */

const MAX_CONCURRENCY = 8;
const MAX_RETRIES = 3;
const COMPLETE_BATCH_SIZE = 100; // UploadFile rows marked COMPLETE per request
const POLICY_REFRESH_MARGIN_MS = 5 * 60 * 1000;
const MANIFEST_KEY = "manifest.json";

export interface ChunkedUploadFile {
  key: string; // path relative to the dataset root, e.g. "expr/chunk_00000.bin.gz"
  blob: Blob;
  size: number;
  contentType: string;
}

interface UploadPolicy {
  url: string;
  fields: Record<string, string>;
  keyPrefix: string; // "datasets/{datasetId}/"
  expiresAt: string; // ISO
}

export interface ChunkedUploadOptions {
  fingerprint: string;
  metadata: {
    title: string;
    numCells: number;
    numGenes: number;
    platform?: string;
    description?: string;
  };
  files: ChunkedUploadFile[];
  asAdmin?: boolean;
  onProgress?: (progress: number, message: string) => void;
}

export interface ChunkedUploadResult {
  datasetId: string;
  uploadId: string;
}

async function withRetries<T>(label: string, fn: () => Promise<T>): Promise<T> {
  let lastErr: unknown = null;

  for (let attempt = 1; attempt <= MAX_RETRIES; attempt++) {
    try {
      return await fn();
    } catch (e) {
      lastErr = e;
      if (attempt < MAX_RETRIES) {
        await new Promise((r) => setTimeout(r, 500 * Math.pow(4, attempt - 1)));
      }
    }
  }
  throw new Error(`${label} failed after ${MAX_RETRIES} attempts: ${lastErr}`);
}

export async function uploadChunkedToS3(
  opts: ChunkedUploadOptions,
): Promise<ChunkedUploadResult> {
  const { fingerprint, metadata, files, asAdmin = false, onProgress } = opts;

  const manifest = files.find((f) => f.key === MANIFEST_KEY);

  if (!manifest) {
    throw new Error(`Upload is missing ${MANIFEST_KEY}`);
  }
  const dataFiles = files.filter((f) => f.key !== MANIFEST_KEY);

  // 1. Initiate: creates Dataset/UploadSession/UploadFile rows, returns a policy.
  const initRes = await fetch("/api/datasets/initiate", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({
      fingerprint,
      metadata,
      files: files.map((f) => ({
        key: f.key,
        size: f.size,
        contentType: f.contentType,
      })),
      asAdmin,
      uploadMode: "post-policy",
    }),
  });

  if (!initRes.ok) {
    const err = await initRes.json().catch(() => ({}));

    throw new Error(
      err.error || `Failed to initiate upload: ${initRes.status}`,
    );
  }

  const init: { datasetId: string; uploadId: string; policy: UploadPolicy } =
    await initRes.json();
  const { datasetId, uploadId } = init;
  let policy = init.policy;

  // Single-flight policy refresh shared by all workers.
  let refreshing: Promise<void> | null = null;
  const refreshPolicy = (): Promise<void> => {
    refreshing ??= (async () => {
      const res = await fetch(`/api/datasets/${datasetId}/upload-policy`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ uploadId }),
      });

      if (!res.ok) {
        const err = await res.json().catch(() => ({}));

        throw new Error(err.error || `Policy refresh failed: ${res.status}`);
      }
      policy = (await res.json()).policy;
    })().finally(() => {
      refreshing = null;
    });

    return refreshing;
  };

  const ensureFreshPolicy = async () => {
    if (Date.parse(policy.expiresAt) - Date.now() < POLICY_REFRESH_MARGIN_MS) {
      await refreshPolicy();
    }
  };

  const postFile = async (file: ChunkedUploadFile): Promise<void> => {
    await withRetries(`Upload of ${file.key}`, async () => {
      await ensureFreshPolicy();

      const fd = new FormData();

      // `fields` carries a placeholder `key`; S3 rejects two `key` parts.
      for (const [k, v] of Object.entries(policy.fields)) {
        if (k === "key") continue;
        fd.append(k, v);
      }
      fd.append("key", `${policy.keyPrefix}${file.key}`);
      fd.append("Content-Type", file.contentType || "application/octet-stream");
      fd.append("file", file.blob, file.key.split("/").pop());

      const res = await fetch(policy.url, { method: "POST", body: fd });

      if (!res.ok) {
        // 403 is what S3 returns for an expired policy — re-sign before the retry.
        if (res.status === 403) await refreshPolicy();
        const text = await res.text().catch(() => "");

        throw new Error(`S3 POST ${res.status}: ${text.slice(0, 200)}`);
      }
    });
  };

  // Batched "mark COMPLETE" for UploadFile rows.
  const pendingKeys: string[] = [];
  const flushCompleted = async (force: boolean) => {
    while (
      pendingKeys.length >= COMPLETE_BATCH_SIZE ||
      (force && pendingKeys.length > 0)
    ) {
      const batch = pendingKeys.splice(0, COMPLETE_BATCH_SIZE);

      await withRetries("Marking files complete", async () => {
        const res = await fetch(`/api/datasets/${datasetId}/files/complete`, {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({ uploadId, fileKeys: batch }),
        });

        if (!res.ok) throw new Error(`files/complete ${res.status}`);
      });
    }
  };

  // 2. Data files with bounded concurrency.
  const total = files.length;
  let completed = 0;
  const reportProgress = () =>
    onProgress?.(
      (completed / total) * 100,
      `Uploaded ${completed}/${total} files`,
    );

  let cursor = 0;
  const worker = async () => {
    while (cursor < dataFiles.length) {
      const file = dataFiles[cursor++];

      await postFile(file);
      pendingKeys.push(file.key);
      completed++;
      reportProgress();
      await flushCompleted(false);
    }
  };

  await Promise.all(
    Array.from({ length: Math.min(MAX_CONCURRENCY, dataFiles.length) }, worker),
  );
  await flushCompleted(true);

  // 3. Manifest last: the dataset only becomes loadable once everything it
  // references is already in S3.
  await postFile(manifest);
  pendingKeys.push(manifest.key);
  completed++;
  reportProgress();
  await flushCompleted(true);

  return { datasetId, uploadId };
}
