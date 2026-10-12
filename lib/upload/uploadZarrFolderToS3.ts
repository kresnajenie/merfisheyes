/**
 * Upload a zarr folder to S3 as a dataset, as it is: the viewer reads the
 * zarr directly, so nothing is processed. Used for an OME-Zarr image and for
 * a SpatialData store.
 *
 * Thousands of small chunk files go up under ONE presigned-POST policy scoped
 * to `datasets/{datasetId}/`, eight at a time, and the policy is re-signed
 * shortly before it expires (or on a 403) so a multi-gigabyte folder is not
 * cut off after an hour. `lastKey` (the root `zarr.json`) goes last: the
 * dataset only becomes openable once everything it lists is in place.
 */

const MAX_CONCURRENCY = 8;
const MAX_RETRIES = 3;
const PROGRESS_BATCH_SIZE = 25; // files counted on the server per request
const POLICY_REFRESH_MARGIN_MS = 5 * 60 * 1000;

export interface ZarrUploadFile {
  /** Key under the dataset's prefix, e.g. "zarr.json" or "0/c/0/3/7". */
  key: string;
  file: File;
}

export interface ZarrUploadProgress {
  loaded: number;
  total: number;
  fileIndex: number;
  fileCount: number;
}

export interface ZarrFolderUploadOptions {
  title: string;
  /** "image" or "spatialdata": becomes the dataset's type. */
  platform: "image" | "spatialdata";
  numCells?: number;
  numGenes?: number;
  files: ZarrUploadFile[];
  /** The file that makes the dataset openable; uploaded after all others. */
  lastKey: string;
  signal?: AbortSignal;
  /** Called with the new dataset id before any bytes are uploaded. */
  onInitiated?: (datasetId: string) => void;
  onPhase?: (phase: "preparing" | "uploading") => void;
  onProgress?: (p: ZarrUploadProgress) => void;
}

interface UploadPolicy {
  url: string;
  fields: Record<string, string>;
  keyPrefix: string;
  expiresAt: string;
}

async function sha256Hex(text: string): Promise<string> {
  const digest = await crypto.subtle.digest(
    "SHA-256",
    new TextEncoder().encode(text),
  );

  return Array.from(new Uint8Array(digest), (b) =>
    b.toString(16).padStart(2, "0"),
  ).join("");
}

export async function uploadZarrFolderToS3(
  opts: ZarrFolderUploadOptions,
): Promise<{ datasetId: string }> {
  const { title, files, signal, onInitiated, onPhase, onProgress } = opts;
  const root = files.find((f) => f.key === opts.lastKey);

  if (!root) throw new Error(`The folder has no ${opts.lastKey}`);

  const total = files.reduce((sum, f) => sum + f.file.size, 0);

  onPhase?.("preparing");
  // The same folder uploaded twice gets the same fingerprint: its metadata
  // plus the size of every file.
  const fingerprint = await sha256Hex(
    (await root.file.text()) +
      files
        .map((f) => `${f.key}:${f.file.size}`)
        .sort()
        .join("\n"),
  );
  const initRes = await fetch("/api/datasets/initiate-zarr", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({
      fingerprint,
      metadata: {
        title,
        numCells: opts.numCells ?? 0,
        numGenes: opts.numGenes ?? 0,
        platform: opts.platform,
      },
      totalFiles: files.length,
    }),
    signal,
  });

  if (!initRes.ok) {
    const err = await initRes.json().catch(() => ({}));

    throw new Error(err.error || `Couldn't start the upload`);
  }
  const init = await initRes.json();
  const datasetId: string = init.datasetId;
  const uploadId: string = init.uploadId;
  let policy: UploadPolicy = init.policy;

  onInitiated?.(datasetId);

  // One refresh at a time, shared by every upload in flight.
  let refreshing: Promise<void> | null = null;
  const refreshPolicy = (): Promise<void> => {
    refreshing ??= (async () => {
      const res = await fetch(`/api/datasets/${datasetId}/upload-policy`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ uploadId }),
        signal,
      });

      if (!res.ok) throw new Error(`Policy refresh failed: ${res.status}`);
      policy = (await res.json()).policy;
    })().finally(() => {
      refreshing = null;
    });

    return refreshing;
  };

  const postFile = async ({ key, file }: ZarrUploadFile) => {
    let lastErr: unknown = null;

    for (let attempt = 1; attempt <= MAX_RETRIES; attempt++) {
      signal?.throwIfAborted();
      try {
        if (
          Date.parse(policy.expiresAt) - Date.now() <
          POLICY_REFRESH_MARGIN_MS
        ) {
          await refreshPolicy();
        }
        const fd = new FormData();

        // `fields` carries a placeholder `key`; S3 rejects two `key` parts.
        for (const [k, v] of Object.entries(policy.fields)) {
          if (k !== "key") fd.append(k, v);
        }
        fd.append("key", `${policy.keyPrefix}${key}`);
        fd.append(
          "Content-Type",
          key.endsWith(".json")
            ? "application/json"
            : "application/octet-stream",
        );
        fd.append("file", file, key.split("/").pop());

        const res = await fetch(policy.url, {
          method: "POST",
          body: fd,
          signal,
        });

        if (res.ok) return;
        // 403 is what S3 answers for an expired policy
        if (res.status === 403) await refreshPolicy();
        lastErr = new Error(`S3 POST ${res.status}`);
      } catch (e) {
        if (signal?.aborted) throw e;
        lastErr = e;
      }
      if (attempt < MAX_RETRIES) {
        await new Promise((r) => setTimeout(r, 500 * 4 ** (attempt - 1)));
      }
    }
    throw new Error(`Upload of ${key} failed: ${lastErr}`);
  };

  // The server counts files, not bytes; told in batches.
  let uncounted = 0;
  const flushCount = async (force: boolean) => {
    if (uncounted === 0 || (!force && uncounted < PROGRESS_BATCH_SIZE)) return;
    const delta = uncounted;

    uncounted = 0;
    const res = await fetch(`/api/datasets/${datasetId}/upload-progress`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ uploadId, delta }),
      signal,
    }).catch(() => null);

    if (!res?.ok) uncounted += delta; // counted again with the next batch
  };

  let loaded = 0;
  let done = 0;
  const finish = async (f: ZarrUploadFile) => {
    await postFile(f);
    loaded += f.file.size;
    done++;
    uncounted++;
    onProgress?.({ loaded, total, fileIndex: done, fileCount: files.length });
    await flushCount(false);
  };

  onPhase?.("uploading");
  const rest = files.filter((f) => f !== root);
  let cursor = 0;

  await Promise.all(
    Array.from({ length: Math.min(MAX_CONCURRENCY, rest.length) }, async () => {
      while (cursor < rest.length) await finish(rest[cursor++]);
    }),
  );
  await finish(root);
  // The count has to be complete before the dataset can be marked so
  for (let i = 0; uncounted > 0 && i < MAX_RETRIES; i++) await flushCount(true);

  const completeRes = await fetch(`/api/datasets/${datasetId}/complete`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ uploadId }),
    signal,
  });

  if (!completeRes.ok) {
    const err = await completeRes.json().catch(() => ({}));

    throw new Error(err.error || `Couldn't finish the upload`);
  }

  return { datasetId };
}
