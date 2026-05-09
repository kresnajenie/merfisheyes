/**
 * Upload a zarr folder (Map<relativePath, File>) to S3 using a single
 * presigned-POST policy. The server returns ONE policy that authorizes
 * uploads to `datasets/{datasetId}/...`; we POST each file against it
 * with bounded concurrency.
 *
 * Progress is tracked at the file level (no per-byte progress — that
 * requires XMLHttpRequest, which we don't need for v1).
 */

const MAX_CONCURRENCY = 8;
const MAX_RETRIES = 3;
const PROGRESS_BATCH_SIZE = 25; // flush completedFiles to server every N files

interface InitiateZarrResponse {
  datasetId: string;
  uploadId: string;
  keyPrefix: string; // e.g. "datasets/ds_abc123/"
  url: string; // S3 POST endpoint
  fields: Record<string, string>;
}

export interface ZarrUploadOptions {
  fileMap: Map<string, File>; // relativePath -> File
  fingerprint: string;
  metadata: {
    title?: string;
    numCells: number;
    numGenes: number;
    platform?: string;
  };
  email: string;
  datasetName: string;
  onProgress?: (progress: number, message: string) => void;
}

export interface ZarrUploadResult {
  datasetId: string;
  uploadId: string;
}

/**
 * 1) Initiate upload (server signs one POST policy, returns datasetId + uploadId).
 * 2) POST each file in fileMap with bounded concurrency + retry.
 * 3) Periodically flush per-file progress to the server.
 * 4) Mark dataset complete.
 */
export async function uploadZarrToS3(
  opts: ZarrUploadOptions,
): Promise<ZarrUploadResult> {
  const { fileMap, fingerprint, metadata, email, datasetName, onProgress } =
    opts;

  if (fileMap.size === 0) {
    throw new Error("Zarr fileMap is empty");
  }

  onProgress?.(0, "Initiating upload...");

  // 1. Initiate
  const initRes = await fetch("/api/datasets/initiate-zarr", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({
      fingerprint,
      metadata,
      totalFiles: fileMap.size,
    }),
  });

  if (!initRes.ok) {
    const err = await initRes.json().catch(() => ({}));

    throw new Error(err.error || `initiate-zarr failed: ${initRes.status}`);
  }

  const init: InitiateZarrResponse = await initRes.json();
  const { datasetId, uploadId, keyPrefix, url, fields } = init;

  console.log(
    `[uploadZarrToS3] datasetId=${datasetId} uploadId=${uploadId} files=${fileMap.size}`,
  );

  // 2. Upload with bounded concurrency
  const entries = Array.from(fileMap.entries());
  const total = entries.length;
  let completed = 0;
  let progressBuffer = 0;

  const flushProgress = async (force = false) => {
    if (progressBuffer === 0) return;
    if (!force && progressBuffer < PROGRESS_BATCH_SIZE) return;
    const delta = progressBuffer;

    progressBuffer = 0;
    try {
      await fetch(`/api/datasets/${datasetId}/upload-progress`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ uploadId, delta }),
      });
    } catch (e) {
      // Progress flushing is best-effort; final flush at the end ensures
      // the server's count matches before /complete is called.
      console.warn("[uploadZarrToS3] progress flush failed:", e);
      progressBuffer += delta; // restore so the next flush retries
    }
  };

  const uploadOne = async (relPath: string, file: File): Promise<void> => {
    let attempt = 0;
    let lastErr: unknown = null;

    while (attempt < MAX_RETRIES) {
      try {
        const fd = new FormData();

        for (const [k, v] of Object.entries(fields)) fd.append(k, v);
        fd.append("key", `${keyPrefix}${relPath}`);
        fd.append("file", file, relPath);

        const res = await fetch(url, { method: "POST", body: fd });

        if (!res.ok) {
          const text = await res.text().catch(() => "");

          throw new Error(`S3 POST ${res.status}: ${text.slice(0, 200)}`);
        }

        return;
      } catch (e) {
        lastErr = e;
        attempt++;
        if (attempt < MAX_RETRIES) {
          const backoff = 500 * Math.pow(4, attempt - 1);

          await new Promise((r) => setTimeout(r, backoff));
        }
      }
    }
    throw new Error(
      `Failed to upload ${relPath} after ${MAX_RETRIES} attempts: ${lastErr}`,
    );
  };

  // Worker pool: keep MAX_CONCURRENCY uploads in flight at all times.
  let cursor = 0;

  const worker = async () => {
    while (true) {
      const idx = cursor++;

      if (idx >= entries.length) return;
      const [relPath, file] = entries[idx];

      await uploadOne(relPath, file);
      completed++;
      progressBuffer++;
      const pct = (completed / total) * 100;

      onProgress?.(pct, `Uploaded ${completed}/${total} files`);
      await flushProgress(false);
    }
  };

  const workers = Array.from(
    { length: Math.min(MAX_CONCURRENCY, entries.length) },
    () => worker(),
  );

  await Promise.all(workers);
  await flushProgress(true);

  // 3. Mark dataset complete
  onProgress?.(99, "Finalizing...");
  const completeRes = await fetch(`/api/datasets/${datasetId}/complete`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ uploadId }),
  });

  if (!completeRes.ok) {
    const err = await completeRes.json().catch(() => ({}));

    throw new Error(err.error || `complete failed: ${completeRes.status}`);
  }

  // 4. Email notification (best-effort)
  try {
    await fetch("/api/send-email", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        email,
        datasetName,
        datasetId,
        metadata: {
          numCells: metadata.numCells,
          numGenes: metadata.numGenes,
          platform: metadata.platform || "h5ad-zarr",
          formatVersion: "zarr",
        },
      }),
    });
  } catch (e) {
    console.warn("[uploadZarrToS3] email notification failed:", e);
  }

  onProgress?.(100, "Upload complete");

  return { datasetId, uploadId };
}
