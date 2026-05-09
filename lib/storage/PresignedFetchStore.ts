import type { AbsolutePath, AsyncReadable } from "@zarrita/storage";

/**
 * Zarrita store backed by a private S3 bucket through our API.
 *
 * For each `get(key)`:
 *   1. Fetch a presigned GET URL from `/api/datasets/{id}/object?key=...`
 *   2. Follow the URL → bytes
 *
 * Presigned URLs are cached in-memory until shortly before they expire so we
 * only round-trip through our API once per key per session.
 *
 * 404 (object missing) returns `undefined` per the AsyncReadable contract,
 * so zarr `open()` can probe for `.zgroup` / `.zarray` without throwing.
 */
export class PresignedFetchStore implements AsyncReadable {
  private datasetId: string;
  private apiBase: string;
  private urlCache = new Map<string, { url: string; expiresAt: number }>();
  // Sub-second margin to avoid using a URL that's about to expire mid-fetch.
  private readonly safetyMarginMs = 60_000;

  constructor(datasetId: string, apiBase = "") {
    this.datasetId = datasetId;
    this.apiBase = apiBase;
  }

  async get(key: AbsolutePath): Promise<Uint8Array | undefined> {
    const relativeKey = key.startsWith("/") ? key.slice(1) : key;
    const presigned = await this.getPresignedUrl(relativeKey);

    if (!presigned) return undefined;

    const res = await fetch(presigned);

    if (res.status === 404 || res.status === 403) return undefined;
    if (!res.ok) {
      throw new Error(
        `S3 GET ${res.status} for ${relativeKey}: ${res.statusText}`,
      );
    }
    const buf = await res.arrayBuffer();

    return new Uint8Array(buf);
  }

  private async getPresignedUrl(relativeKey: string): Promise<string | null> {
    const now = Date.now();
    const cached = this.urlCache.get(relativeKey);

    if (cached && cached.expiresAt - now > this.safetyMarginMs) {
      return cached.url;
    }

    const apiUrl = `${this.apiBase}/api/datasets/${this.datasetId}/object?key=${encodeURIComponent(relativeKey)}`;
    const res = await fetch(apiUrl);

    if (res.status === 404) return null;
    if (!res.ok) {
      throw new Error(`/object failed (${res.status}) for ${relativeKey}`);
    }

    const { url, expiresIn } = (await res.json()) as {
      url: string;
      expiresIn: number;
    };

    this.urlCache.set(relativeKey, {
      url,
      expiresAt: now + expiresIn * 1000,
    });

    return url;
  }
}

/**
 * One-shot helper: fetch the dataset's full key list from /list.
 * The zarr adapter needs this to enumerate obs/obsm/var children since
 * AsyncReadable has no LIST op.
 */
export async function fetchDatasetKeyList(
  datasetId: string,
  apiBase = "",
): Promise<string[]> {
  const res = await fetch(`${apiBase}/api/datasets/${datasetId}/list`);

  if (!res.ok) {
    throw new Error(`/list failed (${res.status}) for ${datasetId}`);
  }
  const { keys } = (await res.json()) as { keys: string[] };

  return keys;
}
