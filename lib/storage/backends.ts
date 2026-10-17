/**
 * What an ElementStore reads from: a public URL, a dataset in the app's
 * bucket (through presigned URLs), or files dropped in the browser. Each
 * answers a missing object with undefined — a public bucket that can't be
 * listed says 403 for one — and takes the read's abort signal to `fetch`.
 */
import type { Backend } from "./ElementStore";
import type { RangeQuery } from "@zarrita/storage";

import { fetchRange } from "./fetchRange";

const missing = (status: number) => status === 404 || status === 403;

/** Objects under a public URL; keys are paths from `root`. */
export function urlBackend(root: string): Backend {
  const base = root.replace(/\/+$/, "");

  return {
    async get(key, signal) {
      const res = await fetch(base + key, { signal });

      if (missing(res.status)) return undefined;
      if (!res.ok) throw new Error(`GET ${res.status} for ${key}`);

      return new Uint8Array(await res.arrayBuffer());
    },
    getRange: (key, range, signal) => fetchRange(base + key, range, signal),
    async open(key, signal) {
      const res = await fetch(base + key, { signal });

      if (missing(res.status)) return undefined;
      if (!res.ok) throw new Error(`GET ${res.status} for ${key}`);

      return res;
    },
  };
}

const URL_TTL_MARGIN_MS = 60_000;
// Presigned URLs by dataset and key, shared by every store on the dataset,
// with the presign call in flight shared too
const presigned = new Map<
  string,
  Promise<{ url: string; expiresAt: number } | null>
>();

/** A presigned URL for an object of an app dataset; null when it has none. */
export async function presignedUrl(
  datasetId: string,
  key: string,
  apiBase = "",
): Promise<string | null> {
  const id = `${datasetId}/${key}`;
  const cached = await presigned.get(id);

  if (cached && cached.expiresAt - Date.now() > URL_TTL_MARGIN_MS) {
    return cached.url;
  }
  const p = (async () => {
    const res = await fetch(
      `${apiBase}/api/datasets/${datasetId}/object?key=${encodeURIComponent(key)}`,
    );

    if (res.status === 404) return null;
    if (!res.ok) throw new Error(`/object failed (${res.status}) for ${key}`);
    const { url, expiresIn } = (await res.json()) as {
      url: string;
      expiresIn: number;
    };

    return { url, expiresAt: Date.now() + expiresIn * 1000 };
  })();

  presigned.set(id, p);
  p.catch(() => presigned.delete(id));

  return (await p)?.url ?? null;
}

/**
 * Objects of a dataset in the app's bucket, under `datasets/{id}/{prefix}`;
 * keys are paths from the prefix.
 */
export function datasetBackend(datasetId: string, prefix = ""): Backend {
  const full = (key: string) =>
    prefix.replace(/^\/+|\/+$/g, "") +
    (prefix.replace(/^\/+|\/+$/g, "") ? "/" : "") +
    key.replace(/^\/+/, "");

  return {
    async get(key, signal) {
      const url = await presignedUrl(datasetId, full(key));

      if (!url) return undefined;
      const res = await fetch(url, { signal });

      if (missing(res.status)) return undefined;
      if (!res.ok) throw new Error(`S3 GET ${res.status} for ${full(key)}`);

      return new Uint8Array(await res.arrayBuffer());
    },
    async getRange(key, range, signal) {
      const url = await presignedUrl(datasetId, full(key));

      return url ? fetchRange(url, range, signal) : undefined;
    },
    async open(key, signal) {
      const url = await presignedUrl(datasetId, full(key));

      if (!url) return undefined;
      const res = await fetch(url, { signal });

      if (missing(res.status)) return undefined;
      if (!res.ok) throw new Error(`S3 GET ${res.status} for ${full(key)}`);

      return res;
    },
  };
}

/** Files picked or dropped in the browser, keyed by path from their root. */
export function fileBackend(files: Map<string, File>): Backend {
  const find = (key: string) => files.get(key.replace(/^\/+/, ""));

  return {
    async get(key) {
      const file = find(key);

      return file && new Uint8Array(await file.arrayBuffer());
    },
    async getRange(key, range: RangeQuery) {
      const file = find(key);

      if (!file) return undefined;
      const part =
        "suffixLength" in range
          ? file.slice(-range.suffixLength)
          : file.slice(range.offset, range.offset + range.length);

      return new Uint8Array(await part.arrayBuffer());
    },
    async open(key) {
      const file = find(key);

      return file && new Response(file.stream());
    },
  };
}
