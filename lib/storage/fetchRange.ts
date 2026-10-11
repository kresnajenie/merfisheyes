import type { RangeQuery } from "@zarrita/storage";

/**
 * GET part of an object over HTTP, for zarr stores whose arrays are sharded
 * (many chunks in one object, read by byte range).
 *
 * A missing object is `undefined`: a public bucket that can't be listed
 * answers 403 for a missing key. A server that ignores Range and sends the
 * whole object still yields the right bytes.
 */
export async function fetchRange(
  url: string,
  range: RangeQuery,
): Promise<Uint8Array | undefined> {
  const res = await fetch(url, {
    headers: {
      Range:
        "suffixLength" in range
          ? `bytes=-${range.suffixLength}`
          : `bytes=${range.offset}-${range.offset + range.length - 1}`,
    },
  });

  if (res.status === 404 || res.status === 403) return undefined;
  if (!res.ok) throw new Error(`GET ${res.status} for ${url.split("?")[0]}`);
  const bytes = new Uint8Array(await res.arrayBuffer());

  if (res.status === 206) return bytes;

  return "suffixLength" in range
    ? bytes.slice(-range.suffixLength)
    : bytes.slice(range.offset, range.offset + range.length);
}
