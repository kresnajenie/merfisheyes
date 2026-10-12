/**
 * What a public URL points at, when it is not a merfisheyes folder: a
 * SpatialData store, an OME-Zarr image, or a folder holding a store.
 *
 * A plain URL can't be listed, so a folder's store is found by name: the
 * `cells` entry of the mapping.json in the folder (a name inside it, or a
 * full URL), else a store named after the folder (`sample/sample.zarr`).
 */
export interface ZarrTarget {
  kind: "spatialdata" | "image";
  /** The store or image itself, without a trailing slash. */
  url: string;
}

async function getJson(url: string): Promise<any | null> {
  try {
    const res = await fetch(url);

    return res.ok ? await res.json() : null;
  } catch {
    return null;
  }
}

async function zarrKind(url: string): Promise<ZarrTarget["kind"] | null> {
  const attrs = (await getJson(`${url}/zarr.json`))?.attributes;

  if (attrs?.spatialdata_attrs) return "spatialdata";
  if (attrs?.ome?.multiscales) return "image";

  return null;
}

export async function resolveZarrUrl(
  input: string,
): Promise<ZarrTarget | null> {
  const base = input.trim().replace(/\/+$/, "");
  const direct = await zarrKind(base);

  if (direct) return { kind: direct, url: base };

  const cells = (await getJson(`${base}/mapping.json`))?.cells;
  const candidates = [`${base}/${base.slice(base.lastIndexOf("/") + 1)}.zarr`];

  if (typeof cells === "string" && cells) {
    candidates.unshift(
      /^https?:\/\//.test(cells)
        ? cells.replace(/\/+$/, "")
        : `${base}/${cells.replace(/^\/+|\/+$/g, "")}`,
    );
  }
  for (const url of candidates) {
    if ((await zarrKind(url)) === "spatialdata") {
      return { kind: "spatialdata", url };
    }
  }

  return null;
}
