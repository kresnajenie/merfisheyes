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
  /** Its root `zarr.json`, as read to tell what it is. */
  root: Record<string, any>;
}

async function getJson(url: string): Promise<any | null> {
  try {
    const res = await fetch(url);

    return res.ok ? await res.json() : null;
  } catch {
    return null;
  }
}

async function zarrKind(url: string): Promise<Omit<ZarrTarget, "url"> | null> {
  const root = await getJson(`${url}/zarr.json`);
  const attrs = root?.attributes;

  if (attrs?.spatialdata_attrs) return { kind: "spatialdata", root };
  if (attrs?.ome?.multiscales) return { kind: "image", root };

  return null;
}

export async function resolveZarrUrl(
  input: string,
): Promise<ZarrTarget | null> {
  const base = input.trim().replace(/\/+$/, "");
  const direct = await zarrKind(base);

  if (direct) return { ...direct, url: base };

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
    const found = await zarrKind(url);

    if (found?.kind === "spatialdata") return { ...found, url };
  }

  return null;
}
