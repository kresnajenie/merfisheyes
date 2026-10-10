import {
  fetchDatasetKeyList,
  PresignedFetchStore,
} from "@/lib/storage/PresignedFetchStore";

export interface SpatialDataStore {
  store: { get(key: string): Promise<Uint8Array | undefined> };
  /**
   * Keys in the store, relative to its root. At least the `zarr.json` of
   * every array and group — enough to enumerate elements and columns.
   */
  keys: string[];
  /**
   * The folder that holds the store, where files that describe it without
   * being part of it live (annotations/, mapping.json).
   */
  sidecar: { get(key: string): Promise<Uint8Array | undefined> };
}

/**
 * Open a SpatialData zarr store on our S3. The store may sit anywhere under
 * `datasets/{id}/` — it is found by its root `zarr.json`.
 */
export async function openSpatialDataStore(
  datasetId: string,
): Promise<SpatialDataStore> {
  const allKeys = await fetchDatasetKeyList(datasetId);
  const rootMarker = allKeys
    .filter((k) => k === "zarr.json" || k.endsWith("/zarr.json"))
    .sort((a, b) => a.length - b.length)[0];

  if (rootMarker === undefined) {
    throw new Error("No zarr v3 store found for this dataset in S3");
  }
  const prefix = rootMarker.slice(0, -"zarr.json".length);

  // "data.zarr/" → "", "a/b.zarr/" → "a/"
  const folder = prefix.replace(/[^/]+\/$/, "");

  return {
    store: new PresignedFetchStore(datasetId, { keyPrefix: prefix }),
    sidecar: new PresignedFetchStore(datasetId, { keyPrefix: folder }),
    keys: allKeys
      .filter((k) => k.startsWith(prefix))
      .map((k) => k.slice(prefix.length)),
  };
}

/**
 * Open a SpatialData zarr store at a public URL (e.g. a public S3 bucket).
 *
 * A plain URL can't be listed, so the store's contents come from the
 * consolidated metadata in its root `zarr.json`.
 */
export async function openSpatialDataStoreFromUrl(
  url: string,
): Promise<SpatialDataStore> {
  const base = url.replace(/\/+$/, "");
  const fetchFrom = (root: string) => ({
    async get(key: string) {
      const res = await fetch(root + key);

      // A public bucket that can't be listed answers 403 for a missing key
      if (res.status === 404 || res.status === 403) return undefined;
      if (!res.ok) throw new Error(`GET ${res.status} for ${key}`);

      return new Uint8Array(await res.arrayBuffer());
    },
  });
  const store = fetchFrom(base);
  const rootBytes = await store.get("/zarr.json");

  if (!rootBytes) throw new Error(`No zarr v3 store found at ${base}`);
  const nodes = JSON.parse(new TextDecoder().decode(rootBytes))
    .consolidated_metadata?.metadata;

  if (!nodes) {
    throw new Error("Store has no consolidated metadata to list it by");
  }

  return {
    store,
    sidecar: fetchFrom(base.slice(0, base.lastIndexOf("/"))),
    keys: ["zarr.json", ...Object.keys(nodes).map((p) => `${p}/zarr.json`)],
  };
}

/** Names of the elements of one type (`images`, `shapes`, `points`, ...). */
export function elementNames(keys: string[], type: string): string[] {
  const names = new Set<string>();

  for (const key of keys) {
    const parts = key.split("/");

    if (parts[0] === type && parts.length > 2) names.add(parts[1]);
  }

  return Array.from(names).sort();
}
