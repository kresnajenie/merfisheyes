import { ElementStore } from "@/lib/storage/ElementStore";
import {
  datasetBackend,
  fileBackend,
  urlBackend,
} from "@/lib/storage/backends";
import { fetchDatasetKeyList } from "@/lib/storage/PresignedFetchStore";

/**
 * A SpatialData Zarr v3 store, open: one ElementStore for its objects (the
 * cache, the queue, cancellation — docs/specs/LOADING.md), its node list,
 * and the folder beside it where files that describe it without being part
 * of it live (annotations/, mapping.json).
 */
export interface SpatialDataStore {
  store: ElementStore;
  /**
   * Keys in the store, relative to its root. At least the `zarr.json` of
   * every array and group — enough to enumerate elements and columns.
   */
  keys: string[];
  sidecar: ElementStore;
  /** The root `zarr.json`, parsed. */
  root: Record<string, any>;
}

const decode = (b: Uint8Array) => JSON.parse(new TextDecoder().decode(b));

/**
 * Finish opening: the root `zarr.json` is read once, its consolidated
 * metadata answers every node's `zarr.json` from then on and names the
 * nodes; a store without it is listed by `listAll` (an app dataset), or
 * cannot be opened (a public URL lists nothing).
 */
async function finishOpen(
  store: ElementStore,
  sidecar: ElementStore,
  listAll: () => Promise<string[]>,
  where: string,
): Promise<SpatialDataStore> {
  const rootBytes = await store.get("/zarr.json");

  if (!rootBytes) throw new Error(`No zarr v3 store found ${where}`);
  let root: Record<string, any> = {};

  try {
    root = decode(rootBytes);
  } catch {
    // Not JSON: listed below, like a store without consolidated metadata
  }
  const nodes = store.seedConsolidated(root);
  const keys = nodes
    ? ["zarr.json", ...nodes.map((p) => `${p}/zarr.json`)]
    : await listAll();

  return { store, keys, sidecar, root };
}

/**
 * Where a Zarr v3 store sits in an app dataset: at the dataset root, or in
 * one folder under it (`data.zarr/`). One or two requests; the whole key
 * list only when neither holds it.
 */
export async function findZarrRoot(
  datasetId: string,
): Promise<{ prefix: string; store: ElementStore; rootBytes: Uint8Array }> {
  const probe = async (prefix: string) => {
    const store = new ElementStore(datasetBackend(datasetId, prefix));
    const rootBytes = await store.get("/zarr.json");

    return rootBytes && { prefix, store, rootBytes };
  };
  const atRoot = await probe("");

  if (atRoot) return atRoot;
  const res = await fetch(`/api/datasets/${datasetId}/list?delimiter=/`);
  const { prefixes = [] } = res.ok
    ? ((await res.json()) as { prefixes?: string[] })
    : {};

  for (const prefix of prefixes.sort((a, b) => a.length - b.length)) {
    const found = await probe(prefix);

    if (found) return found;
  }
  // Deeper than one folder: the shortest zarr.json in the whole listing
  const marker = (await fetchDatasetKeyList(datasetId))
    .filter((k) => k === "zarr.json" || k.endsWith("/zarr.json"))
    .sort((a, b) => a.length - b.length)[0];
  const found = marker && (await probe(marker.slice(0, -"zarr.json".length)));

  if (!found) throw new Error("No zarr v3 store found for this dataset in S3");

  return found;
}

/**
 * Open a SpatialData zarr store on our S3. The store may sit anywhere under
 * `datasets/{id}/` — it is found by its root `zarr.json`.
 */
export async function openSpatialDataStore(
  datasetId: string,
): Promise<SpatialDataStore> {
  const { prefix, store, rootBytes } = await findZarrRoot(datasetId);

  store.seed("/zarr.json", rootBytes);
  // "data.zarr/" → "", "a/b.zarr/" → "a/"
  const folder = prefix.replace(/[^/]+\/$/, "");

  return finishOpen(
    store,
    new ElementStore(datasetBackend(datasetId, folder)),
    async () =>
      (await fetchDatasetKeyList(datasetId))
        .filter((k) => k.startsWith(prefix))
        .map((k) => k.slice(prefix.length)),
    "for this dataset in S3",
  );
}

/** Read-only store over a public URL; keys are paths from `root`. */
export function publicUrlStore(root: string): ElementStore {
  return new ElementStore(urlBackend(root));
}

/**
 * Open a SpatialData zarr store at a public URL (e.g. a public S3 bucket).
 * A plain URL can't be listed, so the store's contents come from the
 * consolidated metadata in its root `zarr.json` — given here when it was
 * already read to find the store (`resolveZarrUrl`).
 */
export async function openSpatialDataStoreFromUrl(
  url: string,
  root?: Record<string, any>,
): Promise<SpatialDataStore> {
  const base = url.replace(/\/+$/, "");
  const store = publicUrlStore(base);

  if (root)
    store.seed("/zarr.json", new TextEncoder().encode(JSON.stringify(root)));

  return finishOpen(
    store,
    publicUrlStore(base.slice(0, base.lastIndexOf("/"))),
    async () => {
      throw new Error("Store has no consolidated metadata to list it by");
    },
    `at ${base}`,
  );
}

/**
 * A SpatialData store from dropped files. `files` is the store, keyed from
 * its root; `beside` is whatever was dropped around it (keyed from the
 * dropped folder), where its sidecars — annotations/ — are looked for.
 * Nothing is read until a key is asked for.
 */
export async function openSpatialDataStoreFromFiles(
  files: Map<string, File>,
  beside: Map<string, File> = new Map(),
): Promise<SpatialDataStore> {
  return finishOpen(
    new ElementStore(fileBackend(files)),
    new ElementStore(fileBackend(beside)),
    async () => Array.from(files.keys()),
    "in these files",
  );
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
