import type { SpatialDataStore } from "./store";
import type { RangeQuery } from "@zarrita/storage";

/**
 * Read-only zarr store over files picked or dropped in the browser. Nothing
 * is read until a key is asked for, and a ranged read (sharded arrays) only
 * touches that slice of its file — so a store of any size opens at once.
 *
 * `files` is keyed by path from the store's root, without a leading slash.
 */
export function fileStore(files: Map<string, File>) {
  const find = (key: string) => files.get(key.replace(/^\/+/, ""));

  return {
    async get(key: string) {
      const file = find(key);

      return file && new Uint8Array(await file.arrayBuffer());
    },
    async getRange(key: string, range: RangeQuery) {
      const file = find(key);

      if (!file) return undefined;
      const part =
        "suffixLength" in range
          ? file.slice(-range.suffixLength)
          : file.slice(range.offset, range.offset + range.length);

      return new Uint8Array(await part.arrayBuffer());
    },
  };
}

/**
 * A SpatialData store among dropped files. `folder` holds everything that was
 * dropped (keyed from the dropped folder); the store is its subfolder `root`
 * ("" when the dropped folder is the store itself), and the rest is where its
 * sidecars — annotations/ — are looked for.
 */
export function openSpatialDataStoreFromFiles(
  folder: Map<string, File>,
  root: string,
): SpatialDataStore {
  const prefix = root ? `${root}/` : "";
  const inStore = new Map<string, File>();

  for (const [key, file] of folder) {
    if (key.startsWith(prefix)) inStore.set(key.slice(prefix.length), file);
  }

  return {
    store: fileStore(inStore),
    keys: Array.from(inStore.keys()),
    sidecar: fileStore(root ? folder : new Map()),
  };
}
