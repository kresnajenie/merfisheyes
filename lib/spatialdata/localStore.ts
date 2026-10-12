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
 * A SpatialData store from dropped files. `files` is the store, keyed from
 * its root; `beside` is whatever was dropped around it (keyed from the
 * dropped folder), where its sidecars — annotations/ — are looked for.
 */
export function openSpatialDataStoreFromFiles(
  files: Map<string, File>,
  beside: Map<string, File> = new Map(),
): SpatialDataStore {
  return {
    store: fileStore(files),
    keys: Array.from(files.keys()),
    sidecar: fileStore(beside),
  };
}
