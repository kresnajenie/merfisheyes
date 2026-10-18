import { ElementStore } from "@/lib/storage/ElementStore";
import { fileBackend } from "@/lib/storage/backends";

/**
 * Read-only store over files picked or dropped in the browser. Nothing is
 * read until a key is asked for, and a ranged read (sharded arrays) only
 * touches that slice of its file — so a store of any size opens at once.
 *
 * `files` is keyed by path from the store's root, without a leading slash.
 */
export function fileStore(files: Map<string, File>): ElementStore {
  return new ElementStore(fileBackend(files));
}

export { openSpatialDataStoreFromFiles } from "./store";
