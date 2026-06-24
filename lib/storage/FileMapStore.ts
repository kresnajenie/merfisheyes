import type { AbsolutePath, AsyncReadable } from "@zarrita/storage";

/**
 * Browser File-backed zarrita store.
 *
 * Wraps a `Map<relativePath, File>` (e.g. from a `webkitdirectory` drop)
 * and exposes the zarrita `AsyncReadable` interface so it can be passed
 * to `zarr.open()` / `readZarr()`.
 *
 * Keys arrive as absolute zarr paths like `/X/data/0.0`; we strip the
 * leading slash and look the file up in the map.
 */
export class FileMapStore implements AsyncReadable {
  private files: Map<string, File>;

  constructor(files: Map<string, File>) {
    this.files = files;
  }

  async get(key: AbsolutePath): Promise<Uint8Array | undefined> {
    const relativeKey = key.startsWith("/") ? key.slice(1) : key;
    const file = this.files.get(relativeKey);

    if (!file) return undefined;

    const buffer = await file.arrayBuffer();

    return new Uint8Array(buffer);
  }
}
