import type { AsyncBuffer } from "hyparquet";

/**
 * hyparquet's package types are its browser build; in Node the same import
 * resolves to `src/node.js`, which also exports a file-backed buffer. Only
 * Node-side code (the catalog sync) may call it.
 */
declare module "hyparquet" {
  export function asyncBufferFromFile(filename: string): Promise<AsyncBuffer>;
}
