/**
 * One-shot worker that densifies a CSR `X` matrix from an h5ad-zarr store.
 *
 * Caller hands us the same `Map<string, File>` the main thread used to open
 * the zarr, plus the matrix shape. We re-open the store inside the worker,
 * read `X/indices`, `X/indptr`, `X/data`, and build a row-major
 * `Float32Array(numCells * numGenes)`. We post the buffer back as a
 * Transferable so there is no clone cost across the worker boundary.
 */
import * as Comlink from "comlink";
import * as zarr from "zarrita";

import { FileMapStore } from "../storage/FileMapStore";

export type DensifyProgress = (progress: number, message: string) => void;

export interface DensifyResult {
  buffer: ArrayBuffer; // backing store of the dense Float32Array
  numCells: number;
  numGenes: number;
}

const api = {
  async densifyCSR(
    fileMap: Map<string, File>,
    numCells: number,
    numGenes: number,
    xPath: string,
    onProgress?: DensifyProgress,
  ): Promise<DensifyResult> {
    onProgress?.(2, "Opening zarr store in worker...");
    const store = new FileMapStore(fileMap);
    const root = await zarr.open(store, { kind: "group" });

    const xGroup = await zarr.open(root.resolve(xPath), { kind: "group" });

    onProgress?.(8, "Loading CSR indptr...");
    const indptrArr = await zarr.open(xGroup.resolve("indptr"), {
      kind: "array",
    });
    const indptrChunk = await zarr.get(indptrArr, [null]);
    const indptr = indptrChunk.data as ArrayLike<number | bigint>;

    onProgress?.(20, "Loading CSR indices...");
    const indicesArr = await zarr.open(xGroup.resolve("indices"), {
      kind: "array",
    });
    const indicesChunk = await zarr.get(indicesArr, [null]);
    const indices = indicesChunk.data as ArrayLike<number | bigint>;

    onProgress?.(45, "Loading CSR data...");
    const dataArr = await zarr.open(xGroup.resolve("data"), { kind: "array" });
    const dataChunk = await zarr.get(dataArr, [null]);
    const data = dataChunk.data as ArrayLike<number | bigint>;

    onProgress?.(65, `Allocating ${numCells.toLocaleString()} × ${numGenes.toLocaleString()} dense matrix...`);
    const denseLen = numCells * numGenes;
    const dense = new Float32Array(denseLen);

    onProgress?.(70, "Densifying CSR matrix...");
    const reportEvery = Math.max(1, Math.floor(numCells / 20));

    for (let r = 0; r < numCells; r++) {
      const start = Number(indptr[r]);
      const stop = Number(indptr[r + 1]);
      const rowOffset = r * numGenes;

      for (let k = start; k < stop; k++) {
        const col = Number(indices[k]);

        dense[rowOffset + col] = Number(data[k]);
      }
      if (r % reportEvery === 0) {
        const pct = 70 + Math.floor((r / numCells) * 28);

        onProgress?.(pct, `Densifying row ${r.toLocaleString()} / ${numCells.toLocaleString()}`);
      }
    }

    onProgress?.(99, "Transferring matrix to main thread...");

    return Comlink.transfer(
      { buffer: dense.buffer, numCells, numGenes },
      [dense.buffer],
    );
  },
};

export type ZarrDensifyWorkerApi = typeof api;

Comlink.expose(api);
