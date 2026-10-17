/**
 * Cell polygons of a SpatialData store, read from a shapes element's
 * GeoParquet as written: the file is fetched once, decoded in a worker, and
 * each polygon linked to its row of the table by the table's instance id.
 */
import type { SpatialDataStore } from "./store";
import type {
  DecodedPolygons,
  ShapesWorkerApi,
} from "@/lib/workers/shapes.worker";

import * as Comlink from "comlink";

import { elementNames } from "./store";

import { hinted } from "@/lib/storage/ElementStore";

export interface CellPolygons {
  /** x, y of every vertex, polygon after polygon. */
  coords: Float32Array;
  /** Polygon i is vertices [starts[i], starts[i + 1]). */
  starts: Uint32Array;
  /** The table row each polygon belongs to; -1 when none. */
  cell: Int32Array;
}

/** Shapes elements that hold cell outlines, best first. */
const CELL_SHAPES = ["cell_boundaries", "cell_polygons", "cells"];

/** The store's cell-outline shapes element, if it has one. */
export function cellShapesName(sd: SpatialDataStore): string | null {
  const names = elementNames(sd.keys, "shapes");

  return CELL_SHAPES.find((n) => names.includes(n)) ?? null;
}

export async function loadCellPolygons(
  sd: SpatialDataStore,
  name: string,
  instanceIds: string[],
): Promise<CellPolygons> {
  const bytes = await sd.store.get(`/shapes/${name}/shapes.parquet`, {
    signal: hinted({ tag: "shapes" }),
  });

  if (!bytes) throw new Error(`shapes/${name} has no shapes.parquet`);
  // A copy the worker can own; the store may hold on to its bytes
  const buffer = bytes.slice().buffer;
  const worker = new Worker(
    new URL("../workers/shapes.worker.ts", import.meta.url),
    { type: "module" },
  );

  try {
    const api = Comlink.wrap<ShapesWorkerApi>(worker);
    const decoded: DecodedPolygons = await api.decode(
      Comlink.transfer(buffer, [buffer]),
    );
    const row = new Map<string, number>();

    instanceIds.forEach((id, i) => row.set(id, i));
    const cell = new Int32Array(decoded.ids.length);

    decoded.ids.forEach((id, i) => {
      cell[i] = row.get(id) ?? -1;
    });

    return { coords: decoded.coords, starts: decoded.starts, cell };
  } finally {
    worker.terminate();
  }
}
