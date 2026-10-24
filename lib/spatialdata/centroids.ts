/**
 * Cell coordinates from the shapes element the table annotates: one point
 * per cell, read from its GeoParquet as written and matched to the table's
 * rows by instance id.
 *
 * In SpatialData the region element is the cells' frame; `obsm/spatial` is
 * a copy the format does not tie to it, and in SpatialData-DB's MALDI
 * stores it is the ion image's pixel index (the spots are at 25 or 50 times
 * it). So when the region has a Point element — itself, or the
 * `<region>_for_vitessce` sibling SpatialData-DB writes (0.1–3.3 MB, 35 MB
 * for 1.46 M CosMx cells) — the cells are placed from it. A table with no
 * `obsm/spatial` and only polygons takes each polygon's vertex mean.
 */
import type { SpatialDataTableAdapter } from "./SpatialDataTableAdapter";
import type { SpatialDataStore } from "./store";
import type { ShapesWorkerApi } from "@/lib/workers/shapes.worker";

import * as Comlink from "comlink";

import { placesAsWritten, readCoordinateSystems } from "./coordinates";
import { pointSiblingOf } from "./conventions";
import { geometryTypes } from "./footprint";
import { elementNames } from "./store";

import { hinted } from "@/lib/storage/ElementStore";

/** The shapes element to take the cells' positions from, if any. */
export function centroidSource(sd: SpatialDataStore, region: string | null) {
  if (!region?.startsWith("shapes/")) return null;
  const name = region.slice("shapes/".length);
  const names = elementNames(sd.keys, "shapes");
  const points = pointSiblingOf(name);

  return names.includes(points) ? `shapes/${points}` : region;
}

/**
 * The Point element holding the cells' centres, when the region has one
 * (the region itself when its geometry is Point, else its `_for_vitessce`
 * sibling); null otherwise.
 */
export async function pointSource(
  sd: SpatialDataStore,
): Promise<string | null> {
  const cs = await readCoordinateSystems(sd);
  const source = centroidSource(sd, cs.region);

  if (!source || !placesAsWritten(cs.place(cs.region!))) return null;
  if (source !== cs.region) return source;
  const types = await geometryTypes(sd, source);

  return types && types.length === 1 && types[0] === "Point" ? source : null;
}

export async function readCellCentroids(
  sd: SpatialDataStore,
  adapter: SpatialDataTableAdapter,
  from?: string,
): Promise<{ coordinates: Float32Array; dimensions: number }> {
  const cs = await readCoordinateSystems(sd);
  const source = from ?? centroidSource(sd, cs.region);

  if (!source) {
    throw new Error(
      "No spatial coordinates: the table has no obsm/spatial and annotates no shapes element",
    );
  }
  // The cells are drawn as written, so the shapes must share their frame
  if (!placesAsWritten(cs.place(cs.region!))) {
    const p = cs.place(cs.region!);

    throw new Error(
      `No spatial coordinates: ${cs.region} is not in the cells' frame${p.ok ? "" : ` (${p.reason})`}`,
    );
  }
  const bytes = await sd.store.get(`/${source}/shapes.parquet`, {
    signal: hinted({ tag: "shapes" }),
  });

  if (!bytes) throw new Error(`${source} has no shapes.parquet`);
  const buffer = bytes.slice().buffer;
  const worker = new Worker(
    new URL("../workers/shapes.worker.ts", import.meta.url),
    { type: "module" },
  );

  try {
    const api = Comlink.wrap<ShapesWorkerApi>(worker);
    const decoded = await api.decodeCentroids(
      Comlink.transfer(buffer, [buffer]),
      adapter.instanceKey ?? undefined,
    );
    const ids = await adapter.readInstanceIds();
    const n = adapter.getDatasetInfo().numCells;
    const coordinates = new Float32Array(n * 2);
    let matched = 0;

    if (ids && ids.length === n) {
      const row = new Map<string, number>();

      ids.forEach((id, i) => row.set(id, i));
      decoded.ids.forEach((id, i) => {
        const r = row.get(id);

        if (r === undefined) return;
        coordinates[r * 2] = decoded.coords[i * 2];
        coordinates[r * 2 + 1] = decoded.coords[i * 2 + 1];
        matched++;
      });
    } else if (decoded.ids.length === n) {
      // No instance ids to match by: the shapes are in table order
      coordinates.set(decoded.coords.subarray(0, n * 2));
      matched = n;
    }
    if (matched < n * 0.99) {
      throw new Error(
        `No spatial coordinates: ${source} places ${matched} of ${n} cells`,
      );
    }
    if (matched < n) {
      console.warn(
        `[centroids] ${n - matched} cells have no shape; drawn at the origin`,
      );
    }
    adapter.spatialDimensions = 2;

    return { coordinates, dimensions: 2 };
  } finally {
    worker.terminate();
  }
}
