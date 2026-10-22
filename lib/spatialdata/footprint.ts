/**
 * The cells' footprint: when a table's cells are Visium spots (circles of
 * one radius) or Visium HD bins (axis-aligned squares tiling a grid), the
 * viewer draws each cell at that real size and shape instead of a dot.
 *
 * Read without the shapes themselves (22 MB for an 8 µm table, 348 MB for a
 * 2 µm one): the region element's GeoParquet footer, one suffix read, says
 * what the geometry is, and for circles holds the radius in its column
 * statistics; a bin's side is the grid's pitch, taken from the cells' own
 * coordinates (the bins tile it). Anything else — polygons, circles of
 * varying radius (Xenium's `cell_circles`) — has no footprint, and the
 * cells stay dots (polygons are drawn by the shapes layer).
 */
import type { SpatialDataStore } from "./store";

import { parquetMetadata } from "hyparquet";

import { placesAsWritten, readCoordinateSystems } from "./coordinates";
import { elementNames } from "./store";

export interface CellFootprint {
  shape: "circle" | "square";
  /** Diameter, or side, in the cells' coordinates. */
  size: number;
}

const TAIL = 65536;

/**
 * The geometry types a shapes element's GeoParquet declares (`["Point"]`,
 * `["Polygon"]`, …), from its footer alone; null when it cannot be read.
 */
export async function geometryTypes(
  sd: SpatialDataStore,
  element: string,
): Promise<string[] | null> {
  const meta = await footer(sd, `/${element}/shapes.parquet`).catch(() => null);

  if (!meta) return null;
  const geo = JSON.parse(
    meta.key_value_metadata?.find((kv) => kv.key === "geo")?.value ?? "{}",
  );

  return Object.values(
    (geo.columns ?? {}) as Record<string, { geometry_types?: string[] }>,
  ).flatMap((c) => c.geometry_types ?? []);
}

/** The footer of a parquet file, read from its end. */
async function footer(sd: SpatialDataStore, key: string) {
  let tail = await sd.store.getRange(key, { suffixLength: TAIL });

  if (!tail || tail.byteLength < 8) return null;
  const view = new DataView(tail.buffer, tail.byteOffset, tail.byteLength);
  const length = view.getUint32(tail.byteLength - 8, true);

  if (length + 8 > tail.byteLength) {
    tail = await sd.store.getRange(key, { suffixLength: length + 8 });
    if (!tail) return null;
  }

  return parquetMetadata(
    tail.buffer.slice(
      tail.byteOffset,
      tail.byteOffset + tail.byteLength,
    ) as ArrayBuffer,
  );
}

/**
 * A grid's pitch: the median distance from a cell to its nearest neighbour,
 * over a sample of up to 2,000 cells. The bins tile the grid, so this is a
 * bin's side; it holds when the grid is rotated in the image's pixels (a
 * Visium HD grid is, by a fraction of a degree, so its columns do not share
 * an x). Neighbours are found through a hash of cells about one spacing
 * wide.
 */
export function gridPitch(coords: Float32Array, dims: number): number | null {
  const n = Math.floor(coords.length / dims);

  if (n < 2) return null;
  let x0 = Infinity;
  let y0 = Infinity;
  let x1 = -Infinity;
  let y1 = -Infinity;

  for (let i = 0; i < n; i++) {
    const x = coords[i * dims];
    const y = coords[i * dims + 1];

    if (x < x0) x0 = x;
    if (x > x1) x1 = x;
    if (y < y0) y0 = y;
    if (y > y1) y1 = y;
  }
  const area = Math.max((x1 - x0) * (y1 - y0), 1e-12);
  const cell = Math.max(Math.sqrt(area / n), 1e-9);
  const key = (cx: number, cy: number) => cx * 73856093 + cy * 19349663;
  const buckets = new Map<number, number[]>();

  for (let i = 0; i < n; i++) {
    const k = key(
      Math.floor((coords[i * dims] - x0) / cell),
      Math.floor((coords[i * dims + 1] - y0) / cell),
    );
    const b = buckets.get(k);

    if (b) b.push(i);
    else buckets.set(k, [i]);
  }
  const step = Math.max(1, Math.floor(n / 2000));
  const nearest: number[] = [];

  for (let i = 0; i < n; i += step) {
    const x = coords[i * dims];
    const y = coords[i * dims + 1];
    const cx = Math.floor((x - x0) / cell);
    const cy = Math.floor((y - y0) / cell);
    let best = Infinity;

    // Rings of buckets outward until the nearest cell found cannot improve
    for (
      let r = 1;
      r <= 4 && (best === Infinity || best > (r - 1) * cell);
      r++
    ) {
      for (let dx = -r; dx <= r; dx++) {
        for (let dy = -r; dy <= r; dy++) {
          if (Math.max(Math.abs(dx), Math.abs(dy)) !== r && r > 1) continue;
          for (const j of buckets.get(key(cx + dx, cy + dy)) ?? []) {
            if (j === i) continue;
            const d = Math.hypot(
              coords[j * dims] - x,
              coords[j * dims + 1] - y,
            );

            if (d > 0 && d < best) best = d;
          }
        }
      }
    }
    if (best < Infinity) nearest.push(best);
  }
  if (!nearest.length) return null;
  nearest.sort((a, b) => a - b);

  return nearest[Math.floor(nearest.length / 2)];
}

export async function readCellFootprint(
  sd: SpatialDataStore,
  coords: Float32Array,
  dims: number,
): Promise<CellFootprint | null> {
  const cs = await readCoordinateSystems(sd);
  const region = cs.region;

  // The cells are drawn as written: their region must share their frame
  if (!region?.startsWith("shapes/") || !placesAsWritten(cs.place(region))) {
    return null;
  }
  const meta = await footer(sd, `/${region}/shapes.parquet`);

  if (!meta) return null;
  const geo = JSON.parse(
    meta.key_value_metadata?.find((kv) => kv.key === "geo")?.value ?? "{}",
  );
  const types: string[] = Object.values(
    (geo.columns ?? {}) as Record<string, { geometry_types?: string[] }>,
  ).flatMap((c) => c.geometry_types ?? []);
  if (types.length === 1 && types[0] === "Point") {
    const radius = oneRadius(meta);

    return radius ? { shape: "circle", size: 2 * radius } : null;
  }
  if (types.length === 1 && types[0] === "Polygon") {
    if (/square/i.test(region)) {
      const side = gridPitch(coords, dims);

      return side ? { shape: "square", size: side } : null;
    }
    // Polygons of one size: SpatialData-DB's `_for_vitessce` sibling gives
    // every cell the same radius when they are a grid's pixels (the squares
    // of a DBiT CUT&Tag or ATAC capture, side twice the radius); cells of
    // varying size (CODEX) have no footprint
    const sibling = `${region}_for_vitessce`;

    if (
      !elementNames(sd.keys, "shapes").includes(sibling.slice("shapes/".length))
    ) {
      return null;
    }
    const siblingMeta = await footer(sd, `/${sibling}/shapes.parquet`);
    const radius = siblingMeta && oneRadius(siblingMeta);

    return radius ? { shape: "square", size: 2 * radius } : null;
  }

  return null;
}

/** The `radius` column's one value, when every row has the same; else null. */
function oneRadius(meta: ReturnType<typeof parquetMetadata>): number | null {
  if (!meta.schema.some((s) => s.name === "radius")) return null;
  const index = meta.row_groups[0]?.columns.findIndex(
    (c) => c.meta_data?.path_in_schema.join(".") === "radius",
  );
  let min = Infinity;
  let max = -Infinity;

  for (const group of meta.row_groups) {
    const stats = group.columns[index]?.meta_data?.statistics;

    if (stats?.min_value == null || stats.max_value == null) return null;
    min = Math.min(min, Number(stats.min_value));
    max = Math.max(max, Number(stats.max_value));
  }

  return max > 0 && max - min <= max * 1e-6 ? max : null;
}
