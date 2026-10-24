/**
 * Decodes a SpatialData shapes element's GeoParquet (geometries as WKB) into
 * flat arrays, off the main thread: the outer ring of each polygon as x, y
 * pairs, where each polygon starts, and each polygon's id; or one point per
 * geometry, for a table whose cells have no coordinates of their own.
 */
import * as Comlink from "comlink";
import { parquetMetadataAsync, parquetReadObjects } from "hyparquet";
import { compressors } from "hyparquet-compressors";

import {
  centroid,
  centroidFromJson,
  outerRing,
  outerRingFromJson,
} from "@/lib/spatialdata/wkb";

export interface DecodedPolygons {
  /** x, y of every vertex, polygon after polygon (closing vertex dropped). */
  coords: Float32Array;
  /** Polygon i is vertices [starts[i], starts[i + 1]). */
  starts: Uint32Array;
  /** Each polygon's id (the element's index column). */
  ids: string[];
}

export interface DecodedCentroids {
  /** x, y per geometry, in file order. */
  coords: Float32Array;
  /** Each geometry's id (the element's index column). */
  ids: string[];
}

/**
 * The geometry column and the id column of a shapes GeoParquet: the
 * table's instance key when the file carries that column, else the
 * dataframe's index (when it is a column, not a stored range), else the
 * first plain column that is not the radius.
 */
async function columnsOf(buffer: ArrayBuffer, preferredId?: string) {
  const metadata = await parquetMetadataAsync(buffer);
  const kv = Object.fromEntries(
    (metadata.key_value_metadata ?? []).map((e) => [e.key, e.value]),
  );
  const geoColumn = kv.geo ? JSON.parse(kv.geo).primary_column : "geometry";
  const plain = metadata.schema.filter(
    (s, i) => i > 0 && s.name !== geoColumn && !s.num_children,
  );
  const pandasIndex = kv.pandas
    ? JSON.parse(kv.pandas).index_columns?.[0]
    : undefined;
  const idColumn: string | undefined =
    (preferredId && plain.some((s) => s.name === preferredId)
      ? preferredId
      : undefined) ??
    (typeof pandasIndex === "string" ? pandasIndex : undefined) ??
    plain.find((s) => s.name !== "radius")?.name;
  const rows = await parquetReadObjects({
    file: buffer,
    metadata,
    compressors,
    columns: idColumn ? [geoColumn, idColumn] : [geoColumn],
  });

  return { rows, geoColumn, idColumn };
}

const api = {
  async decodeCentroids(
    buffer: ArrayBuffer,
    preferredId?: string,
  ): Promise<DecodedCentroids> {
    const { rows, geoColumn, idColumn } = await columnsOf(buffer, preferredId);
    const coords = new Float32Array(rows.length * 2);
    const ids: string[] = new Array(rows.length);

    rows.forEach((row: any, i: number) => {
      const g = row[geoColumn];
      const p =
        g instanceof Uint8Array
          ? centroid(g)
          : typeof g === "object"
            ? centroidFromJson(g)
            : null;

      coords[i * 2] = p ? p[0] : NaN;
      coords[i * 2 + 1] = p ? p[1] : NaN;
      ids[i] = idColumn ? String(row[idColumn]) : String(i);
    });

    return Comlink.transfer({ coords, ids }, [coords.buffer]);
  },

  async decode(buffer: ArrayBuffer): Promise<DecodedPolygons> {
    const { rows, geoColumn, idColumn } = await columnsOf(buffer);
    const rings: Float64Array[] = [];
    const ids: string[] = [];
    let total = 0;

    rows.forEach((row: any, i: number) => {
      const g = row[geoColumn];
      let ring =
        g instanceof Uint8Array
          ? outerRing(g)
          : typeof g === "object"
            ? outerRingFromJson(g)
            : null;

      if (!ring || ring.length < 6) return;
      // Drop the closing vertex, which repeats the first
      const n = ring.length / 2;

      if (ring[0] === ring[(n - 1) * 2] && ring[1] === ring[(n - 1) * 2 + 1]) {
        ring = ring.subarray(0, (n - 1) * 2);
      }
      rings.push(ring);
      ids.push(idColumn ? String(row[idColumn]) : String(i));
      total += ring.length / 2;
    });
    const coords = new Float32Array(total * 2);
    const starts = new Uint32Array(rings.length + 1);
    let v = 0;

    rings.forEach((ring, i) => {
      starts[i] = v;
      coords.set(ring, v * 2);
      v += ring.length / 2;
    });
    starts[rings.length] = v;

    return Comlink.transfer({ coords, starts, ids }, [
      coords.buffer,
      starts.buffer,
    ]);
  },
};

export type ShapesWorkerApi = typeof api;

Comlink.expose(api);
