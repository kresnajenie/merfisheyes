/**
 * Decodes a SpatialData shapes element's GeoParquet (polygons as WKB) into
 * flat arrays, off the main thread: the outer ring of each polygon as x, y
 * pairs, where each polygon starts, and each polygon's id.
 */
import * as Comlink from "comlink";
import { parquetMetadataAsync, parquetReadObjects } from "hyparquet";
import { compressors } from "hyparquet-compressors";

export interface DecodedPolygons {
  /** x, y of every vertex, polygon after polygon (closing vertex dropped). */
  coords: Float32Array;
  /** Polygon i is vertices [starts[i], starts[i + 1]). */
  starts: Uint32Array;
  /** Each polygon's id (the element's index column). */
  ids: string[];
}

/** The outer ring of a WKB Polygon or (largest part of a) MultiPolygon. */
function outerRing(bytes: Uint8Array): Float64Array | null {
  const view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  let offset = 0;

  const header = () => {
    const little = view.getUint8(offset) === 1;

    offset += 1;
    const raw = view.getUint32(offset, little);

    offset += 4;
    // ISO (1003, 2003, 3003) and EWKB (high flag bits) say if there is Z / M
    const type = (raw & 0xffff) % 1000;
    const iso = Math.floor((raw & 0xffff) / 1000);
    const dims =
      2 +
      (iso === 1 || iso === 2 ? 1 : iso === 3 ? 2 : 0) +
      (raw & 0x80000000 ? 1 : 0) +
      (raw & 0x40000000 ? 1 : 0);

    if (raw & 0x20000000) offset += 4; // EWKB SRID

    return { little, type, dims };
  };

  const polygon = (little: boolean, dims: number) => {
    const rings = view.getUint32(offset, little);

    offset += 4;
    let outer: Float64Array | null = null;

    for (let r = 0; r < rings; r++) {
      const n = view.getUint32(offset, little);

      offset += 4;
      if (r === 0) {
        outer = new Float64Array(n * 2);
        for (let i = 0; i < n; i++) {
          outer[i * 2] = view.getFloat64(offset + i * dims * 8, little);
          outer[i * 2 + 1] = view.getFloat64(offset + i * dims * 8 + 8, little);
        }
      }
      offset += n * dims * 8;
    }

    return outer;
  };

  const { little, type, dims } = header();

  if (type === 3) return polygon(little, dims);
  if (type === 6) {
    const parts = view.getUint32(offset, little);
    let best: Float64Array | null = null;

    offset += 4;
    for (let p = 0; p < parts; p++) {
      const h = header();
      const ring = polygon(h.little, h.dims);

      if (ring && (!best || ring.length > best.length)) best = ring;
    }

    return best;
  }

  return null;
}

/** A GeoJSON geometry, as hyparquet may hand a geometry column back. */
function outerRingFromJson(g: any): Float64Array | null {
  const ring =
    g?.type === "Polygon"
      ? g.coordinates[0]
      : g?.type === "MultiPolygon"
        ? g.coordinates.reduce(
            (a: any, p: any) => (p[0].length > a.length ? p[0] : a),
            [],
          )
        : null;

  return ring
    ? Float64Array.from(ring.flatMap((c: number[]) => [c[0], c[1]]))
    : null;
}

const api = {
  async decode(buffer: ArrayBuffer): Promise<DecodedPolygons> {
    const metadata = await parquetMetadataAsync(buffer);
    const kv = Object.fromEntries(
      (metadata.key_value_metadata ?? []).map((e) => [e.key, e.value]),
    );
    const geoColumn = kv.geo ? JSON.parse(kv.geo).primary_column : "geometry";
    const idColumn =
      (kv.pandas && JSON.parse(kv.pandas).index_columns?.[0]) ||
      metadata.schema.find(
        (s, i) => i > 0 && s.name !== geoColumn && !s.num_children,
      )?.name;
    const rows = await parquetReadObjects({
      file: buffer,
      metadata,
      compressors,
      columns: idColumn ? [geoColumn, idColumn] : [geoColumn],
    });
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
