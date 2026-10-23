/**
 * The little of WKB (well-known binary) a shapes element needs: the outer
 * ring of a polygon, and one point per geometry (a Point itself, else the
 * mean of the outer ring's vertices). Pure functions, shared by the shapes
 * worker and the tests.
 */

interface Header {
  little: boolean;
  type: number;
  dims: number;
}

function header(view: DataView, at: number): { h: Header; next: number } {
  const little = view.getUint8(at) === 1;
  const raw = view.getUint32(at + 1, little);
  // ISO (1003, 2003, 3003) and EWKB (high flag bits) say if there is Z / M
  const type = (raw & 0xffff) % 1000;
  const iso = Math.floor((raw & 0xffff) / 1000);
  const dims =
    2 +
    (iso === 1 || iso === 2 ? 1 : iso === 3 ? 2 : 0) +
    (raw & 0x80000000 ? 1 : 0) +
    (raw & 0x40000000 ? 1 : 0);
  const srid = raw & 0x20000000 ? 4 : 0; // EWKB SRID

  return { h: { little, type, dims }, next: at + 5 + srid };
}

/** The outer ring of a Polygon at `at`; returns the offset after the polygon. */
function polygon(
  view: DataView,
  at: number,
  h: Header,
): { outer: Float64Array | null; next: number } {
  const rings = view.getUint32(at, h.little);
  let offset = at + 4;
  let outer: Float64Array | null = null;

  for (let r = 0; r < rings; r++) {
    const n = view.getUint32(offset, h.little);

    offset += 4;
    if (r === 0) {
      outer = new Float64Array(n * 2);
      for (let i = 0; i < n; i++) {
        outer[i * 2] = view.getFloat64(offset + i * h.dims * 8, h.little);
        outer[i * 2 + 1] = view.getFloat64(
          offset + i * h.dims * 8 + 8,
          h.little,
        );
      }
    }
    offset += n * h.dims * 8;
  }

  return { outer, next: offset };
}

/** The outer ring of a WKB Polygon or (largest part of a) MultiPolygon. */
export function outerRing(bytes: Uint8Array): Float64Array | null {
  const view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  const { h, next } = header(view, 0);

  if (h.type === 3) return polygon(view, next, h).outer;
  if (h.type === 6) {
    const parts = view.getUint32(next, h.little);
    let offset = next + 4;
    let best: Float64Array | null = null;

    for (let p = 0; p < parts; p++) {
      const part = header(view, offset);
      const { outer, next: after } = polygon(view, part.next, part.h);

      offset = after;
      if (outer && (!best || outer.length > best.length)) best = outer;
    }

    return best;
  }

  return null;
}

/**
 * One point for a geometry: a Point's own coordinates; a polygon's vertex
 * mean (the closing vertex counted once). Null for anything else.
 */
export function centroid(bytes: Uint8Array): [number, number] | null {
  const view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  const { h, next } = header(view, 0);

  if (h.type === 1) {
    return [
      view.getFloat64(next, h.little),
      view.getFloat64(next + 8, h.little),
    ];
  }
  const ring = outerRing(bytes);

  if (!ring || ring.length < 4) return null;
  let n = ring.length / 2;

  if (
    n > 1 &&
    ring[0] === ring[(n - 1) * 2] &&
    ring[1] === ring[(n - 1) * 2 + 1]
  ) {
    n--;
  }
  let x = 0;
  let y = 0;

  for (let i = 0; i < n; i++) {
    x += ring[i * 2];
    y += ring[i * 2 + 1];
  }

  return [x / n, y / n];
}

/** The same for a GeoJSON geometry, as hyparquet may hand a column back. */
export function centroidFromJson(g: any): [number, number] | null {
  if (g?.type === "Point") return [g.coordinates[0], g.coordinates[1]];
  const ring = outerRingFromJson(g);

  if (!ring || ring.length < 4) return null;
  const n = ring.length / 2;
  let x = 0;
  let y = 0;

  for (let i = 0; i < n; i++) {
    x += ring[i * 2];
    y += ring[i * 2 + 1];
  }

  return [x / n, y / n];
}

/** A GeoJSON Polygon's outer ring (the largest part of a MultiPolygon). */
export function outerRingFromJson(g: any): Float64Array | null {
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
