/**
 * Per-sample 2D transforms (translate + rotate around sample centroid).
 *
 * All distances are in the dataset's raw coordinate units (whatever
 * `spatial.coordinates` stores). Rotation is in radians internally.
 */

export interface SampleTransform {
  dx: number;
  dy: number;
  theta: number; // radians
}

export const IDENTITY_TRANSFORM: SampleTransform = { dx: 0, dy: 0, theta: 0 };

export function isIdentity(t: SampleTransform): boolean {
  return t.dx === 0 && t.dy === 0 && t.theta === 0;
}

interface SpatialLike {
  coordinates: Float32Array | number[][];
  dimensions: number;
}

interface ClusterLike {
  values: any[];
  valueIndices?: Uint16Array | Uint32Array;
  uniqueValues?: string[];
}

/** Read sample value for cell i, preferring fast valueIndices lookup. */
function getSampleId(cluster: ClusterLike, i: number): string {
  if (cluster.valueIndices && cluster.uniqueValues) {
    return cluster.uniqueValues[cluster.valueIndices[i]];
  }
  return String(cluster.values[i]);
}

/** Mean (x, y) per sample value, in raw coord units. */
export function computeSampleCentroids(
  spatial: SpatialLike,
  cluster: ClusterLike,
): Map<string, [number, number]> {
  const dims = spatial.dimensions;
  const coords = spatial.coordinates;
  const n =
    coords instanceof Float32Array
      ? coords.length / dims
      : (coords as number[][]).length;

  const sumX = new Map<string, number>();
  const sumY = new Map<string, number>();
  const counts = new Map<string, number>();

  const getX = coords instanceof Float32Array
    ? (i: number) => (coords as Float32Array)[i * dims]
    : (i: number) => (coords as number[][])[i][0];
  const getY = coords instanceof Float32Array
    ? (i: number) => (coords as Float32Array)[i * dims + 1]
    : (i: number) => (coords as number[][])[i][1];

  for (let i = 0; i < n; i++) {
    const s = getSampleId(cluster, i);
    sumX.set(s, (sumX.get(s) ?? 0) + getX(i));
    sumY.set(s, (sumY.get(s) ?? 0) + getY(i));
    counts.set(s, (counts.get(s) ?? 0) + 1);
  }

  const centroids = new Map<string, [number, number]>();
  for (const [s, c] of counts) {
    centroids.set(s, [(sumX.get(s) ?? 0) / c, (sumY.get(s) ?? 0) / c]);
  }
  return centroids;
}

/**
 * Apply per-sample transforms to a Float32Array of scene-space positions
 * (xyz interleaved, length = n*3). `basePositions` holds the untransformed
 * scene-space coords; `out` is written in place. Z values pass through.
 *
 * Centroids and dx/dy are in *raw* coord units — they are multiplied by
 * `scale` here so we can stay in scene space throughout.
 */
export function applyTransformsToPositions(
  basePositions: Float32Array,
  out: Float32Array,
  cluster: ClusterLike | null,
  centroids: Map<string, [number, number]> | null,
  transforms: Map<string, SampleTransform>,
  scale: number,
): void {
  // Fast path: no transforms → just copy
  if (!cluster || !centroids || transforms.size === 0) {
    out.set(basePositions);
    return;
  }

  // Pre-resolve transform + centroid per unique value to avoid per-cell map
  // lookups (huge win for 700k+ points).
  const uniqueValues = cluster.uniqueValues;
  const valueIndices = cluster.valueIndices;
  const n = basePositions.length / 3;

  if (uniqueValues && valueIndices) {
    const cosT = new Float32Array(uniqueValues.length);
    const sinT = new Float32Array(uniqueValues.length);
    const cxS = new Float32Array(uniqueValues.length);
    const cyS = new Float32Array(uniqueValues.length);
    const dxS = new Float32Array(uniqueValues.length);
    const dyS = new Float32Array(uniqueValues.length);
    const hasT = new Uint8Array(uniqueValues.length);

    for (let k = 0; k < uniqueValues.length; k++) {
      const t = transforms.get(uniqueValues[k]);
      if (!t || isIdentity(t)) continue;
      const c = centroids.get(uniqueValues[k]);
      if (!c) continue;
      cosT[k] = Math.cos(t.theta);
      sinT[k] = Math.sin(t.theta);
      cxS[k] = c[0] * scale;
      cyS[k] = c[1] * scale;
      dxS[k] = t.dx * scale;
      dyS[k] = t.dy * scale;
      hasT[k] = 1;
    }

    for (let i = 0; i < n; i++) {
      const k = valueIndices[i];
      const baseX = basePositions[i * 3];
      const baseY = basePositions[i * 3 + 1];
      if (hasT[k]) {
        const ox = baseX - cxS[k];
        const oy = baseY - cyS[k];
        out[i * 3] = ox * cosT[k] - oy * sinT[k] + cxS[k] + dxS[k];
        out[i * 3 + 1] = ox * sinT[k] + oy * cosT[k] + cyS[k] + dyS[k];
      } else {
        out[i * 3] = baseX;
        out[i * 3 + 1] = baseY;
      }
      out[i * 3 + 2] = basePositions[i * 3 + 2];
    }
    return;
  }

  // Fallback: read sample id per-cell, look up transform/centroid in maps.
  for (let i = 0; i < n; i++) {
    const s = getSampleId(cluster, i);
    const t = transforms.get(s);
    const c = centroids.get(s);
    const baseX = basePositions[i * 3];
    const baseY = basePositions[i * 3 + 1];
    if (t && c && !isIdentity(t)) {
      const cx = c[0] * scale;
      const cy = c[1] * scale;
      const cosT = Math.cos(t.theta);
      const sinT = Math.sin(t.theta);
      const ox = baseX - cx;
      const oy = baseY - cy;
      out[i * 3] = ox * cosT - oy * sinT + cx + t.dx * scale;
      out[i * 3 + 1] = ox * sinT + oy * cosT + cy + t.dy * scale;
    } else {
      out[i * 3] = baseX;
      out[i * 3 + 1] = baseY;
    }
    out[i * 3 + 2] = basePositions[i * 3 + 2];
  }
}

/**
 * Build a gzipped meyes-format spatial.bin.gz Blob with transforms baked in.
 * Output matches process_spatial_data.py write_coordinate_binary():
 *   header: <num_points: u32, dimensions: u32> (little-endian)
 *   body:   float32[num_points * dimensions], row-major
 */
export async function buildSpatialBinGz(
  spatial: SpatialLike,
  cluster: ClusterLike | null,
  centroids: Map<string, [number, number]> | null,
  transforms: Map<string, SampleTransform>,
): Promise<Blob> {
  const dims = spatial.dimensions;
  const coords = spatial.coordinates;
  const n =
    coords instanceof Float32Array
      ? coords.length / dims
      : (coords as number[][]).length;

  // Build a fresh Float32Array of transformed raw coords.
  const out = new Float32Array(n * dims);
  const getX = coords instanceof Float32Array
    ? (i: number) => (coords as Float32Array)[i * dims]
    : (i: number) => (coords as number[][])[i][0];
  const getY = coords instanceof Float32Array
    ? (i: number) => (coords as Float32Array)[i * dims + 1]
    : (i: number) => (coords as number[][])[i][1];
  const getZ =
    dims === 3
      ? coords instanceof Float32Array
        ? (i: number) => (coords as Float32Array)[i * dims + 2]
        : (i: number) => (coords as number[][])[i][2]
      : null;

  const haveTransforms = !!(cluster && centroids && transforms.size > 0);
  const uniqueValues = cluster?.uniqueValues;
  const valueIndices = cluster?.valueIndices;

  // Per-unique-value resolved transform (raw units, not scaled).
  let cosT: Float32Array | null = null;
  let sinT: Float32Array | null = null;
  let cxArr: Float32Array | null = null;
  let cyArr: Float32Array | null = null;
  let dxArr: Float32Array | null = null;
  let dyArr: Float32Array | null = null;
  let hasT: Uint8Array | null = null;
  if (haveTransforms && uniqueValues && valueIndices) {
    cosT = new Float32Array(uniqueValues.length);
    sinT = new Float32Array(uniqueValues.length);
    cxArr = new Float32Array(uniqueValues.length);
    cyArr = new Float32Array(uniqueValues.length);
    dxArr = new Float32Array(uniqueValues.length);
    dyArr = new Float32Array(uniqueValues.length);
    hasT = new Uint8Array(uniqueValues.length);
    for (let k = 0; k < uniqueValues.length; k++) {
      const t = transforms.get(uniqueValues[k]);
      if (!t || isIdentity(t)) continue;
      const c = centroids!.get(uniqueValues[k]);
      if (!c) continue;
      cosT[k] = Math.cos(t.theta);
      sinT[k] = Math.sin(t.theta);
      cxArr[k] = c[0];
      cyArr[k] = c[1];
      dxArr[k] = t.dx;
      dyArr[k] = t.dy;
      hasT[k] = 1;
    }
  }

  for (let i = 0; i < n; i++) {
    let x = getX(i);
    let y = getY(i);
    if (haveTransforms) {
      let cx = 0,
        cy = 0,
        dx = 0,
        dy = 0,
        ct = 1,
        st = 0,
        hit = false;
      if (uniqueValues && valueIndices && hasT) {
        const k = valueIndices[i];
        if (hasT[k]) {
          cx = cxArr![k];
          cy = cyArr![k];
          dx = dxArr![k];
          dy = dyArr![k];
          ct = cosT![k];
          st = sinT![k];
          hit = true;
        }
      } else if (cluster && centroids) {
        const s = getSampleId(cluster, i);
        const t = transforms.get(s);
        const c = centroids.get(s);
        if (t && c && !isIdentity(t)) {
          cx = c[0];
          cy = c[1];
          dx = t.dx;
          dy = t.dy;
          ct = Math.cos(t.theta);
          st = Math.sin(t.theta);
          hit = true;
        }
      }
      if (hit) {
        const ox = x - cx;
        const oy = y - cy;
        x = ox * ct - oy * st + cx + dx;
        y = ox * st + oy * ct + cy + dy;
      }
    }
    out[i * dims] = x;
    out[i * dims + 1] = y;
    if (getZ) out[i * dims + 2] = getZ(i);
  }

  // Pack header (8 bytes) + float32 body
  const header = new ArrayBuffer(8);
  const headerView = new DataView(header);
  headerView.setUint32(0, n, true);
  headerView.setUint32(4, dims, true);
  const raw = new Uint8Array(8 + out.byteLength);
  raw.set(new Uint8Array(header), 0);
  raw.set(new Uint8Array(out.buffer), 8);

  // Gzip via CompressionStream (Chromium, Firefox 113+, Safari 16.4+)
  if (typeof CompressionStream === "undefined") {
    throw new Error(
      "CompressionStream API not available in this browser — cannot gzip.",
    );
  }
  const cs = new CompressionStream("gzip");
  const stream = new Response(raw).body!.pipeThrough(cs);
  return new Response(stream).blob();
}
