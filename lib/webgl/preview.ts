/**
 * Tiny renderable previews, written by scripts/spiralia/export_preview.py.
 *
 *   preview/points.bin.gz  u32 version | u32 nPoints | u32 dims
 *                          f32 min[dims] | f32 scale[dims]
 *                          u16 coords[nPoints * dims]
 *                          u8  colorIndex[nPoints]
 *   preview/index.json     palette in colour-index order
 *
 * ~5,000 points and 32 KB, so a rail of them costs less than a tenth of one
 * real dataset. Colour is a palette index rather than RGB: a third of the
 * bytes, and it keeps a preview consistent with the viewer's own colours.
 */

export interface DatasetPreview {
  /** Centred on the origin and scaled to fit a unit-ish box, ready to draw. */
  positions: Float32Array;
  /** Flat RGB triples, one per point, resolved through the palette. */
  colors: Float32Array;
  count: number;
}

function hexToRgb(hex: string): [number, number, number] {
  const h = hex.replace("#", "");
  const full =
    h.length === 3 ? h[0] + h[0] + h[1] + h[1] + h[2] + h[2] : h.padEnd(6, "0");
  const n = parseInt(full.slice(0, 6), 16);

  return [((n >> 16) & 255) / 255, ((n >> 8) & 255) / 255, (n & 255) / 255];
}

/**
 * Fetch one dataset's preview, or null when it has none.
 *
 * Returns null rather than throwing: a dataset without a preview is a normal
 * state, and a rail should simply skip it.
 */
export async function loadPreview(
  baseUrl: string,
  signal?: AbortSignal,
): Promise<DatasetPreview | null> {
  const base = baseUrl.replace(/\/+$/, "");

  try {
    const [idxRes, binRes] = await Promise.all([
      fetch(`${base}/preview/index.json`, { signal }),
      fetch(`${base}/preview/points.bin.gz`, { signal }),
    ]);

    if (!idxRes.ok || !binRes.ok) return null;

    const index = (await idxRes.json()) as { palette: string[] };
    const stream = binRes.body!.pipeThrough(new DecompressionStream("gzip"));
    const buf = await new Response(stream).arrayBuffer();

    const head = new Uint32Array(buf, 0, 3);
    const count = head[1];
    const dims = head[2];

    let off = 12;
    const min = new Float32Array(buf, off, dims);

    off += dims * 4;
    const scale = new Float32Array(buf, off, dims);

    off += dims * 4;
    const codes = new Uint16Array(buf, off, count * dims);

    off += count * dims * 2;
    const colorIndex = new Uint8Array(buf, off, count);

    // Dequantise, then centre and normalise so every preview draws at the same
    // size in a thumbnail regardless of the embryo's real extent.
    const positions = new Float32Array(count * 3);
    const lo = [Infinity, Infinity, Infinity];
    const hi = [-Infinity, -Infinity, -Infinity];

    for (let i = 0; i < count; i++) {
      for (let d = 0; d < 3; d++) {
        const v = d < dims ? min[d] + codes[i * dims + d] * scale[d] : 0;

        positions[i * 3 + d] = v;
        if (v < lo[d]) lo[d] = v;
        if (v > hi[d]) hi[d] = v;
      }
    }

    const centre = lo.map((l, d) => (l + hi[d]) / 2);
    const extent = Math.max(...hi.map((h, d) => h - lo[d]), 1e-6);

    for (let i = 0; i < count; i++) {
      for (let d = 0; d < 3; d++) {
        positions[i * 3 + d] = (positions[i * 3 + d] - centre[d]) / extent;
      }
    }

    const palette = index.palette.map(hexToRgb);
    const colors = new Float32Array(count * 3);

    for (let i = 0; i < count; i++) {
      const [r, g, b] = palette[colorIndex[i]] ?? [0.5, 0.5, 0.5];

      colors[i * 3] = r;
      colors[i * 3 + 1] = g;
      colors[i * 3 + 2] = b;
    }

    return { positions, colors, count };
  } catch {
    return null;
  }
}
