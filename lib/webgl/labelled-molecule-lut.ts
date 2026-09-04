/**
 * Lookup-table construction for the labelled single-molecule viewer.
 *
 * These are the only things that change when a checkbox is toggled: three
 * selection LUTs (one texel per category) and one palette LUT. Everything
 * per-point stays on the GPU untouched.
 */

import { SIZE_ENCODE_RANGE } from "./labelled-molecule-shaders";

/**
 * Selection mask for one column, one byte per category.
 *
 * An EMPTY selection means "no constraint from this menu", so the LUT is all
 * 255 — the shader multiplies the three masks together and an unconstrained
 * menu must contribute 1.0, not 0.0.
 */
export function buildSelectionLut(
  uniqueValues: string[],
  selection: Set<string>,
): Uint8Array {
  const lut = new Uint8Array(uniqueValues.length);

  if (selection.size === 0) {
    lut.fill(255);

    return lut;
  }

  for (let i = 0; i < uniqueValues.length; i++) {
    if (selection.has(uniqueValues[i])) lut[i] = 255;
  }

  return lut;
}

/**
 * Colour + size for every category of the colouring column, as RGBA8:
 * RGB is the colour, A is the size multiplier encoded over
 * [0, SIZE_ENCODE_RANGE].
 */
export function buildPaletteLut(
  uniqueValues: string[],
  colorOf: (value: string) => string,
  sizeOf: (value: string) => number,
): Uint8Array {
  const lut = new Uint8Array(uniqueValues.length * 4);

  for (let i = 0; i < uniqueValues.length; i++) {
    const [r, g, b] = hexToRgb255(colorOf(uniqueValues[i]));
    const size = Math.min(
      SIZE_ENCODE_RANGE,
      Math.max(0, sizeOf(uniqueValues[i])),
    );

    lut[i * 4] = r;
    lut[i * 4 + 1] = g;
    lut[i * 4 + 2] = b;
    lut[i * 4 + 3] = Math.round((size / SIZE_ENCODE_RANGE) * 255);
  }

  return lut;
}

export function hexToRgb255(hex: string): [number, number, number] {
  const h = hex.replace("#", "");
  const full =
    h.length === 3
      ? h[0] + h[0] + h[1] + h[1] + h[2] + h[2]
      : h.padEnd(6, "0").slice(0, 6);
  const n = parseInt(full, 16);

  return [(n >> 16) & 255, (n >> 8) & 255, n & 255];
}

/**
 * How many molecules pass every menu — the same product the shader computes,
 * evaluated on the CPU so the UI can show a count.
 *
 * One pass over the point count with three typed-array lookups. This is far
 * cheaper than the colour/size/alpha rebuild it replaces, but it is still O(N),
 * so call it for display only, not per frame.
 */
export function countVisible(
  columns: { indices: ArrayLike<number>; lut: Uint8Array }[],
  pointCount: number,
): number {
  const active = columns.filter((c) => !c.lut.every((v) => v === 255));

  if (active.length === 0) return pointCount;

  let n = 0;

  outer: for (let i = 0; i < pointCount; i++) {
    for (let c = 0; c < active.length; c++) {
      if (active[c].lut[active[c].indices[i]] === 0) continue outer;
    }
    n++;
  }

  return n;
}
