/**
 * Fixed palette of 10 distinct colors for single molecule gene visualization.
 * Inspired by fluorescence microscopy imaging conventions.
 *
 * The first 5 colors (cyan, magenta, yellow, green, orange) are maximally
 * distinct from each other. The remaining 5 fill in hue gaps.
 *
 * Colors are assigned by selection order using lowest-available-slot filling.
 * After 10 genes, the palette cycles with progressively lighter variants.
 */

const GENE_COLOR_PALETTE = [
  "#00FFFF", // 0  Cyan        (DAPI / nuclei)
  "#FF00FF", // 1  Magenta     (Cy5 / far-red)
  "#FFFF00", // 2  Yellow      (YFP)
  "#00FF00", // 3  Green       (GFP)
  "#FF8800", // 4  Orange      (mOrange)
  "#FF4444", // 5  Red         (RFP / Texas Red)
  "#4488FF", // 6  Blue        (CFP)
  "#00FF88", // 7  Spring Green
  "#FF0088", // 8  Hot Pink
  "#AA44FF", // 9  Violet
];

/**
 * Returns the color for a given slot index.
 * Slots 0-9 return the base palette colors.
 * Slots 10+ cycle through the palette with progressively lighter tints
 * (blended toward white by 25% per cycle, capped at 60%).
 */
export function getColorForSlot(slot: number): string {
  const baseIndex = slot % GENE_COLOR_PALETTE.length;
  const cycle = Math.floor(slot / GENE_COLOR_PALETTE.length);
  const baseColor = GENE_COLOR_PALETTE[baseIndex];

  if (cycle === 0) return baseColor;

  // Blend toward white for subsequent cycles
  const r = parseInt(baseColor.slice(1, 3), 16);
  const g = parseInt(baseColor.slice(3, 5), 16);
  const b = parseInt(baseColor.slice(5, 7), 16);

  const blend = Math.min(cycle * 0.25, 0.6);
  const nr = Math.round(r + (255 - r) * blend);
  const ng = Math.round(g + (255 - g) * blend);
  const nb = Math.round(b + (255 - b) * blend);

  return `#${nr.toString(16).padStart(2, "0")}${ng.toString(16).padStart(2, "0")}${nb.toString(16).padStart(2, "0")}`;
}

/**
 * Finds the lowest slot index not currently in use.
 */
export function findLowestAvailableSlot(usedSlots: Set<number>): number {
  let slot = 0;

  while (usedSlots.has(slot)) slot++;

  return slot;
}
