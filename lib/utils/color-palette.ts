/**
 * Shared color palette for cluster visualization
 * Bright, distinct colors optimized for black backgrounds
 */

export const DEFAULT_COLOR_PALETTE = [
  "#1f77b4", // Blue
  "#ff7f0e", // Orange
  "#2ca02c", // Green
  "#d62728", // Red
  "#9467bd", // Purple
  "#00d9ff", // Bright Cyan
  "#e377c2", // Pink
  "#ffeb3b", // Bright Yellow
  "#bcbd22", // Yellow-green
  "#17becf", // Cyan
  "#ff006e", // Bright Magenta
  "#00ff00", // Lime Green
  "#ff1744", // Bright Red
  "#00e5ff", // Bright Aqua
  "#ffc107", // Amber
  "#e91e63", // Deep Pink
  "#4caf50", // Medium Green
  "#ff9800", // Deep Orange
  "#9c27b0", // Deep Purple
  "#00bcd4", // Light Cyan
  "#ffeb3b", // Yellow
  "#f44336", // Red
  "#3f51b5", // Indigo
  "#8bc34a", // Light Green
  "#ff5722", // Orange Red
  "#673ab7", // Deep Purple
  "#03a9f4", // Light Blue
  "#cddc39", // Lime
  "#ff9100", // Orange
  "#7c4dff", // Deep Purple
  "#00e676", // Green
  "#ff3d00", // Deep Orange
  "#651fff", // Deep Purple
  "#1de9b6", // Teal
  "#ff6e40", // Deep Orange
  "#d500f9", // Purple
  "#00b0ff", // Light Blue
  "#76ff03", // Light Green
  "#ff1744", // Red
  "#00e5ff", // Cyan
] as const;

/**
 * Get a color from the palette by index
 * Cycles through the palette if index exceeds palette length
 */
export function getColorFromPalette(index: number): string {
  return DEFAULT_COLOR_PALETTE[index % DEFAULT_COLOR_PALETTE.length];
}

/**
 * Generate a color palette for a list of unique values
 * Returns a map of value -> color
 */
export function generateColorPalette(values: string[]): Record<string, string> {
  const palette: Record<string, string> = {};

  values.forEach((value, index) => {
    palette[value] = getColorFromPalette(index);
  });

  return palette;
}

/**
 * The order a column's labels are listed and dealt colours in, everywhere:
 * as text, with numbers by value (cluster_2 before cluster_10). The browser
 * preview and the SpatialData viewer both sort with it, so a dataset keeps
 * its colours when it is saved (docs/plans/PREVIEW-SAVE.md, D5).
 */
export function compareLabels(a: string, b: string): number {
  return a.localeCompare(b, undefined, { numeric: true });
}

/** A stored colour as a 7-character hex: scanpy may store "#rrggbbaa". */
export function normalizeHexColor(hex: unknown): string {
  if (typeof hex !== "string") return "#808080";

  return hex.length > 7 ? hex.substring(0, 7) : hex;
}

/**
 * A categorical column's colours: the default palette dealt over `labels`
 * (already in compareLabels order), replaced by the colours the file stored
 * for the column — scanpy's `uns[<column>_colors]`, one per category in the
 * column's stored category order — when there is one per category and they
 * are not all the same. `categories` defaults to `labels`: scanpy builds a
 * plain column's categories sorted.
 */
export function columnPalette(
  labels: string[],
  stored?: { colors: unknown[] | null; categories?: string[] | null } | null,
): Record<string, string> {
  const palette = generateColorPalette(labels);
  const colors = stored?.colors?.map(normalizeHexColor);
  const categories = stored?.categories ?? labels;

  if (
    colors &&
    colors.length === categories.length &&
    (new Set(colors).size > 1 || categories.length === 1)
  ) {
    categories.forEach((category, i) => {
      palette[String(category)] = colors[i];
    });
  }

  return palette;
}
