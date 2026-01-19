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
