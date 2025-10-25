import type { StandardizedDataset } from "../StandardizedDataset";

/**
 * Calculates the value at the specified percentile of the given array, ignoring NaN values.
 * @param arr - The array of numerical values.
 * @param percentile - The percentile to calculate (between 0 and 1).
 * @returns The value at the specified percentile.
 */
export function calculateGenePercentile(
  arr: number[],
  percentile: number,
): number {
  // Filter out NaN values and create a sorted copy
  const sortedArr = arr
    .filter((value) => !isNaN(value) && value !== null)
    .sort((a, b) => a - b);

  if (sortedArr.length === 0) {
    return NaN; // Return NaN if all values are NaN
  }

  // Calculate the index for the xth percentile
  const index = Math.floor(sortedArr.length * percentile);

  // Ensure index is within bounds
  const boundedIndex = Math.min(Math.max(0, index), sortedArr.length - 1);

  return sortedArr[boundedIndex];
}

/**
 * Normalizes the values in the array to a range between 0 and 1.
 * Preserves NaN values and handles division by zero.
 * @param arr - The array of numerical values.
 * @param nmax - The maximum value in the array.
 * @returns The array with normalized values.
 */
export function normalizeArray(arr: number[], nmax: number): number[] {
  // Handle edge case where nmax is 0, NaN, or null
  if (!nmax || isNaN(nmax)) {
    return arr.map(() => 0);
  }

  return arr.map((value) => {
    // Preserve NaN and null values
    if (isNaN(value) || value === null) {
      return NaN;
    }

    // Clip and normalize
    return Math.min(value / nmax, 1);
  });
}

/**
 * Generates a color value in the coolwarm colormap based on the input value.
 * Returns RGB values normalized to 0-1 range.
 * @param value - The value for which to generate the color (between 0 and 1).
 * @returns The color as [r, g, b] tuple (0-1 range).
 */
export function coolwarm(value: number): [number, number, number] {
  // Check for NaN values - return white
  if (isNaN(value)) {
    return [1, 1, 1]; // White for NaN
  }

  // Define start and end colors (cool: blue, warm: red)
  // Values in 0-1 range for WebGL
  const blue = { r: 0, g: 0, b: 1 }; // Blue
  const white = { r: 1, g: 1, b: 1 }; // White
  const red = { r: 1, g: 0, b: 0 }; // Red

  if (value < 0.5) {
    // Blue to white
    const t = value * 2; // 0 to 1

    return [white.r * t, white.g * t, blue.b];
  } else if (value === 0.5) {
    // White
    return [white.r, white.g, white.b];
  } else {
    // White to red
    const t = (value - 0.5) * 2; // 0 to 1

    return [red.r, white.g - white.g * t, white.b - white.b * t];
  }
}

/**
 * Converts hex color to RGB array (0-1 range)
 */

const grey = "#808080";

export function hexToRgb(hex: string): [number, number, number] {
  const result = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex);

  if (!result) return [0, 0, 0];

  return [
    parseInt(result[1], 16) / 255,
    parseInt(result[2], 16) / 255,
    parseInt(result[3], 16) / 255,
  ];
}

/**
 * Generates a color gradient from min (blue) to max (red)
 */
export function valueToColor(
  value: number,
  min: number,
  max: number,
): [number, number, number] {
  const normalized = max > min ? (value - min) / (max - min) : 0;
  const clamped = Math.max(0, Math.min(1, normalized));

  // Blue (low) -> Cyan -> Green -> Yellow -> Red (high)
  if (clamped < 0.25) {
    const t = clamped / 0.25;

    return [0, t, 1]; // Blue to Cyan
  } else if (clamped < 0.5) {
    const t = (clamped - 0.25) / 0.25;

    return [0, 1, 1 - t]; // Cyan to Green
  } else if (clamped < 0.75) {
    const t = (clamped - 0.5) / 0.25;

    return [t, 1, 0]; // Green to Yellow
  } else {
    const t = (clamped - 0.75) / 0.25;

    return [1, 1 - t, 0]; // Yellow to Red
  }
}

/**
 * Update visualization for gene expression mode
 */
export async function updateGeneVisualization(
  dataset: StandardizedDataset,
  selectedGene: string | null,
  alphaScale: number,
  sizeScale: number,
): Promise<{
  colors: Float32Array;
  sizes: Float32Array;
  alphas: Float32Array;
} | null> {
  const count = dataset.getPointCount();

  if (!selectedGene) {
    console.log("No gene selected");

    return null;
  }

  // Fetch gene expression data
  const expression = await dataset.getGeneExpression(selectedGene);

  if (!expression) {
    console.warn(`Gene expression data not found for: ${selectedGene}`);

    return null;
  }

  console.log("Expression data loaded for gene:", selectedGene);

  // Calculate 95th percentile for normalization
  const percentile95 = calculateGenePercentile(expression, 0.95);

  console.log("95th percentile:", percentile95);

  // Normalize expression values to 0-1 range
  const normalizedExpression = normalizeArray(expression, percentile95);

  console.log("Normalized expression:", normalizedExpression.slice(0, 10));

  // Initialize arrays
  const colors = new Float32Array(count * 3);
  const sizes = new Float32Array(count);
  const alphas = new Float32Array(count);

  const baseSize = 2.0;
  const baseAlpha = 1.0;

  // Apply coolwarm colormap and use normalized values for size
  for (let i = 0; i < count; i++) {
    const normalizedValue = normalizedExpression[i];

    // Colors from coolwarm
    const [r, g, b] = coolwarm(normalizedValue);

    colors[i * 3] = r;
    colors[i * 3 + 1] = g;
    colors[i * 3 + 2] = b;

    // Sizes based on expression level (higher expression = bigger)
    // Use normalized value to scale size (0.5 to 2x base size)
    const sizeMultiplier = isNaN(normalizedValue)
      ? 1.0
      : 0.5 + normalizedValue * 1.5;

    sizes[i] = baseSize * sizeMultiplier * sizeScale;

    // Alphas
    alphas[i] = baseAlpha * alphaScale;
  }

  console.log("Gene visualization updated successfully");

  return { colors, sizes, alphas };
}

/**
 * Update visualization for numerical celltype mode (uses same logic as gene expression)
 */
export function updateNumericalCelltypeVisualization(
  dataset: StandardizedDataset,
  selectedColumn: string | null,
  alphaScale: number,
  sizeScale: number,
): {
  colors: Float32Array;
  sizes: Float32Array;
  alphas: Float32Array;
} | null {
  const count = dataset.getPointCount();

  if (!dataset.clusters || !selectedColumn) {
    console.log("No cluster data or column selected");

    return null;
  }

  // Find the selected cluster column data
  const selectedCluster = dataset.clusters.find(
    (c) => c.column === selectedColumn,
  );

  if (!selectedCluster) {
    console.warn(`Cluster column not found: ${selectedColumn}`);

    return null;
  }

  // Get the numerical values
  const values = selectedCluster.values.map((v) => Number(v));

  console.log("Numerical celltype values loaded:", values.slice(0, 10));

  // Calculate 95th percentile for normalization (same as gene expression)
  const percentile95 = calculateGenePercentile(values, 0.95);

  console.log("95th percentile:", percentile95);

  // Normalize values to 0-1 range (same as gene expression)
  const normalizedValues = normalizeArray(values, percentile95);

  console.log("Normalized values:", normalizedValues.slice(0, 10));

  // Initialize arrays
  const colors = new Float32Array(count * 3);
  const sizes = new Float32Array(count);
  const alphas = new Float32Array(count);

  const baseSize = 2.0;
  const baseAlpha = 1.0;

  // Apply coolwarm colormap and use normalized values for size (same as gene expression)
  for (let i = 0; i < count; i++) {
    const normalizedValue = normalizedValues[i];

    // Colors from coolwarm
    const [r, g, b] = coolwarm(normalizedValue);

    colors[i * 3] = r;
    colors[i * 3 + 1] = g;
    colors[i * 3 + 2] = b;

    // Sizes based on value level (higher value = bigger)
    // Use normalized value to scale size (0.5 to 2x base size)
    const sizeMultiplier = isNaN(normalizedValue)
      ? 1.0
      : 0.5 + normalizedValue * 1.5;

    sizes[i] = baseSize * sizeMultiplier * sizeScale;

    // Alphas
    alphas[i] = baseAlpha * alphaScale;
  }

  console.log("Numerical celltype visualization updated successfully");

  return { colors, sizes, alphas };
}

/**
 * Update visualization for celltype mode
 */
export function updateCelltypeVisualization(
  dataset: StandardizedDataset,
  selectedColumn: string | null,
  selectedCelltypes: Set<string>,
  colorPalette: Record<string, string>,
  alphaScale: number,
  sizeScale: number,
): {
  colors: Float32Array;
  sizes: Float32Array;
  alphas: Float32Array;
} | null {
  const count = dataset.getPointCount();

  if (!dataset.clusters || !selectedColumn) {
    console.log("No cluster data or column selected");

    return null;
  }

  // Find the selected cluster column data
  const selectedCluster = dataset.clusters.find(
    (c) => c.column === selectedColumn,
  );

  if (!selectedCluster) {
    console.warn(`Cluster column not found: ${selectedColumn}`);

    return null;
  }

  const values = selectedCluster.values;
  const palette = selectedCluster.palette || colorPalette;

  // Initialize arrays
  const colors = new Float32Array(count * 3);
  const sizes = new Float32Array(count);
  const alphas = new Float32Array(count);

  const baseSize = 2.0;
  const baseAlpha = 1.0;

  // Apply colors and sizes based on celltype selection
  for (let i = 0; i < count; i++) {
    const category = String(values[i]);
    const isSelected =
      selectedCelltypes.size === 0 || selectedCelltypes.has(category);

    // Colors
    const hex = isSelected ? palette[category] || grey : grey;
    const [r, g, b] = hexToRgb(hex);

    colors[i * 3] = r;
    colors[i * 3 + 1] = g;
    colors[i * 3 + 2] = b;

    // Sizes - selected items are bigger
    const sizeMultiplier =
      selectedCelltypes.size === 0 || selectedCelltypes.has(category)
        ? 2.0
        : 1.0;

    sizes[i] = baseSize * sizeMultiplier * sizeScale;

    // Alphas
    alphas[i] = baseAlpha * alphaScale;
  }

  console.log("Celltype visualization updated successfully");

  return { colors, sizes, alphas };
}

/**
 * Updates visualization with combined gene expression + celltype filtering
 * Shows gene expression gradient on selected celltypes, greys out others
 */
export async function updateCombinedVisualization(
  dataset: StandardizedDataset,
  selectedGene: string | null,
  selectedColumn: string | null,
  selectedCelltypes: Set<string>,
  alphaScale: number,
  sizeScale: number,
): Promise<{
  colors: Float32Array;
  sizes: Float32Array;
  alphas: Float32Array;
} | null> {
  const count = dataset.getPointCount();

  if (!selectedGene) {
    console.log("No gene selected for combined visualization");
    return null;
  }

  if (!dataset.clusters || !selectedColumn) {
    console.log("No cluster data or column selected for combined visualization");
    return null;
  }

  // Fetch gene expression data
  const expression = await dataset.getGeneExpression(selectedGene);

  if (!expression) {
    console.warn(`Gene expression data not found for: ${selectedGene}`);
    return null;
  }

  // Find the selected cluster column data
  const selectedCluster = dataset.clusters.find(
    (c) => c.column === selectedColumn,
  );

  if (!selectedCluster) {
    console.warn(`Cluster column not found: ${selectedColumn}`);
    return null;
  }

  const clusterValues = selectedCluster.values;

  // Calculate 95th percentile for normalization
  const percentile95 = calculateGenePercentile(expression, 0.95);

  // Normalize expression values to 0-1 range
  const normalizedExpression = normalizeArray(expression, percentile95);

  // Initialize arrays
  const colors = new Float32Array(count * 3);
  const sizes = new Float32Array(count);
  const alphas = new Float32Array(count);

  const baseSize = 2.0;
  const baseAlpha = 1.0;

  // Apply combined visualization
  for (let i = 0; i < count; i++) {
    const category = String(clusterValues[i]);
    const isSelected =
      selectedCelltypes.size === 0 || selectedCelltypes.has(category);

    if (isSelected) {
      // Selected celltypes: show gene expression gradient
      const normalizedValue = normalizedExpression[i];

      // Colors from coolwarm gradient
      const [r, g, b] = coolwarm(normalizedValue);

      colors[i * 3] = r;
      colors[i * 3 + 1] = g;
      colors[i * 3 + 2] = b;

      // Sizes based on expression level (higher expression = bigger)
      const sizeMultiplier = isNaN(normalizedValue)
        ? 1.0
        : 0.5 + normalizedValue * 1.5;

      sizes[i] = baseSize * sizeMultiplier * sizeScale;

      // Alphas based on gene expression
      alphas[i] = baseAlpha * alphaScale;
    } else {
      // Non-selected celltypes: show grey with reduced size/alpha
      const [r, g, b] = hexToRgb(grey);

      colors[i * 3] = r;
      colors[i * 3 + 1] = g;
      colors[i * 3 + 2] = b;

      // Smaller size for non-selected
      sizes[i] = baseSize * 1.0 * sizeScale;

      // Same alpha as celltype mode
      alphas[i] = baseAlpha * alphaScale;
    }
  }

  console.log("Combined gene + celltype visualization updated successfully");

  return { colors, sizes, alphas };
}
