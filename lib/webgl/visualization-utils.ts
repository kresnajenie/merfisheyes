import type { StandardizedDataset } from "../StandardizedDataset";

import { getClusterValue } from "../StandardizedDataset";
import {
  VISUALIZATION_CONFIG,
} from "../config/visualization.config";
import { colormapRgb } from "../utils/colormaps";

export interface AdvancedVizSettings {
  selectedSizeMultiplier: number;
  greyedOutSizeMultiplier: number;
  greyedOutAlpha: number;
  expressionAlphaMin: number;
  expressionAlphaMax: number;
  pointSizeMultiplierMin: number;
  pointSizeMultiplierMax: number;
}

const defaultAdvanced: AdvancedVizSettings = {
  selectedSizeMultiplier: VISUALIZATION_CONFIG.SELECTED_SIZE_MULTIPLIER as number,
  greyedOutSizeMultiplier: VISUALIZATION_CONFIG.GREYED_OUT_SIZE_MULTIPLIER as number,
  greyedOutAlpha: VISUALIZATION_CONFIG.GREYED_OUT_ALPHA as number,
  expressionAlphaMin: VISUALIZATION_CONFIG.EXPRESSION_ALPHA_MIN as number,
  expressionAlphaMax: VISUALIZATION_CONFIG.EXPRESSION_ALPHA_MAX as number,
  pointSizeMultiplierMin: VISUALIZATION_CONFIG.POINT_SIZE_MULTIPLIER_MIN as number,
  pointSizeMultiplierMax: VISUALIZATION_CONFIG.POINT_SIZE_MULTIPLIER_MAX as number,
};

function calcSizeMultiplier(normalizedValue: number, adv: AdvancedVizSettings): number {
  if (isNaN(normalizedValue)) return 1.0;
  return adv.pointSizeMultiplierMin + normalizedValue * (adv.pointSizeMultiplierMax - adv.pointSizeMultiplierMin);
}

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
  scaleMin: number = 0,
  scaleMax: number = 3,
  setScaleMin?: (min: number) => void,
  setScaleMax?: (max: number) => void,
  adv: AdvancedVizSettings = defaultAdvanced,
  colormap: string = "bwr",
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

  // Calculate 95th percentile for auto-scaling
  const percentile95 = calculateGenePercentile(
    expression,
    VISUALIZATION_CONFIG.GENE_EXPRESSION_PERCENTILE,
  );

  console.log("95th percentile:", percentile95);

  // Auto-update scale range: min = 0, max = 95th percentile
  if (setScaleMin && setScaleMax) {
    setScaleMin(0);
    setScaleMax(percentile95);
    console.log("Auto-scaled gene range to:", 0, "-", percentile95);
  }

  // Use the manual scale values from store (scaleMin and scaleMax parameters)
  console.log("Using scale range for visualization:", scaleMin, "-", scaleMax);

  // Normalize expression values to 0-1 range using manual scale values
  const normalizedExpression = expression.map((value) => {
    if (isNaN(value) || value === null) {
      return NaN;
    }
    // Map [scaleMin, scaleMax] to [0, 1]
    const normalized = (value - scaleMin) / (scaleMax - scaleMin);

    // Clamp to [0, 1]
    return Math.max(0, Math.min(1, normalized));
  });

  console.log("Normalized expression:", normalizedExpression.slice(0, 10));

  // Initialize arrays
  const colors = new Float32Array(count * 3);
  const sizes = new Float32Array(count);
  const alphas = new Float32Array(count);

  const baseSize = VISUALIZATION_CONFIG.POINT_BASE_SIZE;
  const baseAlpha = VISUALIZATION_CONFIG.POINT_BASE_ALPHA;

  // Apply colormap with expression-based alpha
  for (let i = 0; i < count; i++) {
    const normalizedValue = normalizedExpression[i];

    const [r, g, b] = colormapRgb(colormap, normalizedValue);
    colors[i * 3] = r;
    colors[i * 3 + 1] = g;
    colors[i * 3 + 2] = b;

    sizes[i] = baseSize * calcSizeMultiplier(normalizedValue, adv);

    const expressionAlpha = isNaN(normalizedValue)
      ? adv.expressionAlphaMin
      : adv.expressionAlphaMin + (adv.expressionAlphaMax - adv.expressionAlphaMin) * normalizedValue;
    alphas[i] = expressionAlpha * alphaScale;
  }

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
  scaleMin: number = 0,
  scaleMax: number = 3,
  setScaleMin?: (min: number) => void,
  setScaleMax?: (max: number) => void,
  adv: AdvancedVizSettings = defaultAdvanced,
  colormap: string = "bwr",
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

  // Get the numerical values using indexed lookup
  const numCount = selectedCluster.valueIndices
    ? selectedCluster.valueIndices.length
    : selectedCluster.values.length;
  const values = new Float32Array(numCount);

  for (let i = 0; i < numCount; i++) {
    values[i] = Number(getClusterValue(selectedCluster, i));
  }

  console.log("Numerical celltype values loaded:", values.slice(0, 10));

  // Calculate 95th percentile for auto-scaling
  const percentile95 = calculateGenePercentile(
    Array.from(values),
    VISUALIZATION_CONFIG.NUMERICAL_CLUSTER_PERCENTILE,
  );

  console.log("95th percentile:", percentile95);

  // Auto-update scale range: min = 0, max = 95th percentile
  if (setScaleMin && setScaleMax) {
    setScaleMin(0);
    setScaleMax(percentile95);
    console.log("Auto-scaled numerical range to:", 0, "-", percentile95);
  }

  // Use the manual scale values from store (scaleMin and scaleMax parameters)
  console.log("Using scale range for visualization:", scaleMin, "-", scaleMax);

  // Normalize values to 0-1 range using manual scale values
  const normalizedValues = new Float32Array(numCount);

  for (let i = 0; i < numCount; i++) {
    const value = values[i];

    if (isNaN(value)) {
      normalizedValues[i] = NaN;
    } else {
      const normalized = (value - scaleMin) / (scaleMax - scaleMin);

      normalizedValues[i] = Math.max(0, Math.min(1, normalized));
    }
  }

  console.log("Normalized values:", normalizedValues.slice(0, 10));

  // Initialize arrays
  const colors = new Float32Array(count * 3);
  const sizes = new Float32Array(count);
  const alphas = new Float32Array(count);

  const baseSize = VISUALIZATION_CONFIG.POINT_BASE_SIZE;
  const baseAlpha = VISUALIZATION_CONFIG.POINT_BASE_ALPHA;

  // Apply colormap with value-based alpha
  for (let i = 0; i < count; i++) {
    const normalizedValue = normalizedValues[i];

    const [r, g, b] = colormapRgb(colormap, normalizedValue);
    colors[i * 3] = r;
    colors[i * 3 + 1] = g;
    colors[i * 3 + 2] = b;

    sizes[i] = baseSize * calcSizeMultiplier(normalizedValue, adv);

    const valueAlpha = isNaN(normalizedValue)
      ? adv.expressionAlphaMin
      : adv.expressionAlphaMin + (adv.expressionAlphaMax - adv.expressionAlphaMin) * normalizedValue;
    alphas[i] = valueAlpha * alphaScale;
  }

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
  adv: AdvancedVizSettings = defaultAdvanced,
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

  const palette = selectedCluster.palette || colorPalette;

  // Initialize arrays
  const colors = new Float32Array(count * 3);
  const sizes = new Float32Array(count);
  const alphas = new Float32Array(count);

  const baseSize = VISUALIZATION_CONFIG.POINT_BASE_SIZE;
  const baseAlpha = VISUALIZATION_CONFIG.POINT_BASE_ALPHA;

  // Apply colors and sizes based on celltype selection
  for (let i = 0; i < count; i++) {
    const category = getClusterValue(selectedCluster, i);
    const isSelected =
      selectedCelltypes.size === 0 || selectedCelltypes.has(category);

    // Colors
    const hex = isSelected ? palette[category] || grey : grey;
    const [r, g, b] = hexToRgb(hex);

    colors[i * 3] = r;
    colors[i * 3 + 1] = g;
    colors[i * 3 + 2] = b;

    // Sizes - selected items are bigger
    const sizeMultiplier = isSelected
      ? adv.selectedSizeMultiplier
      : adv.greyedOutSizeMultiplier;

    sizes[i] = baseSize * sizeMultiplier;

    // Alphas
    alphas[i] = isSelected
      ? baseAlpha * alphaScale
      : adv.greyedOutAlpha * alphaScale;
  }

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
  scaleMin: number = 0,
  scaleMax: number = 3,
  setScaleMin?: (min: number) => void,
  setScaleMax?: (max: number) => void,
  adv: AdvancedVizSettings = defaultAdvanced,
  colormap: string = "bwr",
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
    console.log(
      "No cluster data or column selected for combined visualization",
    );

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

  // Calculate 95th percentile for auto-scaling
  const percentile95 = calculateGenePercentile(
    expression,
    VISUALIZATION_CONFIG.GENE_EXPRESSION_PERCENTILE,
  );

  console.log("95th percentile:", percentile95);

  // Auto-update scale range: min = 0, max = 95th percentile
  if (setScaleMin && setScaleMax) {
    setScaleMin(0);
    setScaleMax(percentile95);
    console.log("Auto-scaled gene range to:", 0, "-", percentile95);
  }

  // Use the manual scale values from store (scaleMin and scaleMax parameters)
  console.log("Using scale range for visualization:", scaleMin, "-", scaleMax);

  // Normalize expression values to 0-1 range using manual scale values
  const normalizedExpression = expression.map((value) => {
    if (isNaN(value) || value === null) {
      return NaN;
    }
    // Map [scaleMin, scaleMax] to [0, 1]
    const normalized = (value - scaleMin) / (scaleMax - scaleMin);

    // Clamp to [0, 1]
    return Math.max(0, Math.min(1, normalized));
  });

  // Initialize arrays
  const colors = new Float32Array(count * 3);
  const sizes = new Float32Array(count);
  const alphas = new Float32Array(count);

  const baseSize = VISUALIZATION_CONFIG.POINT_BASE_SIZE;
  const baseAlpha = VISUALIZATION_CONFIG.POINT_BASE_ALPHA;

  // Apply combined visualization
  for (let i = 0; i < count; i++) {
    const category = getClusterValue(selectedCluster, i);
    const isSelected =
      selectedCelltypes.size === 0 || selectedCelltypes.has(category);

    if (isSelected) {
      // Selected celltypes: show gene expression gradient
      const normalizedValue = normalizedExpression[i];

      // Colors from active colormap
      const [r, g, b] = colormapRgb(colormap, normalizedValue);

      colors[i * 3] = r;
      colors[i * 3 + 1] = g;
      colors[i * 3 + 2] = b;

      // Sizes based on expression level (higher expression = bigger)
      sizes[i] = baseSize * calcSizeMultiplier(normalizedValue, adv);

      // Alpha based on expression level
      const expressionAlpha = isNaN(normalizedValue)
        ? adv.expressionAlphaMin
        : adv.expressionAlphaMin + (adv.expressionAlphaMax - adv.expressionAlphaMin) * normalizedValue;
      alphas[i] = expressionAlpha * alphaScale;
    } else {
      // Non-selected celltypes: show grey with reduced alpha
      const [r, g, b] = hexToRgb(grey);

      colors[i * 3] = r;
      colors[i * 3 + 1] = g;
      colors[i * 3 + 2] = b;

      sizes[i] = baseSize * adv.greyedOutSizeMultiplier;

      alphas[i] = adv.greyedOutAlpha * alphaScale;
    }
  }

  return { colors, sizes, alphas };
}
