/**
 * Visualization Configuration
 *
 * This file contains all configurable constants for the visualization system.
 * Modify these values to adjust the appearance and behavior of the 3D visualization.
 */

export const VISUALIZATION_CONFIG = {
  /**
   * Gene Expression Normalization
   * Percentile used to normalize gene expression values (0.0 - 1.0)
   * Higher percentile = more lenient scaling, lower percentile = more aggressive scaling
   * Default: 0.95 (95th percentile)
   */
  GENE_EXPRESSION_PERCENTILE: 0.95,

  /**
   * Numerical Cluster Normalization
   * Percentile used to normalize numerical cluster values (0.0 - 1.0)
   * Default: 0.95 (95th percentile)
   */
  NUMERICAL_CLUSTER_PERCENTILE: 0.95,

  /**
   * Point Sizes
   * Base size for rendered points before any scaling is applied (single cell visualization)
   */
  POINT_BASE_SIZE: 0.5,

  /**
   * Single Molecule Point Base Size
   * Base size for single molecule points (separate from cell visualization)
   * Default: 2.0
   */
  SINGLE_MOLECULE_POINT_BASE_SIZE: 5.0,

  /**
   * Point Size Multiplier Range
   * Controls how much the point size varies based on expression level
   * Formula: baseSize * (sizeMin + normalizedValue * (sizeMax - sizeMin)) * sizeScale
   *
   * sizeMin: Minimum multiplier for low expression (e.g., 0.5 = half size)
   * sizeMax: Maximum multiplier for high expression (e.g., 2.0 = double size)
   *
   * Default: min=0.5, max=2.0 (points range from 50% to 200% of base size)
   */
  POINT_SIZE_MULTIPLIER_MIN: 0.5,
  POINT_SIZE_MULTIPLIER_MAX: 3.0,

  /**
   * Point Alpha (Opacity)
   * Base alpha/opacity for rendered points (0.0 - 1.0)
   * 1.0 = fully opaque, 0.0 = fully transparent
   */
  POINT_BASE_ALPHA: 1.0,

  /**
   * Scale Bar Settings
   * Default min/max values for gene and numerical cluster scales
   */
  SCALE_BAR_DEFAULT_MIN: 0,
  SCALE_BAR_DEFAULT_MAX: 3,

  /**
   * Scale Bar Step Size
   * Percentage of max value used for step size when scrubbing (0.0 - 1.0)
   * Default: 0.001 (0.1% of max value)
   */
  SCALE_BAR_STEP_PERCENTAGE: 0.001,

  /**
   * Scale Bar Minimum Step
   * Minimum step size to prevent too-fine control on very small values
   * Default: 0.001
   */
  SCALE_BAR_MIN_STEP: 0.001,

  /**
   * Scale Bar Display Decimals
   * Number of decimal places to show on scale bar numbers
   * Default: 3
   */
  SCALE_BAR_DECIMALS: 3,
} as const;

/**
 * Helper function to calculate size multiplier based on normalized value
 * @param normalizedValue - Value between 0 and 1
 * @returns Size multiplier
 */
export function calculateSizeMultiplier(normalizedValue: number): number {
  if (isNaN(normalizedValue)) {
    return 1.0;
  }
  return (
    VISUALIZATION_CONFIG.POINT_SIZE_MULTIPLIER_MIN +
    normalizedValue *
      (VISUALIZATION_CONFIG.POINT_SIZE_MULTIPLIER_MAX -
        VISUALIZATION_CONFIG.POINT_SIZE_MULTIPLIER_MIN)
  );
}
