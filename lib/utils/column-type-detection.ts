/**
 * Utility functions for detecting column types (categorical vs numerical)
 * Used across all data adapters and preprocessing pipelines
 */

/**
 * Check if a value looks like a number (can be parsed as float)
 */
function looksLikeNumber(value: any): boolean {
  if (typeof value === "number") return true;
  if (typeof value === "string") {
    const trimmed = value.trim();
    if (trimmed === "" || trimmed === "nan" || trimmed === "NaN") return false;
    const num = parseFloat(trimmed);
    return !isNaN(num) && isFinite(num);
  }
  return false;
}

/**
 * Check if a value looks like a float (has decimal point or scientific notation)
 */
function looksLikeFloat(value: any): boolean {
  if (typeof value === "number") {
    // Check if it has a fractional part
    return !Number.isInteger(value);
  }
  if (typeof value === "string") {
    const trimmed = value.trim();
    // Check for decimal point or scientific notation (e.g., "1.5", "1e-10")
    return /\.|e|E/.test(trimmed);
  }
  return false;
}

/**
 * Determine if data is categorical or numerical
 *
 * Algorithm:
 * 1. Always treat columns named "leiden" or "louvain" as categorical
 * 2. Check if values look like floats → numerical
 * 3. For integer-like values: if ≥80% are unique → numerical, otherwise categorical
 * 4. For non-numeric strings → categorical
 *
 * @param values - Array of values from the column
 * @param columnName - Optional column name (for special cases like "leiden")
 * @returns true if categorical, false if numerical
 */
export function isCategorical(
  values: any[],
  columnName?: string,
): boolean {
  if (!values || values.length === 0) return false;

  // Special case: leiden/louvain columns are always categorical
  if (columnName) {
    const lowerName = columnName.toLowerCase();
    if (lowerName.includes("leiden") || lowerName.includes("louvain")) {
      return true;
    }
  }

  // Filter out null/undefined values for analysis
  const validValues = values.filter((v) => v != null);
  if (validValues.length === 0) return false;

  // Sample up to 1000 values for performance on large datasets
  const sampleSize = Math.min(1000, validValues.length);
  const step = Math.max(1, Math.floor(validValues.length / sampleSize));
  const sample = validValues.filter((_, i) => i % step === 0).slice(0, sampleSize);

  // Check if any values look like floats
  let floatCount = 0;
  let numberCount = 0;

  for (const value of sample) {
    if (looksLikeFloat(value)) {
      floatCount++;
    }
    if (looksLikeNumber(value)) {
      numberCount++;
    }
  }

  const floatRatio = floatCount / sample.length;
  const numberRatio = numberCount / sample.length;

  // If majority (>50%) look like floats → numerical
  if (floatRatio > 0.5) {
    return false;
  }

  // If not numbers at all → categorical (string labels)
  if (numberRatio < 0.5) {
    return true;
  }

  // At this point, values are integer-like numbers or integer strings
  // Check uniqueness ratio: if ≥80% unique → numerical, otherwise categorical
  const uniqueValues = new Set(validValues.map((v) => String(v)));
  const uniqueRatio = uniqueValues.size / validValues.length;

  return uniqueRatio < 0.8;
}

/**
 * Legacy function name for backward compatibility
 */
export function isCategoricalData(
  values: any[],
  columnName?: string,
): boolean {
  return isCategorical(values, columnName);
}
