/**
 * Normalize coordinates to the range [-1, 1] centered at 0
 * First centers the data by subtracting the mean, then scales by the maximum absolute value
 * @param coordinates - Array of coordinate arrays
 * @param precomputedScalingFactor - Optional pre-computed scaling factor (max abs of centered coords).
 *   If provided, skips the expensive max-abs calculation over all coordinates.
 * @returns Object containing normalized coordinates (rounded to 2 decimal places),
 *   the scaling factor, the center point, and the max raw coordinate
 */
export function normalizeCoordinates(
  coordinates: number[][],
  precomputedScalingFactor?: number,
): {
  normalized: number[][];
  scalingFactor: number;
  center?: number[];
  maxRawCoordinate?: number;
} | null {
  if (!coordinates || coordinates.length === 0) {
    return null;
  }

  const n = coordinates[0].length;

  // Track max raw coordinate before any transforms
  let maxRawCoordinate = 0;

  // Calculate the mean (center) of each dimension
  const center = new Array(n).fill(0);

  for (const coord of coordinates) {
    for (let i = 0; i < n; i++) {
      center[i] += coord[i];
      const absVal = Math.abs(coord[i]);

      if (absVal > maxRawCoordinate) {
        maxRawCoordinate = absVal;
      }
    }
  }
  for (let i = 0; i < n; i++) {
    center[i] /= coordinates.length;
  }

  // Center the coordinates by subtracting the mean
  const centeredCoords = coordinates.map((coord) =>
    coord.slice(0, n).map((c, i) => c - center[i]),
  );

  // Use pre-computed scaling factor or calculate from data
  let maxAbs: number;

  if (precomputedScalingFactor && precomputedScalingFactor > 0) {
    maxAbs = precomputedScalingFactor;
  } else {
    // Find the maximum absolute value after centering
    maxAbs = centeredCoords.reduce((max, coord) => {
      return Math.max(max, ...coord.map((c) => Math.abs(c)));
    }, 0);
  }

  if (maxAbs === 0) {
    return {
      normalized: centeredCoords,
      scalingFactor: 1,
      center,
      maxRawCoordinate,
    };
  }

  // Scale to [-1, 1] and round to 2 decimal places
  const normalized = centeredCoords.map((coord) =>
    coord.map((c) => Math.round((c / maxAbs) * 100) / 100),
  );

  return {
    normalized,
    scalingFactor: maxAbs,
    center,
    maxRawCoordinate,
  };
}
