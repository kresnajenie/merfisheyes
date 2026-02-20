/**
 * Normalize coordinates to the range [-1, 1] centered at 0
 * First centers the data by subtracting the mean, then scales by the maximum absolute value
 * @param coordinates - Array of coordinate arrays
 * @returns Object containing normalized coordinates, the scaling factor, and the center point
 */
export function normalizeCoordinates(coordinates: number[][]): {
  normalized: number[][];
  scalingFactor: number;
  center?: number[];
} | null {
  if (!coordinates || coordinates.length === 0) {
    return null;
  }

  const n = coordinates[0].length;

  // Calculate the mean (center) of each dimension
  const center = new Array(n).fill(0);

  for (const coord of coordinates) {
    for (let i = 0; i < n; i++) {
      center[i] += coord[i];
    }
  }
  for (let i = 0; i < n; i++) {
    center[i] /= coordinates.length;
  }

  // Center the coordinates by subtracting the mean
  const centeredCoords = coordinates.map((coord) =>
    coord.slice(0, n).map((c, i) => c - center[i]),
  );

  // Find the maximum absolute value after centering
  const maxAbs = centeredCoords.reduce((max, coord) => {
    return Math.max(max, ...coord.map((c) => Math.abs(c)));
  }, 0);

  if (maxAbs === 0) {
    return {
      normalized: centeredCoords,
      scalingFactor: 1,
      center,
    };
  }

  // Scale to [-1, 1]
  const normalized = centeredCoords.map((coord) =>
    coord.map((c) => c / maxAbs),
  );

  return {
    normalized,
    scalingFactor: maxAbs,
    center,
  };
}
