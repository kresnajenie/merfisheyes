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

/**
 * Normalize a flat Float32Array of coordinates to [-1, 1] range.
 * Much more memory efficient than the number[][] version for large datasets.
 */
export function normalizeCoordinatesFlat(
  coords: Float32Array,
  dimensions: number,
): {
  normalized: Float32Array;
  scalingFactor: number;
} | null {
  const numPoints = coords.length / dimensions;

  if (numPoints === 0) return null;

  // Calculate center (mean of each dimension)
  const center = new Float64Array(dimensions);

  for (let i = 0; i < numPoints; i++) {
    for (let d = 0; d < dimensions; d++) {
      center[d] += coords[i * dimensions + d];
    }
  }
  for (let d = 0; d < dimensions; d++) {
    center[d] /= numPoints;
  }

  // Find max absolute value after centering
  let maxAbs = 0;

  for (let i = 0; i < coords.length; i++) {
    const centered = Math.abs(coords[i] - center[i % dimensions]);

    if (centered > maxAbs) maxAbs = centered;
  }

  if (maxAbs === 0) maxAbs = 1;

  // Normalize in-place to a new Float32Array
  const normalized = new Float32Array(coords.length);

  for (let i = 0; i < coords.length; i++) {
    normalized[i] = (coords[i] - center[i % dimensions]) / maxAbs;
  }

  return { normalized, scalingFactor: maxAbs };
}
