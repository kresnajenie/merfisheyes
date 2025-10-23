/**
 * Normalize coordinates to the range [-1, 1] based on the maximum absolute value
 * @param coordinates - Array of coordinate arrays
 * @returns Object containing normalized coordinates and the scaling factor
 */
export function normalizeCoordinates(coordinates: number[][]): {
  normalized: number[][];
  scalingFactor: number;
} | null {
  if (!coordinates || coordinates.length === 0) {
    return null;
  }

  const n = coordinates[0].length;
  const maxAbs = coordinates.reduce((max, coord) => {
    return Math.max(max, ...coord.slice(0, n).map((c) => Math.abs(c)));
  }, 0);

  if (maxAbs === 0) {
    return {
      normalized: coordinates,
      scalingFactor: 1,
    };
  }

  const normalized = coordinates.map((coord) =>
    coord.slice(0, n).map((c) => c / maxAbs),
  );

  return {
    normalized,
    scalingFactor: maxAbs,
  };
}
