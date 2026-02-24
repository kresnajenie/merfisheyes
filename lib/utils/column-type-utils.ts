import type { StandardizedDataset } from "@/lib/StandardizedDataset";

/**
 * Resolves the effective column type, accounting for user overrides.
 * Priority: override > loaded cluster type > allClusterColumnTypes > "categorical" default.
 */
export function getEffectiveColumnType(
  column: string,
  dataset: StandardizedDataset,
  overrides: Record<string, "categorical" | "numerical">,
): "categorical" | "numerical" {
  if (overrides[column]) return overrides[column];

  const loadedCluster = dataset.clusters?.find((c) => c.column === column);

  if (loadedCluster) return loadedCluster.type as "categorical" | "numerical";

  if (dataset.allClusterColumnTypes?.[column]) {
    return dataset.allClusterColumnTypes[column] as
      | "categorical"
      | "numerical";
  }

  return "categorical";
}

/**
 * Checks whether a column can be interpreted as numerical.
 * Returns true if the auto-detected type is already "numerical",
 * or if the loaded values are all parseable as finite numbers.
 * For unloaded columns, falls back to allClusterColumnTypes.
 */
export function canTreatAsNumerical(
  column: string,
  dataset: StandardizedDataset,
): boolean {
  // If the auto-detected type is already numerical, yes
  if (dataset.allClusterColumnTypes?.[column] === "numerical") return true;

  const loadedCluster = dataset.clusters?.find((c) => c.column === column);

  if (!loadedCluster) {
    // Not loaded — can't inspect values, fall back to stored type
    return dataset.allClusterColumnTypes?.[column] === "numerical";
  }

  if (loadedCluster.type === "numerical") return true;

  // Sample up to 50 values to check if they parse as numbers
  const values = loadedCluster.values;
  const sampleSize = Math.min(values.length, 50);

  for (let i = 0; i < sampleSize; i++) {
    const v = values[i];

    if (v === null || v === undefined || v === "") continue;
    if (typeof v === "number") continue;
    const n = Number(v);

    if (!Number.isFinite(n)) return false;
  }

  return true;
}
