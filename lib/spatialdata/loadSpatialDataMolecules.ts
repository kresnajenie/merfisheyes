import type { SpatialDataStore } from "./store";
import type { DatasetLinkConfig } from "@/lib/config/dataset-links";
import type { SingleMoleculeDataset } from "@/lib/SingleMoleculeDataset";

/**
 * The molecules of a SpatialData store that has no cells to draw them over:
 * the single molecule dataset linked by the mapping.json beside the store
 * (its points, converted with scripts/spatialdata_points_to_sm.py). Null when
 * the store has none. With cells, the scene finds the same link itself.
 */
export async function loadSpatialDataMolecules({
  sidecar,
}: SpatialDataStore): Promise<SingleMoleculeDataset | null> {
  const bytes = await sidecar.get("/mapping.json");

  if (!bytes) return null;
  const mapping: Partial<DatasetLinkConfig> = JSON.parse(
    new TextDecoder().decode(bytes),
  );
  const molecules = mapping.links?.__all__;

  if (!molecules) return null;
  const { SingleMoleculeDataset } = await import("@/lib/SingleMoleculeDataset");

  return mapping.source === "app"
    ? SingleMoleculeDataset.fromS3(molecules)
    : SingleMoleculeDataset.fromCustomS3(molecules);
}
