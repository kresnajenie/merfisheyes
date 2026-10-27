import type { PointsElement } from "./points";

import { SingleMoleculeDataset } from "@/lib/SingleMoleculeDataset";

/**
 * A single molecule dataset standing for a points element in the Zarr
 * points format: the genes and their counts, for the gene picker and the
 * legends. Its molecules aren't loaded here — the scene's PointsLayer reads
 * them for the view (`metadata.pointsElement`).
 */
export async function pointsMoleculeDataset(
  element: PointsElement,
  id: string,
  name: string,
): Promise<SingleMoleculeDataset> {
  const counts = await element.featureCounts();
  const moleculeCounts: Record<string, { assigned: number }> = {};

  element.features.forEach((f, i) => {
    moleculeCounts[f] = { assigned: counts[i] };
  });
  const ds = new SingleMoleculeDataset({
    id,
    name,
    type: "spatialdata_points",
    uniqueGenes: element.features,
    geneIndex: new Map(),
    dimensions: element.dims,
    scalingFactor: 1,
    moleculeCounts,
    // Molecules in no cell are marked per row, so the legends and picker
    // offer their switches; the PointsLayer draws them apart
    hasUnassigned: !!element.readUnassigned,
    metadata: { pointsElement: element },
  });

  ds.hasAssignedMolecules = (gene: string) =>
    (moleculeCounts[gene]?.assigned ?? 0) > 0;
  // Drawn by the PointsLayer, a view at a time; nothing to hand out whole.
  (ds as any).getCoordinatesByGene = async () => new Float32Array(0);

  return ds;
}

/** The points element behind a dataset made by pointsMoleculeDataset. */
export function pointsElementOf(
  ds: SingleMoleculeDataset | null,
): PointsElement | null {
  return (ds?.metadata?.pointsElement as PointsElement | undefined) ?? null;
}
