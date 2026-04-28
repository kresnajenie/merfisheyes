import type { StandardizedDataset } from "@/lib/StandardizedDataset";

/**
 * Lazy-load a cluster column on demand. No-ops if the column is already
 * present on the dataset. Mutates `dataset.clusters` via `addClusters`.
 *
 * @returns true if a fetch occurred (caller may want to bump clusterVersion).
 */
export async function loadClusterColumn(
  dataset: StandardizedDataset,
  column: string,
): Promise<boolean> {
  if (!dataset.adapter) return false;
  if (dataset.clusters?.some((c) => c.column === column)) return false;

  let newClusters: Array<{
    column: string;
    type: string;
    values: any[];
    palette: Record<string, string> | null;
    uniqueValues?: string[];
  }> | null = null;

  if (dataset.adapter.mode === "local") {
    newClusters = await dataset.adapter.loadClusters([column]);
  } else {
    const { getStandardizedDatasetWorker } = await import(
      "@/lib/workers/standardizedDatasetWorkerManager"
    );
    const worker = await getStandardizedDatasetWorker();
    newClusters = await worker.loadClusterFromS3(
      dataset.id,
      [column],
      dataset.metadata?.customS3BaseUrl,
    );
  }

  if (newClusters && newClusters.length > 0) {
    dataset.addClusters(newClusters);
    return true;
  }
  return false;
}
