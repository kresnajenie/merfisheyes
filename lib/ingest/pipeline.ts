/**
 * Which pipeline the server's worker runs on an upload.
 *
 * "prepare" is scripts/prepare.py, the one a person runs on their own machine
 * (docs/plans/ONE-PIPELINE.md): one SpatialData store per upload, cells and
 * molecules together. Unset, the worker runs the legacy chunking processors,
 * exactly as before. Read by the upload UI (one dataset, not a cells dataset
 * plus a molecules dataset) and stamped on the dataset's processingParams,
 * which is what the worker, finalize and reconcile go by — so a dataset keeps
 * the pipeline it was uploaded under whatever the setting says later.
 */
export const PREPARE_PIPELINE = "prepare";

export function preparePipelineEnabled(): boolean {
  return process.env.NEXT_PUBLIC_INGEST_PIPELINE === PREPARE_PIPELINE;
}

export function usesPreparePipeline(processingParams: unknown): boolean {
  return (
    (processingParams as { pipeline?: unknown } | null)?.pipeline ===
    PREPARE_PIPELINE
  );
}

/** The counts prepare.py writes to OUT/stats.json, as a dataset row holds them. */
export function statsFromPrepare(
  stats: {
    cells?: number;
    genes?: number;
    molecules?: number;
    molecule_genes?: number;
  } | null,
): { numCells?: number; numGenes?: number } {
  return {
    numCells: stats?.cells || stats?.molecules || undefined,
    numGenes: stats?.genes || stats?.molecule_genes || undefined,
  };
}
