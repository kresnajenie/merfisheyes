import {
  COMPUTE_TIER_ORDER,
  COMPUTE_TIERS,
  normalizeComputeTier,
  type ComputeTier,
} from "@/lib/batch";

const GB = 1024 ** 3;

/**
 * Peak memory of prepare.py per GB uploaded, by file, from the runs in
 * docs/plans/ONE-PIPELINE.md §1 (each measured in a container):
 *
 * - a molecule parquet: Beaumont set1 6.1 GB → 15.6 GB, set3 3.4 → 12.8,
 *   set4 2.3 → 8.5
 * - a molecule CSV: MERSCOPE uterine 32 GB → 10.8 GB, colon 20 GB (with
 *   704k cells) → 12.5 GB
 * - a cell table (MERSCOPE's cell_by_gene / cell_metadata CSVs, .h5ad,
 *   Xenium's matrix): a few times its size once dense in memory
 * - an image: read a tile at a time
 */
function peakPerGb(key: string): number {
  const name = key.toLowerCase().split("/").pop() ?? "";

  if (name.endsWith(".parquet")) return 3.5;
  if (/^(cell_by_gene|cell_metadata|cells)\./.test(name)) return 4;
  if (name.endsWith(".csv")) return 0.3;
  if (name.endsWith(".csv.gz")) return 1.5;
  if (/\.(ome\.)?tiff?$/.test(name) || key.includes("/images/")) return 0.5;

  return 3;
}

/** What the run is expected to need, in GB: the files plus the interpreter. */
export function estimatePeakGb(files: { key: string; size: number }[]): number {
  return files.reduce((sum, f) => sum + (f.size / GB) * peakPerGb(f.key), 1.5);
}

/**
 * The compute tier for a prepare.py run, chosen up front from what was
 * uploaded (decision f): the smallest tier whose memory the estimate stays
 * under 75% of, never below the owner's own tier. A run that still runs out
 * of memory is retried a tier up, as before.
 */
export function tierForUpload(
  files: { key: string; size: number }[],
  floor?: unknown,
): ComputeTier {
  const need = estimatePeakGb(files);
  const start = COMPUTE_TIER_ORDER.indexOf(normalizeComputeTier(floor));

  for (const tier of COMPUTE_TIER_ORDER.slice(start)) {
    if (need <= (0.75 * Number(COMPUTE_TIERS[tier].memoryMiB)) / 1024) {
      return tier;
    }
  }

  return COMPUTE_TIER_ORDER[COMPUTE_TIER_ORDER.length - 1];
}
