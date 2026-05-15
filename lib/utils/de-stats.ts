/**
 * Per-celltype expression stats precomputed during in-browser parsing.
 *
 * For every gene × every celltype value of the priority categorical cluster
 * column, we store:
 *   - mean expression
 *   - pct expressing (fraction of cells with nonzero expression)
 *
 * Plus per-celltype cell counts so consumers can reconstruct global per-gene
 * sums (means × counts) and derive fold-change-vs-rest later.
 *
 * Values are computed on whatever the adapter stored in the matrix — raw or
 * normalized — no transformation here.
 */
export interface DeStats {
  column: string;
  celltypes: string[];
  cellCounts: number[];
  genes: string[];
  means: Float32Array;
  pctExpressing: Float32Array;
}

interface ClusterLike {
  column: string;
  type: string;
  valueIndices?: Uint16Array | Uint32Array;
  uniqueValues?: string[];
}

const PRIORITY_COLUMNS = [
  "class_name",
  "subclass_name",
  "supertype_name",
  "cluster_name",
  "subcluster_name",
  "super_cluster_name",
  "leiden",
  "Cluster",
];

/**
 * Pick the priority categorical column. Walks the same priority order as
 * selectBestClusterColumn() but skips any column that isn't categorical.
 * Returns null if no categorical column exists.
 */
export function selectPriorityCategoricalCluster<T extends ClusterLike>(
  clusters: T[] | null | undefined,
): T | null {
  if (!clusters || clusters.length === 0) return null;

  for (const name of PRIORITY_COLUMNS) {
    const hit = clusters.find(
      (c) => c.column === name && c.type === "categorical",
    );
    if (hit) return hit;
  }

  const celltype = clusters.find(
    (c) =>
      c.type === "categorical" && c.column.toLowerCase().includes("celltype"),
  );
  if (celltype) return celltype;

  const anyCat = clusters.find((c) => c.type === "categorical");
  if (anyCat) return anyCat;

  return null;
}

export type DeStatsProgress = (fraction: number) => void | Promise<void>;

interface ComputeOpts {
  onProgress?: DeStatsProgress;
  yieldEvery?: number; // rows between event-loop yields
}

const DEFAULT_YIELD_EVERY = 100_000;

/**
 * Compute per-celltype mean expression and pct expressing for every gene.
 * Single linear pass over the matrix, yielding to the event loop every
 * `yieldEvery` cells so the UI stays responsive on multi-million-cell datasets.
 * Returns null if inputs are incomplete.
 *
 * Supported matrix shapes:
 *   - Flat TypedArray, row-major [cell × gene]: index = cellIdx * G + geneIdx
 *   - number[][] nested row-major: matrix[cellIdx][geneIdx]
 *
 * Map<gene, Float32Array> (Xenium/MERSCOPE) is not handled yet — H5AD only.
 */
export async function computeDeStats(
  cluster: ClusterLike,
  genes: string[],
  matrix: any,
  opts: ComputeOpts = {},
): Promise<DeStats | null> {
  if (!cluster || !matrix || !genes || genes.length === 0) return null;
  if (!cluster.valueIndices || !cluster.uniqueValues) return null;

  const celltypes = cluster.uniqueValues;
  const C = celltypes.length;
  const G = genes.length;
  if (C === 0) return null;

  const valueIndices = cluster.valueIndices;
  const N = valueIndices.length;
  const yieldEvery = opts.yieldEvery ?? DEFAULT_YIELD_EVERY;

  const cellCounts = new Array<number>(C).fill(0);
  for (let i = 0; i < N; i++) {
    cellCounts[valueIndices[i]]++;
  }

  // Accumulate in Float64 to keep precision over millions of cells; cast at end.
  const meanSum = new Float64Array(G * C);
  const nonzeroCounts = new Uint32Array(G * C);

  if (ArrayBuffer.isView(matrix)) {
    const m = matrix as unknown as { [i: number]: number; length: number };
    const expectedLen = N * G;
    if (m.length < expectedLen) {
      console.warn(
        `[deStats] matrix length ${m.length} < expected ${expectedLen} (N=${N} × G=${G}); skipping`,
      );
      return null;
    }
    for (let i = 0; i < N; i++) {
      const ct = valueIndices[i];
      const rowOffset = i * G;
      for (let g = 0; g < G; g++) {
        const v = m[rowOffset + g];
        if (!v) continue;
        const idx = g * C + ct;
        meanSum[idx] += v;
        nonzeroCounts[idx]++;
      }
      if (i > 0 && i % yieldEvery === 0) {
        await opts.onProgress?.(i / N);
        await yieldToEventLoop();
      }
    }
  } else if (Array.isArray(matrix) && Array.isArray(matrix[0])) {
    for (let i = 0; i < N; i++) {
      const ct = valueIndices[i];
      const row = matrix[i];
      if (!row) continue;
      for (let g = 0; g < G; g++) {
        const v = row[g];
        if (!v) continue;
        const idx = g * C + ct;
        meanSum[idx] += v;
        nonzeroCounts[idx]++;
      }
      if (i > 0 && i % yieldEvery === 0) {
        await opts.onProgress?.(i / N);
        await yieldToEventLoop();
      }
    }
  } else {
    console.warn("[deStats] Unsupported matrix shape; skipping");
    return null;
  }

  const means = new Float32Array(G * C);
  const pctExpressing = new Float32Array(G * C);
  for (let g = 0; g < G; g++) {
    for (let ct = 0; ct < C; ct++) {
      const idx = g * C + ct;
      const cnt = cellCounts[ct];
      if (cnt > 0) {
        means[idx] = meanSum[idx] / cnt;
        pctExpressing[idx] = nonzeroCounts[idx] / cnt;
      }
    }
  }

  await opts.onProgress?.(1);

  return {
    column: cluster.column,
    celltypes,
    cellCounts,
    genes,
    means,
    pctExpressing,
  };
}

function yieldToEventLoop(): Promise<void> {
  return new Promise((resolve) => setTimeout(resolve, 0));
}

export interface RankedDeg {
  gene: string;
  log2FC: number;
  meanIn: number;
  meanOut: number;
  pctIn: number;
  pctOut: number;
}

/**
 * Rank genes by fold change of mean expression for a target celltype vs the
 * rest. Reconstructs per-gene global sums from `means × cellCounts`, so we
 * never re-iterate the matrix.
 *
 *   mean_out = (Σ_c means[g,c] × n_c  −  means[g,t] × n_t) / (N − n_t)
 *   log2FC   = log2((mean_in + ε) / (mean_out + ε))
 *
 * Returns genes sorted by log2FC desc.
 */
export function rankDegsForCelltype(
  deStats: DeStats,
  targetCelltype: string,
  opts: { pseudocount?: number } = {},
): RankedDeg[] {
  const pseudo = opts.pseudocount ?? 1e-9;
  const t = deStats.celltypes.indexOf(targetCelltype);
  if (t === -1) return [];

  const C = deStats.celltypes.length;
  const G = deStats.genes.length;
  const nIn = deStats.cellCounts[t];
  const nTotal = deStats.cellCounts.reduce((a, b) => a + b, 0);
  const nOut = nTotal - nIn;
  if (nIn === 0 || nOut === 0) return [];

  const results: RankedDeg[] = new Array(G);
  const LOG2 = Math.log(2);

  for (let g = 0; g < G; g++) {
    const base = g * C;
    let sumMeanWeighted = 0;
    let sumPctWeighted = 0;
    for (let c = 0; c < C; c++) {
      const nC = deStats.cellCounts[c];
      sumMeanWeighted += deStats.means[base + c] * nC;
      sumPctWeighted += deStats.pctExpressing[base + c] * nC;
    }
    const meanIn = deStats.means[base + t];
    const pctIn = deStats.pctExpressing[base + t];
    const meanOut = (sumMeanWeighted - meanIn * nIn) / nOut;
    const pctOut = (sumPctWeighted - pctIn * nIn) / nOut;

    const log2FC = Math.log((meanIn + pseudo) / (meanOut + pseudo)) / LOG2;

    results[g] = {
      gene: deStats.genes[g],
      log2FC,
      meanIn,
      meanOut,
      pctIn,
      pctOut,
    };
  }

  results.sort((a, b) => b.log2FC - a.log2FC);
  return results;
}

/**
 * Minimal dataset shape we depend on for on-demand recompute. Defined locally
 * to avoid a circular import with StandardizedDataset.
 */
interface DatasetForDeStats {
  id: string;
  genes: string[];
  matrix: any;
  clusters: ClusterLike[] | null;
  deStatsByColumn?: Map<string, DeStats>;
}

const inFlight = new Map<string, Promise<DeStats | null>>();

/**
 * Track which (dataset, column) pairs are currently computing. Components can
 * read this synchronously to render a loading state. The Set lives at module
 * scope and is mutated in place by ensureDeStatsForColumn.
 */
export const deStatsInFlight = new Set<string>();

function flightKey(datasetId: string, column: string): string {
  return `${datasetId}::${column}`;
}

/**
 * Return the deStats for a given cluster column on this dataset, computing it
 * once if needed. Cached on `dataset.deStatsByColumn`. Concurrent calls for
 * the same (dataset, column) are deduplicated.
 *
 *  - Returns null and skips compute if the column isn't loaded as a
 *    categorical cluster on the dataset.
 *  - Returns null on unsupported matrix shapes (e.g. Map-based for
 *    Xenium/MERSCOPE — not yet wired up).
 */
export async function ensureDeStatsForColumn(
  dataset: DatasetForDeStats,
  column: string,
  onProgress?: DeStatsProgress,
): Promise<DeStats | null> {
  if (!dataset || !column) return null;

  const cache = dataset.deStatsByColumn;
  if (cache?.has(column)) return cache.get(column) ?? null;

  const key = flightKey(dataset.id, column);
  const existing = inFlight.get(key);
  if (existing) return existing;

  const cluster = dataset.clusters?.find((c) => c.column === column);
  if (!cluster || cluster.type !== "categorical") return null;
  if (!cluster.valueIndices || !cluster.uniqueValues) return null;
  if (!dataset.matrix) return null;

  deStatsInFlight.add(key);
  const promise = computeDeStats(cluster, dataset.genes, dataset.matrix, {
    onProgress,
  })
    .then((result) => {
      if (result && dataset.deStatsByColumn) {
        dataset.deStatsByColumn.set(column, result);
      }
      return result;
    })
    .finally(() => {
      inFlight.delete(key);
      deStatsInFlight.delete(key);
    });

  inFlight.set(key, promise);
  return promise;
}

export function isDeStatsInFlight(datasetId: string, column: string): boolean {
  return deStatsInFlight.has(flightKey(datasetId, column));
}
