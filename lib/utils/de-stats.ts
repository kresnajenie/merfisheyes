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

/**
 * Parse the gzip-decompressed binary layout written by Python's
 * write_de_stats_binary (scripts/process_spatial_data.py).
 *
 *   Header (16 bytes):
 *     version       uint32
 *     num_genes     uint32
 *     num_celltypes uint32
 *     reserved      uint32
 *   Celltype names (variable):
 *     for each celltype: { name_len: uint32; name_bytes: utf-8 }
 *   Cell counts:        uint32[C]
 *   Means:              float32[G * C]   row-major gene-major (g * C + c)
 *   Pct expressing:     float32[G * C]
 *
 * Gene order is positional and must match `genes` here (which mirrors
 * expr/index.json on disk).
 */
export function parseDeStatsBuffer(
  buffer: ArrayBuffer,
  column: string,
  genes: string[],
): DeStats {
  const view = new DataView(buffer);
  const version = view.getUint32(0, true);
  if (version !== 1) {
    throw new Error(`Unsupported de-stats binary version: ${version}`);
  }
  const numGenes = view.getUint32(4, true);
  const numCelltypes = view.getUint32(8, true);

  if (numGenes !== genes.length) {
    throw new Error(
      `de-stats numGenes ${numGenes} != dataset.genes.length ${genes.length} for column "${column}"`,
    );
  }

  let offset = 16;
  const decoder = new TextDecoder("utf-8");
  const celltypes: string[] = new Array(numCelltypes);
  for (let c = 0; c < numCelltypes; c++) {
    const nameLen = view.getUint32(offset, true);
    offset += 4;
    celltypes[c] = decoder.decode(
      new Uint8Array(buffer, offset, nameLen),
    );
    offset += nameLen;
  }

  // Cell counts: uint32[C]. Aligned reads via DataView (avoid alignment
  // assumptions of typed-array views over the source buffer).
  const cellCounts: number[] = new Array(numCelltypes);
  for (let c = 0; c < numCelltypes; c++) {
    cellCounts[c] = view.getUint32(offset, true);
    offset += 4;
  }

  const total = numGenes * numCelltypes;
  const meansBytes = total * 4;

  // Float32 views need 4-byte alignment relative to the source buffer's
  // origin. The buffer may not be aligned here, so copy into fresh arrays.
  const means = new Float32Array(total);
  const meansDV = new DataView(buffer, offset, meansBytes);
  for (let i = 0; i < total; i++) {
    means[i] = meansDV.getFloat32(i * 4, true);
  }
  offset += meansBytes;

  const pctExpressing = new Float32Array(total);
  const pctDV = new DataView(buffer, offset, meansBytes);
  for (let i = 0; i < total; i++) {
    pctExpressing[i] = pctDV.getFloat32(i * 4, true);
  }

  return {
    column,
    celltypes,
    cellCounts,
    genes,
    means,
    pctExpressing,
  };
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
 * Rank genes by fold change of mean expression for a target celltype.
 *
 * Reference is either:
 *   - All other cells (default): `mean_out` reconstructed from the rest of
 *     the per-celltype table via `Σ_c means[g,c] × n_c` minus the target.
 *   - A specific celltype (opts.reference): read directly from the table.
 *
 *   log2FC = log2((mean_in + ε) / (mean_out + ε))
 *
 * Returns genes sorted by log2FC desc.
 */
export function rankDegsForCelltype(
  deStats: DeStats,
  targetCelltype: string,
  opts: { pseudocount?: number; reference?: string | null } = {},
): RankedDeg[] {
  const pseudo = opts.pseudocount ?? 1e-9;
  const t = deStats.celltypes.indexOf(targetCelltype);
  if (t === -1) return [];

  const C = deStats.celltypes.length;
  const G = deStats.genes.length;
  const nIn = deStats.cellCounts[t];

  let r = -1;
  let nOut = 0;
  if (opts.reference) {
    r = deStats.celltypes.indexOf(opts.reference);
    if (r === -1 || r === t) return [];
    nOut = deStats.cellCounts[r];
  } else {
    const nTotal = deStats.cellCounts.reduce((a, b) => a + b, 0);
    nOut = nTotal - nIn;
  }
  if (nIn === 0 || nOut === 0) return [];

  const results: RankedDeg[] = new Array(G);
  const LOG2 = Math.log(2);
  const useSpecificRef = r !== -1;

  for (let g = 0; g < G; g++) {
    const base = g * C;
    const meanIn = deStats.means[base + t];
    const pctIn = deStats.pctExpressing[base + t];

    let meanOut: number;
    let pctOut: number;

    if (useSpecificRef) {
      meanOut = deStats.means[base + r];
      pctOut = deStats.pctExpressing[base + r];
    } else {
      let sumMeanWeighted = 0;
      let sumPctWeighted = 0;
      for (let c = 0; c < C; c++) {
        const nC = deStats.cellCounts[c];
        sumMeanWeighted += deStats.means[base + c] * nC;
        sumPctWeighted += deStats.pctExpressing[base + c] * nC;
      }
      meanOut = (sumMeanWeighted - meanIn * nIn) / nOut;
      pctOut = (sumPctWeighted - pctIn * nIn) / nOut;
    }

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
 * Minimal dataset shape we depend on for on-demand recompute or fetch.
 * Defined locally to avoid a circular import with StandardizedDataset.
 */
interface DatasetForDeStats {
  id: string;
  genes: string[];
  matrix: any;
  clusters: ClusterLike[] | null;
  deStatsByColumn?: Map<string, DeStats>;
  availableDeStatsColumns?: string[];
  adapter?: {
    loadDeStats?: (column: string, genes: string[]) => Promise<DeStats | null>;
  } | null;
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

  // Decide between fetch-from-disk (precomputed by Python) and main-thread
  // recompute. Fetch is preferred when the adapter advertises this column.
  const canFetch =
    !!dataset.adapter?.loadDeStats &&
    !!dataset.availableDeStatsColumns?.includes(column);

  if (!canFetch) {
    const cluster = dataset.clusters?.find((c) => c.column === column);
    if (!cluster || cluster.type !== "categorical") return null;
    if (!cluster.valueIndices || !cluster.uniqueValues) return null;
    if (!dataset.matrix) return null;
  }

  deStatsInFlight.add(key);
  const compute = async (): Promise<DeStats | null> => {
    if (canFetch) {
      try {
        const fetched = await dataset.adapter!.loadDeStats!(
          column,
          dataset.genes,
        );
        if (fetched) {
          await onProgress?.(1);
          return fetched;
        }
      } catch (err) {
        console.warn(
          `[deStats] adapter.loadDeStats("${column}") failed; falling back to compute:`,
          err,
        );
      }
    }
    const cluster = dataset.clusters?.find((c) => c.column === column);
    if (!cluster || cluster.type !== "categorical") return null;
    if (!cluster.valueIndices || !cluster.uniqueValues) return null;
    if (!dataset.matrix) return null;
    return computeDeStats(cluster, dataset.genes, dataset.matrix, {
      onProgress,
    });
  };

  const promise = compute()
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
