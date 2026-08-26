// lib/workers/standardized-dataset.worker.ts
import * as Comlink from "comlink";

import { H5adAdapter } from "../adapters/H5adAdapter";
import { XeniumAdapter } from "../adapters/XeniumAdapter";
import { MerscopeAdapter } from "../adapters/MerscopeAdapter";
import { ChunkedDataAdapter } from "../adapters/ChunkedDataAdapter";
import { normalizeCoordinates, normalizeCoordinatesFlat, round2 } from "../utils/coordinates";
import { shouldFilterGene } from "../utils/gene-filters";
import { isCategorical } from "../utils/column-type-detection";
import { selectBestClusterColumnByName } from "../utils/dataset-utils";
import {
  computeDeStats,
  selectPriorityCategoricalCluster,
  type DeStats,
} from "../utils/de-stats";

/**
 * Progress callback type that can be proxied by Comlink
 */
type ProgressCallback = (
  progress: number,
  message: string,
) => void | Promise<void>;

/**
 * Determine if data is categorical or numerical
 * Now uses centralized detection logic from column-type-detection utility
 */
function isCategoricalData(values: any[], columnName?: string): boolean {
  return isCategorical(values, columnName);
}

/**
 * Serializable cluster data
 */
interface SerializableClusterData {
  column: string;
  type: string;
  values: any[];
  valueIndices?: Uint16Array | Uint32Array;
  palette: Record<string, string> | null;
  uniqueValues?: string[];
}

/**
 * Serializable spatial data
 */
interface SerializableSpatialData {
  coordinates: Float32Array | number[][];
  dimensions: number;
  scalingFactor: number;
}

/**
 * Serializable dataset data that can be transferred to main thread
 */
export interface SerializableStandardizedDataset {
  id: string;
  name: string;
  type: string;
  spatial: SerializableSpatialData;
  embeddings: Record<string, number[][]>;
  genes: string[];
  clusters: SerializableClusterData[] | null;
  metadata: Record<string, any>;
  matrix: any; // Pre-loaded expression matrix for gene visualization
  allClusterColumnNames?: string[];
  allClusterColumnTypes?: Record<string, string>;
  allEmbeddingNames?: string[];
  normalized?: boolean; // false = raw coordinates, true/undefined = normalized [-1,1]
  deStats?: DeStats | null;
  availableDeStatsColumns?: string[];
}

/**
 * Worker API exposed via Comlink
 */
// Exported so the parity harness (scripts/parity/run-js.mts) can run the exact
// browser parsing path headlessly in Node; the browser only ever reaches it
// through Comlink.
export const workerApi = {
  /**
   * Parse H5AD file and return serializable dataset data
   */
  async parseH5ad(
    file: File,
    onProgress?: ProgressCallback,
  ): Promise<SerializableStandardizedDataset> {
    console.log("[Worker] Starting H5AD parsing:", file.name);

    const adapter = new H5adAdapter();

    await adapter.initialize(file, onProgress);

    // Load all data through adapter
    await onProgress?.(92, "Loading spatial coordinates...");
    const spatial = adapter.loadSpatialCoordinates();

    console.log("[Worker] Spatial data:", spatial);
    await onProgress?.(94, "Loading embeddings...");
    const embeddings = adapter.loadEmbeddings();

    // Embeddings are rounded to 2 dp like the spatial coordinates (the server
    // pipeline rounds every coordinate file the same way).
    for (const [name, coords] of Object.entries(embeddings)) {
      embeddings[name] = coords.map((point) => point.map(round2));
    }

    console.log("[Worker] Embeddings:", embeddings);
    await onProgress?.(96, "Loading genes...");
    const genes = await adapter.loadGenes();

    console.log("[Worker] Genes:", genes.length, "genes loaded");
    await onProgress?.(98, "Loading clusters...");
    const clusters = await adapter.loadClusters();

    console.log("[Worker] Clusters:", clusters);
    const dataInfo = adapter.getDatasetInfo();

    console.log("[Worker] Dataset info:", dataInfo);

    // Load expression matrix for gene visualization
    await onProgress?.(99, "Loading expression matrix...");
    let matrix = adapter.fetchFullMatrix();

    console.log("[Worker] Expression matrix loaded");

    // Control probes / blanks are dropped in every format (shared rule with
    // the server pipeline). The dense row-major matrix is compacted in place —
    // kept columns only move left — so no second copy is needed.
    const keptIdx = genes.map((g, i) => (shouldFilterGene(g) ? -1 : i)).filter((i) => i >= 0);

    if (keptIdx.length !== genes.length) {
      const numCellsTotal = dataInfo.numCells;
      const src = genes.length;
      const dst = keptIdx.length;

      if (ArrayBuffer.isView(matrix)) {
        const flat = matrix as any;

        for (let r = 0; r < numCellsTotal; r++) {
          for (let k = 0; k < dst; k++) flat[r * dst + k] = flat[r * src + keptIdx[k]];
        }
        matrix = flat.subarray(0, numCellsTotal * dst);
      } else if (matrix instanceof Map) {
        for (const g of genes) if (shouldFilterGene(g)) matrix.delete(g);
      } else if (Array.isArray(matrix)) {
        matrix = matrix.map((row: any[]) => keptIdx.map((i) => row[i]));
      }
      console.log(`[Worker] Filtered ${src - dst} control/blank genes`);
      genes.splice(0, genes.length, ...keptIdx.map((i) => genes[i]));
      dataInfo.numGenes = genes.length;
    }

    // Precompute per-celltype expression stats for the priority categorical
    // cluster column. Skipped gracefully if no categorical column exists.
    const priorityCluster = selectPriorityCategoricalCluster(clusters);
    const deStats = priorityCluster
      ? await computeDeStats(priorityCluster, genes, matrix)
      : null;

    if (deStats) {
      console.log(
        `[Worker] deStats computed for "${deStats.column}": ${deStats.celltypes.length} celltypes × ${deStats.genes.length} genes`,
      );
    } else {
      console.log("[Worker] deStats skipped (no categorical cluster column)");
    }

    // Round spatial coordinates to 2 decimal places (no normalization)
    const coords = spatial.coordinates;
    const roundedCoords: number[][] = [];

    for (let i = 0; i < coords.length; i++) {
      const point = coords[i];

      roundedCoords.push(
        point.map(round2),
      );
    }

    // Generate dataset ID
    const timestamp = Date.now();
    const random = Math.random().toString(36).substring(2, 11);
    const id = `h5ad_${file.name.replace(".h5ad", "")}_${timestamp}_${random}`;

    const serializable: SerializableStandardizedDataset = {
      id: id,
      name: file.name.replace(".h5ad", ""),
      type: "h5ad",
      spatial: {
        coordinates: roundedCoords,
        dimensions: spatial.dimensions,
        scalingFactor: 1,
      },
      embeddings: embeddings,
      genes: genes,
      clusters: clusters,
      metadata: {
        ...adapter.metadata,
        originalFileName: file.name,
        numCells: dataInfo.numCells,
        numGenes: dataInfo.numGenes,
        spatialScalingFactor: 1,
      },
      matrix: matrix,
      normalized: false,
      deStats: deStats,
    };

    console.log("[Worker] H5AD parsing complete");

    return serializable;
  },

  /**
   * Parse Xenium files and return serializable dataset data
   */
  async parseXenium(
    files: File[],
    onProgress?: ProgressCallback,
  ): Promise<SerializableStandardizedDataset> {
    console.log("[Worker] Starting Xenium parsing:", files.length, "files");

    const adapter = new XeniumAdapter();

    await adapter.initialize(files, onProgress);

    // Load all data through adapter
    await onProgress?.(92, "Loading spatial coordinates...");
    const spatial = adapter.loadSpatialCoordinates();

    console.log("[Worker] Spatial data:", spatial);
    await onProgress?.(94, "Loading embeddings...");
    const embeddings = adapter.loadEmbeddings();

    console.log("[Worker] Embeddings:", embeddings);
    await onProgress?.(96, "Loading genes...");
    const genes = await adapter.loadGenes();

    console.log("[Worker] Genes:", genes.length, "genes loaded");
    await onProgress?.(98, "Loading clusters...");
    // Every metadata column (same as the server pipeline), typed per column.
    const clusterList = await adapter.loadAllClusters();

    console.log("[Worker] Clusters:", clusterList.map((c) => c.column));

    const clusters = clusterList.length
      ? clusterList.map((clusterData) => {
          // Detect on per-cell values: a unique-values list always looks
          // "numerical" (unique ratio 1). "" is a missing label, like NaN.
          const detectValues =
            clusterData.valueIndices && clusterData.uniqueValues
              ? Array.from(clusterData.valueIndices, (i) => clusterData.uniqueValues![i] || null)
              : clusterData.values ?? [];
          const categorical = isCategoricalData(detectValues, clusterData.column);

          return {
            column: clusterData.column,
            type: categorical ? "categorical" : "numerical",
            values: clusterData.values,
            valueIndices: clusterData.valueIndices,
            palette: categorical ? clusterData.palette : null,
            uniqueValues: clusterData.uniqueValues,
          };
        })
      : null;

    const dataInfo = adapter.getDatasetInfo();

    console.log("[Worker] Dataset info:", dataInfo);

    // Load expression matrix for gene visualization
    await onProgress?.(99, "Loading expression matrix...");
    const matrix = adapter.fetchFullMatrix();

    console.log("[Worker] Expression matrix loaded");

    // Precompute per-celltype DE stats from the loaded matrix (mirrors the
    // H5AD branch). Skipped gracefully if there's no categorical column.
    const priorityCluster = selectPriorityCategoricalCluster(clusters);
    const deStats = priorityCluster
      ? await computeDeStats(priorityCluster, genes, matrix)
      : null;

    if (deStats) {
      console.log(
        `[Worker] deStats computed for "${deStats.column}": ${deStats.celltypes.length} celltypes × ${deStats.genes.length} genes`,
      );
    } else {
      console.log("[Worker] deStats skipped (no categorical cluster column)");
    }

    // Round spatial coordinates to 2 decimal places (no normalization)
    const coords = spatial.coordinates;
    const roundedCoords: number[][] = [];

    for (let i = 0; i < coords.length; i++) {
      const point = coords[i];

      roundedCoords.push(
        point.map(round2),
      );
    }

    // Generate dataset ID and name from folder
    const timestamp = Date.now();
    const random = Math.random().toString(36).substring(2, 11);
    const folderName =
      files[0]?.webkitRelativePath?.split("/")[0] || "xenium_data";
    const id = `xenium_${folderName}_${timestamp}_${random}`;

    const serializable: SerializableStandardizedDataset = {
      id: id,
      name: folderName,
      type: "xenium",
      spatial: {
        coordinates: roundedCoords,
        dimensions: spatial.dimensions,
        scalingFactor: 1,
      },
      embeddings: embeddings,
      genes: genes,
      clusters: clusters,
      metadata: {
        ...((adapter as any).metadata || {}),
        fileCount: files.length,
        numCells: dataInfo.numCells,
        numGenes: dataInfo.numGenes,
        spatialScalingFactor: 1,
      },
      matrix: matrix,
      normalized: false,
      deStats: deStats,
    };

    console.log("[Worker] Xenium parsing complete");

    return serializable;
  },

  /**
   * Parse MERSCOPE files and return serializable dataset data
   */
  async parseMerscope(
    files: File[],
    onProgress?: ProgressCallback,
  ): Promise<SerializableStandardizedDataset> {
    console.log("[Worker] Starting MERSCOPE parsing:", files.length, "files");

    const adapter = new MerscopeAdapter();

    await adapter.initialize(files, onProgress);

    // Load all data through adapter
    await onProgress?.(92, "Loading spatial coordinates...");
    const spatial = adapter.loadSpatialCoordinates();

    console.log("[Worker] Spatial data:", spatial);
    await onProgress?.(94, "Loading embeddings...");
    const rawEmbeddings = adapter.loadEmbeddings();

    console.log("[Worker] Raw embeddings:", rawEmbeddings);

    // Convert embeddings to proper format (filter out undefined values)
    const embeddings: Record<string, number[][]> = {};

    if (rawEmbeddings?.umap !== undefined) {
      embeddings.umap = rawEmbeddings.umap;
    }

    await onProgress?.(96, "Loading genes...");
    const genes = await adapter.loadGenes();

    console.log("[Worker] Genes:", genes.length, "genes loaded");
    await onProgress?.(98, "Loading clusters...");
    // Every metadata column (same as the server pipeline), typed per column.
    const clusterList = await adapter.loadAllClusters();

    console.log("[Worker] Clusters:", clusterList.map((c) => c.column));

    const clusters = clusterList.length
      ? clusterList.map((clusterData) => {
          // Detect on per-cell values: a unique-values list always looks
          // "numerical" (unique ratio 1). "" is a missing label, like NaN.
          const detectValues =
            clusterData.valueIndices && clusterData.uniqueValues
              ? Array.from(clusterData.valueIndices, (i) => clusterData.uniqueValues![i] || null)
              : clusterData.values ?? [];
          const categorical = isCategoricalData(detectValues, clusterData.column);

          return {
            column: clusterData.column,
            type: categorical ? "categorical" : "numerical",
            values: clusterData.values,
            valueIndices: clusterData.valueIndices,
            palette: categorical ? clusterData.palette : null,
            uniqueValues: clusterData.uniqueValues,
          };
        })
      : null;

    const dataInfo = adapter.getDatasetInfo();

    console.log("[Worker] Dataset info:", dataInfo);

    // Load expression matrix for gene visualization
    await onProgress?.(99, "Loading expression matrix...");
    const matrix = adapter.fetchFullMatrix();

    console.log("[Worker] Expression matrix loaded");

    // Precompute per-celltype DE stats from the loaded matrix (mirrors the
    // H5AD branch). Skipped gracefully if there's no categorical column.
    const priorityCluster = selectPriorityCategoricalCluster(clusters);
    const deStats = priorityCluster
      ? await computeDeStats(priorityCluster, genes, matrix)
      : null;

    if (deStats) {
      console.log(
        `[Worker] deStats computed for "${deStats.column}": ${deStats.celltypes.length} celltypes × ${deStats.genes.length} genes`,
      );
    } else {
      console.log("[Worker] deStats skipped (no categorical cluster column)");
    }

    // Round spatial coordinates to 2 decimal places (no normalization)
    const coords = spatial.coordinates;
    const roundedCoords: number[][] = [];

    for (let i = 0; i < coords.length; i++) {
      const point = coords[i];

      roundedCoords.push(
        point.map(round2),
      );
    }

    // Generate dataset ID and name from folder
    const timestamp = Date.now();
    const random = Math.random().toString(36).substring(2, 11);
    const folderName =
      files[0]?.webkitRelativePath?.split("/")[0] || "merscope_data";
    const id = `merscope_${folderName}_${timestamp}_${random}`;

    const serializable: SerializableStandardizedDataset = {
      id: id,
      name: folderName,
      type: "merscope",
      spatial: {
        coordinates: roundedCoords,
        dimensions: spatial.dimensions,
        scalingFactor: 1,
      },
      embeddings: embeddings,
      genes: genes,
      clusters: clusters,
      metadata: {
        ...((adapter as any).metadata || {}),
        fileCount: files.length,
        numCells: dataInfo.numCells,
        numGenes: dataInfo.numGenes,
        spatialScalingFactor: 1,
      },
      matrix: matrix,
      normalized: false,
      deStats: deStats,
    };

    console.log("[Worker] MERSCOPE parsing complete");

    return serializable;
  },

  /**
   * Load dataset from S3 and return serializable dataset data
   */
  async parseS3(
    datasetId: string,
    onProgress?: ProgressCallback,
    priorityColumnHint?: string,
  ): Promise<SerializableStandardizedDataset> {
    console.log("[Worker] Starting S3 loading:", datasetId);

    const adapter = new ChunkedDataAdapter(datasetId);

    await adapter.initialize();

    // Get cluster column info for deferred loading
    const clusterColumnInfo = adapter.getClusterColumnInfo();

    // Use URL-specified column if it exists in the dataset, otherwise auto-detect
    let priorityColumn: string | null = null;

    if (
      priorityColumnHint &&
      clusterColumnInfo.names.includes(priorityColumnHint)
    ) {
      priorityColumn = priorityColumnHint;
    } else {
      priorityColumn = selectBestClusterColumnByName(
        clusterColumnInfo.names,
        clusterColumnInfo.types,
      );
    }

    // Load all data through adapter (all methods return Promises)
    await onProgress?.(30, "Loading spatial coordinates...");
    const spatial = await adapter.loadSpatialCoordinates();

    console.log("[Worker] Spatial data:", spatial);

    // Skip eager embedding loading — embeddings are loaded on demand in the UI
    const embeddings: Record<string, number[][]> = {};

    await onProgress?.(60, "Loading genes...");
    const genes = await adapter.loadGenes();

    console.log("[Worker] Genes:", genes.length, "genes loaded");
    await onProgress?.(75, "Loading priority cluster column...");
    // Only load the priority cluster column (rest loaded in background on main thread)
    const clusters = priorityColumn
      ? await adapter.loadClusters([priorityColumn])
      : null;

    console.log("[Worker] Priority cluster loaded:", priorityColumn);
    const dataInfo = adapter.getDatasetInfo();

    console.log("[Worker] Dataset info:", dataInfo);

    // Load expression matrix for gene visualization
    await onProgress?.(90, "Loading expression matrix...");
    const matrix = await adapter.fetchFullMatrix();

    console.log("[Worker] Expression matrix loaded");

    // All datasets now use raw coordinates.
    // Old datasets (normalized !== false): coords are in [-1,1], denormalize using scalingFactor.
    // New datasets (normalized: false): coords are already raw microns.
    const wasNormalized = dataInfo.normalized !== false;
    const coords = spatial.coordinates;
    let finalCoords: Float32Array | number[][] = coords;
    let scalingFactor = (spatial as any).scalingFactor || 1;

    if (wasNormalized && scalingFactor > 1) {
      // Old dataset: denormalize back to raw coordinates
      console.log(`[Worker] Denormalizing old coords (×${scalingFactor})`);
      if (coords instanceof Float32Array) {
        const denorm = new Float32Array(coords.length);
        for (let i = 0; i < coords.length; i++) {
          denorm[i] = coords[i] * scalingFactor;
        }
        finalCoords = denorm;
      } else {
        finalCoords = (coords as number[][]).map((point) =>
          point.map((v) => v * scalingFactor),
        );
      }
    } else if (!wasNormalized) {
      console.log("[Worker] Using raw coordinates (normalized: false)");
    }

    const serializable: SerializableStandardizedDataset = {
      id: dataInfo.id,
      name: dataInfo.name,
      type: dataInfo.type,
      spatial: {
        coordinates: finalCoords,
        dimensions: spatial.dimensions,
        scalingFactor: 1, // always 1 now — coords are raw
      },
      embeddings: embeddings,
      genes: genes,
      clusters: clusters,
      metadata: {
        numCells: dataInfo.numCells,
        numGenes: dataInfo.numGenes,
        spatialDimensions: dataInfo.spatialDimensions,
        availableEmbeddings: dataInfo.availableEmbeddings,
        clusterCount: dataInfo.clusterCount,
        spatialScalingFactor: 1,
        wasNormalized,
      },
      matrix: matrix,
      allClusterColumnNames: clusterColumnInfo.names,
      allClusterColumnTypes: clusterColumnInfo.types,
      allEmbeddingNames: dataInfo.availableEmbeddings || [],
      normalized: false, // always false now — coords are always raw
      availableDeStatsColumns: adapter.getAvailableDeStatsColumns(),
    };

    console.log("[Worker] S3 loading complete");

    return serializable;
  },

  /**
   * Load a single embedding on demand (off main thread).
   * Creates a temporary adapter, initializes it, loads the embedding.
   */
  async loadEmbeddingFromS3(
    datasetId: string,
    embeddingName: string,
    customS3BaseUrl?: string,
  ): Promise<{ name: string; data: number[][] } | null> {
    const adapter = customS3BaseUrl
      ? new ChunkedDataAdapter("custom", undefined, customS3BaseUrl)
      : new ChunkedDataAdapter(datasetId);

    await adapter.initialize();

    return adapter.loadEmbedding(embeddingName);
  },

  /**
   * Load cluster column(s) on demand (off main thread).
   * Creates a temporary adapter, initializes it, loads the cluster data.
   */
  async loadClusterFromS3(
    datasetId: string,
    columnNames: string[],
    customS3BaseUrl?: string,
  ): Promise<Array<{
    column: string;
    type: string;
    values: any[];
    palette: Record<string, string> | null;
    uniqueValues?: string[];
  }> | null> {
    const adapter = customS3BaseUrl
      ? new ChunkedDataAdapter("custom", undefined, customS3BaseUrl)
      : new ChunkedDataAdapter(datasetId);

    await adapter.initialize();

    return adapter.loadClusters(columnNames);
  },
};

// Expose worker API via Comlink (only inside a real worker — the module is also
// imported by the parity harness in Node, which has no `self`).
if (typeof self !== "undefined") {
  Comlink.expose(workerApi);
}

export type StandardizedDatasetWorkerApi = typeof workerApi;
