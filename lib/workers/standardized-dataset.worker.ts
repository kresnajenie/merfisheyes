// lib/workers/standardized-dataset.worker.ts
import * as Comlink from "comlink";

import { H5adAdapter } from "../adapters/H5adAdapter";
import { XeniumAdapter } from "../adapters/XeniumAdapter";
import { MerscopeAdapter } from "../adapters/MerscopeAdapter";
import { ChunkedDataAdapter } from "../adapters/ChunkedDataAdapter";
import { normalizeCoordinates } from "../utils/coordinates";

/**
 * Progress callback type that can be proxied by Comlink
 */
type ProgressCallback = (
  progress: number,
  message: string,
) => void | Promise<void>;

/**
 * Serializable cluster data
 */
interface SerializableClusterData {
  column: string;
  type: string;
  values: any[];
  palette: Record<string, string> | null;
}

/**
 * Serializable spatial data
 */
interface SerializableSpatialData {
  coordinates: number[][];
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
}

/**
 * Worker API exposed via Comlink
 */
const workerApi = {
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

    console.log("[Worker] Embeddings:", embeddings);
    await onProgress?.(96, "Loading genes...");
    const genes = await adapter.loadGenes();

    console.log("[Worker] Genes:", genes.length, "genes loaded");
    await onProgress?.(98, "Loading clusters...");
    const clusters = await adapter.loadClusters();

    console.log("[Worker] Clusters:", clusters);
    const dataInfo = adapter.getDatasetInfo();

    console.log("[Worker] Dataset info:", dataInfo);

    // Normalize spatial coordinates
    const normalizedSpatial = normalizeCoordinates(spatial.coordinates);

    // Generate dataset ID
    const timestamp = Date.now();
    const random = Math.random().toString(36).substring(2, 11);
    const id = `h5ad_${file.name.replace(".h5ad", "")}_${timestamp}_${random}`;

    const serializable: SerializableStandardizedDataset = {
      id: id,
      name: file.name.replace(".h5ad", ""),
      type: "h5ad",
      spatial: {
        coordinates: normalizedSpatial?.normalized || spatial.coordinates,
        dimensions: spatial.dimensions,
        scalingFactor: normalizedSpatial?.scalingFactor || 1,
      },
      embeddings: embeddings,
      genes: genes,
      clusters: clusters,
      metadata: {
        ...adapter.metadata,
        originalFileName: file.name,
        numCells: dataInfo.numCells,
        numGenes: dataInfo.numGenes,
        spatialScalingFactor: normalizedSpatial?.scalingFactor || 1,
      },
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
    const clusterData = await adapter.loadClusters();

    console.log("[Worker] Clusters:", clusterData);

    // Wrap single cluster object in array and add type field for StandardizedDataset format
    const clusters = clusterData
      ? [
          {
            column: clusterData.column,
            type: "categorical",
            values: clusterData.values,
            palette: clusterData.palette,
          },
        ]
      : null;

    const dataInfo = adapter.getDatasetInfo();

    console.log("[Worker] Dataset info:", dataInfo);

    // Normalize spatial coordinates
    const normalizedSpatial = normalizeCoordinates(spatial.coordinates);

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
        coordinates: normalizedSpatial?.normalized || spatial.coordinates,
        dimensions: spatial.dimensions,
        scalingFactor: normalizedSpatial?.scalingFactor || 1,
      },
      embeddings: embeddings,
      genes: genes,
      clusters: clusters,
      metadata: {
        ...((adapter as any).metadata || {}),
        fileCount: files.length,
        numCells: dataInfo.numCells,
        numGenes: dataInfo.numGenes,
        spatialScalingFactor: normalizedSpatial?.scalingFactor || 1,
      },
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
    const clusterData = await adapter.loadClusters();

    console.log("[Worker] Clusters:", clusterData);

    // Wrap single cluster object in array and add type field for StandardizedDataset format
    const clusters = clusterData
      ? [
          {
            column: clusterData.column,
            type: "categorical",
            values: clusterData.values,
            palette: clusterData.palette,
          },
        ]
      : null;

    const dataInfo = adapter.getDatasetInfo();

    console.log("[Worker] Dataset info:", dataInfo);

    // Normalize spatial coordinates
    const normalizedSpatial = normalizeCoordinates(spatial.coordinates);

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
        coordinates: normalizedSpatial?.normalized || spatial.coordinates,
        dimensions: spatial.dimensions,
        scalingFactor: normalizedSpatial?.scalingFactor || 1,
      },
      embeddings: embeddings,
      genes: genes,
      clusters: clusters,
      metadata: {
        ...((adapter as any).metadata || {}),
        fileCount: files.length,
        numCells: dataInfo.numCells,
        numGenes: dataInfo.numGenes,
        spatialScalingFactor: normalizedSpatial?.scalingFactor || 1,
      },
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
  ): Promise<SerializableStandardizedDataset> {
    console.log("[Worker] Starting S3 loading:", datasetId);

    const adapter = new ChunkedDataAdapter(datasetId);

    await adapter.initialize();

    // Load all data through adapter (all methods return Promises)
    await onProgress?.(30, "Loading spatial coordinates...");
    const spatial = await adapter.loadSpatialCoordinates();

    console.log("[Worker] Spatial data:", spatial);
    await onProgress?.(50, "Loading embeddings...");
    const embeddings = await adapter.loadEmbeddings();

    console.log("[Worker] Embeddings:", embeddings);
    await onProgress?.(70, "Loading genes...");
    const genes = await adapter.loadGenes();

    console.log("[Worker] Genes:", genes.length, "genes loaded");
    await onProgress?.(90, "Loading clusters...");
    const clusters = await adapter.loadClusters();

    console.log("[Worker] Clusters:", clusters);
    const dataInfo = adapter.getDatasetInfo();

    console.log("[Worker] Dataset info:", dataInfo);

    // Normalize spatial coordinates
    const normalizedSpatial = normalizeCoordinates(spatial.coordinates);

    const serializable: SerializableStandardizedDataset = {
      id: dataInfo.id,
      name: dataInfo.name,
      type: dataInfo.type,
      spatial: {
        coordinates: normalizedSpatial?.normalized || spatial.coordinates,
        dimensions: spatial.dimensions,
        scalingFactor: normalizedSpatial?.scalingFactor || 1,
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
        spatialScalingFactor: normalizedSpatial?.scalingFactor || 1,
      },
    };

    console.log("[Worker] S3 loading complete");

    return serializable;
  },
};

// Expose worker API via Comlink
Comlink.expose(workerApi);

export type StandardizedDatasetWorkerApi = typeof workerApi;
