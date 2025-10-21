// lib/workers/single-molecule.worker.ts
import * as Comlink from "comlink";
import { SingleMoleculeDataset } from "../SingleMoleculeDataset";
import type { MoleculeDatasetType } from "../config/moleculeColumnMappings";

/**
 * Progress callback type that can be proxied by Comlink
 */
type ProgressCallback = (progress: number, message: string) => void;

/**
 * Serializable dataset data that can be transferred to main thread
 */
interface SerializableDatasetData {
  id: string;
  name: string;
  type: string;
  uniqueGenes: string[];
  geneIndexEntries: [string, number[]][]; // Map entries as array for serialization
  dimensions: 2 | 3;
  scalingFactor: number;
  metadata: Record<string, any>;
}

/**
 * Worker API exposed via Comlink
 */
const workerApi = {
  /**
   * Parse parquet file and return serializable dataset data
   */
  async parseParquet(
    file: File,
    datasetType: MoleculeDatasetType,
    onProgress?: ProgressCallback
  ): Promise<SerializableDatasetData> {
    console.log("[Worker] Starting parquet parsing:", file.name);

    // Parse using SingleMoleculeDataset
    // onProgress is already a Comlink proxy, just pass it directly
    const dataset = await SingleMoleculeDataset.fromParquet(
      file,
      datasetType,
      onProgress
    );

    console.log("[Worker] Parsing complete, serializing data...");

    // Convert to serializable format
    const serializable: SerializableDatasetData = {
      id: dataset.id,
      name: dataset.name,
      type: dataset.type,
      uniqueGenes: dataset.uniqueGenes,
      geneIndexEntries: Array.from(dataset.getGeneIndexEntries()), // Convert Map to array
      dimensions: dataset.dimensions,
      scalingFactor: dataset.scalingFactor,
      metadata: dataset.metadata,
    };

    console.log("[Worker] Serialization complete");
    return serializable;
  },

  /**
   * Parse CSV file and return serializable dataset data
   */
  async parseCSV(
    file: File,
    datasetType: MoleculeDatasetType,
    onProgress?: ProgressCallback
  ): Promise<SerializableDatasetData> {
    console.log("[Worker] Starting CSV parsing:", file.name);

    // onProgress is already a Comlink proxy, just pass it directly
    const dataset = await SingleMoleculeDataset.fromCSV(
      file,
      datasetType,
      onProgress
    );

    console.log("[Worker] Parsing complete, serializing data...");

    const serializable: SerializableDatasetData = {
      id: dataset.id,
      name: dataset.name,
      type: dataset.type,
      uniqueGenes: dataset.uniqueGenes,
      geneIndexEntries: Array.from(dataset.getGeneIndexEntries()),
      dimensions: dataset.dimensions,
      scalingFactor: dataset.scalingFactor,
      metadata: dataset.metadata,
    };

    console.log("[Worker] Serialization complete");
    return serializable;
  },
};

// Expose worker API via Comlink
Comlink.expose(workerApi);

export type SingleMoleculeWorkerApi = typeof workerApi;
