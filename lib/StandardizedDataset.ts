import { normalizeCoordinates } from "./utils/coordinates";

interface SpatialData {
  coordinates: number[][];
  dimensions: number;
  scalingFactor: number;
}

interface ClusterData {
  column: string;
  type: string;
  values: any[];
  palette: Record<string, string> | null;
}

interface StandardizedDatasetParams {
  id: string;
  name: string;
  type: string;
  spatial: {
    coordinates: number[][];
    dimensions: number;
  };
  embeddings?: Record<string, number[][]>;
  genes?: string[];
  clusters?: ClusterData[] | null;
  metadata?: Record<string, any>;
  rawData?: any;
  adapter?: any;
}

/**
 * Standardized dataset format for all data types
 * This class ensures consistent data structure across different parsers
 */
export class StandardizedDataset {
  id: string;
  name: string;
  type: string;
  spatial: SpatialData;
  embeddings: Record<string, number[][]>;
  genes: string[];
  clusters: ClusterData[] | null;
  metadata: Record<string, any>;
  rawData: any;
  adapter: any;
  matrix: any | null;

  constructor({
    id,
    name,
    type,
    spatial,
    embeddings = {},
    genes = [],
    clusters = null,
    metadata = {},
    rawData = null,
    adapter = null,
  }: StandardizedDatasetParams) {
    this.id = id;
    this.name = name;
    this.type = type;

    // Normalize spatial coordinates to [-1, 1]
    const normalizedSpatial = normalizeCoordinates(spatial.coordinates);

    this.spatial = {
      coordinates: normalizedSpatial?.normalized || spatial.coordinates,
      dimensions: spatial.dimensions,
      scalingFactor: normalizedSpatial?.scalingFactor || 1,
    };

    this.embeddings = embeddings;
    this.genes = genes;
    this.clusters = clusters;
    this.metadata = {
      ...metadata,
      spatialScalingFactor: normalizedSpatial?.scalingFactor || 1,
    };
    this.rawData = rawData;
    this.adapter = adapter;
    this.matrix = null;

    this.validateStructure();
  }

  /**
   * Validate the dataset structure
   */
  validateStructure() {
    if (!this.id || typeof this.id !== "string") {
      throw new Error("Dataset must have a valid string ID");
    }

    if (!this.name || typeof this.name !== "string") {
      throw new Error("Dataset must have a valid string name");
    }

    if (!this.type || typeof this.type !== "string") {
      throw new Error("Dataset must have a valid string type");
    }

    if (
      !this.spatial ||
      !this.spatial.coordinates ||
      !Array.isArray(this.spatial.coordinates)
    ) {
      throw new Error("Dataset must have valid spatial coordinates");
    }

    if (![2, 3].includes(this.spatial.dimensions)) {
      throw new Error("Spatial dimensions must be 2 or 3");
    }

    if (this.spatial.coordinates.length > 0) {
      const firstCoord = this.spatial.coordinates[0];

      if (
        !Array.isArray(firstCoord) ||
        firstCoord.length < this.spatial.dimensions
      ) {
        throw new Error("Spatial coordinates format is invalid");
      }
    }
  }

  /**
   * Get the number of data points
   */
  getPointCount(): number {
    return this.spatial.coordinates.length;
  }

  /**
   * Get dataset summary
   */
  getSummary() {
    return {
      id: this.id,
      name: this.name,
      type: this.type,
      pointCount: this.getPointCount(),
      spatialDimensions: this.spatial.dimensions,
      availableEmbeddings: Object.keys(this.embeddings),
      geneCount: this.genes.length,
      clusterColumns: this.clusters ? this.clusters.map((c) => c.column) : [],
      numClusterColumns: this.clusters ? this.clusters.length : 0,
    };
  }

  /**
   * Get gene expression data for a specific gene
   * @param gene - Gene name
   * @returns Array of expression values or null if not found
   */
  async getGeneExpression(gene: string | null): Promise<number[] | null> {
    if (!gene) {
      return null;
    }
    if (!this.adapter) {
      throw new Error("No adapter available for gene expression data access");
    }

    // Cache matrix for subsequent queries
    if (!this.matrix) {
      this.matrix = this.adapter.fetchFullMatrix();
    }

    // Use cached genes array instead of re-fetching
    const geneIndex = this.genes.indexOf(gene);

    if (geneIndex === -1) {
      return null;
    }

    return this.adapter.fetchColumn(this.matrix, geneIndex);
  }

  /**
   * Create StandardizedDataset from H5AD file
   */
  /**
   * Reconstruct StandardizedDataset from serializable data (used after web worker processing)
   */
  static fromSerializedData(data: {
    id: string;
    name: string;
    type: string;
    spatial: {
      coordinates: number[][];
      dimensions: number;
      scalingFactor: number;
    };
    embeddings: Record<string, number[][]>;
    genes: string[];
    clusters:
      | {
          column: string;
          type: string;
          values: any[];
          palette: Record<string, string> | null;
        }[]
      | null;
    metadata: Record<string, any>;
  }): StandardizedDataset {
    return new StandardizedDataset({
      id: data.id,
      name: data.name,
      type: data.type,
      spatial: data.spatial,
      embeddings: data.embeddings,
      genes: data.genes,
      clusters: data.clusters,
      metadata: data.metadata,
      rawData: null,
      adapter: null,
    });
  }

  static async fromH5ad(
    file: File,
    onProgress?: (progress: number, message: string) => Promise<void> | void,
  ): Promise<StandardizedDataset> {
    // Use web worker for parsing
    const { getStandardizedDatasetWorker } = await import(
      "./workers/standardizedDatasetWorkerManager"
    );
    const worker = await getStandardizedDatasetWorker();

    // Import Comlink dynamically
    const Comlink = await import("comlink");

    // Parse in worker with proxied progress callback
    const serializedData = await worker.parseH5ad(
      file,
      onProgress ? Comlink.proxy(onProgress) : undefined,
    );

    // Reconstruct StandardizedDataset from serialized data
    return StandardizedDataset.fromSerializedData(serializedData);
  }

  /**
   * Create StandardizedDataset from Xenium files
   */
  static async fromXenium(
    files: File[],
    onProgress?: (progress: number, message: string) => Promise<void> | void,
  ): Promise<StandardizedDataset> {
    // Use web worker for parsing
    const { getStandardizedDatasetWorker } = await import(
      "./workers/standardizedDatasetWorkerManager"
    );
    const worker = await getStandardizedDatasetWorker();

    // Import Comlink dynamically
    const Comlink = await import("comlink");

    // Parse in worker with proxied progress callback
    const serializedData = await worker.parseXenium(
      files,
      onProgress ? Comlink.proxy(onProgress) : undefined,
    );

    // Reconstruct StandardizedDataset from serialized data
    return StandardizedDataset.fromSerializedData(serializedData);
  }

  /**
   * Create StandardizedDataset from MERSCOPE files
   */
  static async fromMerscope(
    files: File[],
    onProgress?: (progress: number, message: string) => Promise<void> | void,
  ): Promise<StandardizedDataset> {
    // Use web worker for parsing
    const { getStandardizedDatasetWorker } = await import(
      "./workers/standardizedDatasetWorkerManager"
    );
    const worker = await getStandardizedDatasetWorker();

    // Import Comlink dynamically
    const Comlink = await import("comlink");

    // Parse in worker with proxied progress callback
    const serializedData = await worker.parseMerscope(
      files,
      onProgress ? Comlink.proxy(onProgress) : undefined,
    );

    // Reconstruct StandardizedDataset from serialized data
    return StandardizedDataset.fromSerializedData(serializedData);
  }

  /**
   * Create StandardizedDataset from S3 chunked data
   */
  static async fromS3(
    datasetId: string,
    onProgress?: (progress: number, message: string) => Promise<void> | void,
  ): Promise<StandardizedDataset> {
    // Use web worker for parsing
    const { getStandardizedDatasetWorker } = await import(
      "./workers/standardizedDatasetWorkerManager"
    );
    const worker = await getStandardizedDatasetWorker();

    // Import Comlink dynamically
    const Comlink = await import("comlink");

    // Parse in worker with proxied progress callback
    const serializedData = await worker.parseS3(
      datasetId,
      onProgress ? Comlink.proxy(onProgress) : undefined,
    );

    // Reconstruct StandardizedDataset from serialized data
    return StandardizedDataset.fromSerializedData(serializedData);
  }
}
