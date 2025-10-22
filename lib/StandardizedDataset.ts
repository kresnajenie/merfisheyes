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

    // Find gene index
    const geneIndex = this.genes.indexOf(gene);
    if (geneIndex === -1) {
      return null;
    }

    // If matrix is already cached (from worker), use it directly
    if (this.matrix) {
      return this.extractColumnFromMatrix(this.matrix, geneIndex);
    }

    // Otherwise, need adapter to fetch matrix
    if (!this.adapter) {
      throw new Error("No adapter available for gene expression data access");
    }

    // Cache matrix for subsequent queries
    this.matrix = this.adapter.fetchFullMatrix();

    return this.adapter.fetchColumn(this.matrix, geneIndex);
  }

  /**
   * Extract a column from a matrix (works with different matrix formats)
   */
  private extractColumnFromMatrix(matrix: any, column: number): number[] {
    // Case 1: Map<string, Float32Array> (Xenium/MERSCOPE format)
    if (matrix instanceof Map) {
      const gene = this.genes[column];
      if (!gene || !matrix.has(gene)) {
        return [];
      }
      return Array.from(matrix.get(gene)!);
    }

    // Case 2: Array of arrays (row-major)
    if (Array.isArray(matrix) && Array.isArray(matrix[0])) {
      return matrix.map((row: any) => row[column]);
    }

    // Case 3: TypedArray (flattened row-major - H5AD format)
    if (ArrayBuffer.isView(matrix)) {
      const typedArray = matrix as any;
      // Assume it's flattened row-major: [row0col0, row0col1, ..., row1col0, row1col1, ...]
      const numCells = this.spatial.coordinates.length;
      const numGenes = this.genes.length;
      
      if (column >= numGenes) {
        throw new Error("Column index out of bounds");
      }

      return Array.from(
        { length: numCells },
        (_, i) => typedArray[i * numGenes + column],
      );
    }

    throw new Error("Unsupported matrix format");
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
    matrix?: any;
  }): StandardizedDataset {
    const dataset = new StandardizedDataset({
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

    // Pre-cache the matrix if provided (from worker)
    if (data.matrix) {
      dataset.matrix = data.matrix;
    }

    return dataset;
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
