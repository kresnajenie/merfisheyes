import { normalizeCoordinates } from "./utils/coordinates";
import { selectBestClusterColumnByName } from "./utils/dataset-utils";

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
  allClusterColumnNames: string[];
  allClusterColumnTypes: Record<string, string>;
  clustersFullyLoaded: boolean;
  allEmbeddingNames: string[];
  embeddingsFullyLoaded: boolean;

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
    this.allClusterColumnNames = [];
    this.allClusterColumnTypes = {};
    this.clustersFullyLoaded = true; // Default true; S3/chunked paths set false
    this.allEmbeddingNames = [];
    this.embeddingsFullyLoaded = true; // Default true; S3/chunked paths set false

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
   * Adds newly loaded cluster columns, deduplicating by column name.
   * Used by the background cluster loader to append remaining columns.
   */
  addClusters(
    newClusters: Array<{
      column: string;
      type: string;
      values: any[];
      palette: Record<string, string> | null;
    }>,
  ) {
    if (!newClusters || newClusters.length === 0) return;

    const existing = this.clusters || [];
    const existingNames = new Set(existing.map((c) => c.column));
    const toAdd = newClusters.filter((c) => !existingNames.has(c.column));

    if (toAdd.length > 0) {
      this.clusters = [...existing, ...toAdd];
    }
  }

  /**
   * Add a single embedding that was loaded on demand.
   */
  addEmbedding(name: string, data: number[][]) {
    this.embeddings[name] = data;
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

    // If adapter has fetchGeneExpression method (ChunkedDataAdapter), use it
    if (
      this.adapter &&
      typeof this.adapter.fetchGeneExpression === "function"
    ) {
      return await this.adapter.fetchGeneExpression(gene);
    }

    // Otherwise, need adapter to fetch full matrix
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
    allClusterColumnNames?: string[];
    allClusterColumnTypes?: Record<string, string>;
    allEmbeddingNames?: string[];
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

    // Set deferred cluster loading info if provided
    if (data.allClusterColumnNames) {
      dataset.allClusterColumnNames = data.allClusterColumnNames;
      dataset.allClusterColumnTypes = data.allClusterColumnTypes || {};
      // If not all columns are loaded yet, mark as not fully loaded
      const loadedColumns = new Set(
        (data.clusters || []).map((c) => c.column),
      );
      dataset.clustersFullyLoaded = data.allClusterColumnNames.every((name) =>
        loadedColumns.has(name),
      );
    }

    // Set deferred embedding loading info if provided
    if (data.allEmbeddingNames && data.allEmbeddingNames.length > 0) {
      dataset.allEmbeddingNames = data.allEmbeddingNames;
      const loadedEmbeddings = new Set(Object.keys(data.embeddings || {}));
      dataset.embeddingsFullyLoaded = data.allEmbeddingNames.every((name) =>
        loadedEmbeddings.has(name),
      );
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
    priorityColumn?: string,
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
      priorityColumn,
    );

    // Reconstruct StandardizedDataset from serialized data
    const dataset = StandardizedDataset.fromSerializedData(serializedData);

    // Override manifest's dataset_id with the actual database ID
    // (manifest may contain a client-generated ID that differs from the db-assigned ds_... ID)
    dataset.id = datasetId;

    // For S3 datasets, create a fresh adapter in the main thread
    // This allows on-demand gene expression loading and background cluster loading
    const { ChunkedDataAdapter } = await import(
      "./adapters/ChunkedDataAdapter"
    );
    const adapter = new ChunkedDataAdapter(datasetId);

    await adapter.initialize();

    // Attach adapter for on-demand gene expression queries
    dataset.adapter = adapter;

    // Read column info from adapter (already cached from initialize())
    const columnInfo = adapter.getClusterColumnInfo();
    dataset.allClusterColumnNames = columnInfo.names;
    dataset.allClusterColumnTypes = columnInfo.types;

    return dataset;
  }

  /**
   * Load dataset from custom S3 URL (user-owned bucket)
   * @param customS3BaseUrl - Base S3 URL to dataset folder (e.g., https://bucket.s3.region.amazonaws.com/path/to/folder)
   * @param onProgress - Optional progress callback
   * @returns StandardizedDataset
   */
  static async fromCustomS3(
    customS3BaseUrl: string,
    onProgress?: (progress: number, message: string) => Promise<void> | void,
    priorityColumnHint?: string,
  ): Promise<StandardizedDataset> {
    await onProgress?.(10, "Initializing custom S3 adapter...");

    // Create adapter with custom S3 base URL
    const { ChunkedDataAdapter } = await import(
      "./adapters/ChunkedDataAdapter"
    );
    const adapter = new ChunkedDataAdapter(
      "custom", // Dummy dataset ID
      undefined, // No local files
      customS3BaseUrl, // Custom S3 base URL
    );

    await onProgress?.(30, "Loading manifest from custom S3...");

    await adapter.initialize();

    await onProgress?.(50, "Loading spatial coordinates...");
    const spatial = await adapter.loadSpatialCoordinates();

    // Skip eager embedding loading — embeddings are loaded on demand
    const embeddings: Record<string, number[][]> = {};

    await onProgress?.(70, "Loading genes...");
    const genes = await adapter.loadGenes();

    // Deferred cluster loading: use URL hint if valid, otherwise auto-detect
    const columnInfo = adapter.getClusterColumnInfo();
    let priorityColumn: string | null = null;

    if (
      priorityColumnHint &&
      columnInfo.names.includes(priorityColumnHint)
    ) {
      priorityColumn = priorityColumnHint;
    } else {
      priorityColumn = selectBestClusterColumnByName(
        columnInfo.names,
        columnInfo.types,
      );
    }

    await onProgress?.(80, "Loading priority cluster column...");
    const clusters = priorityColumn
      ? await adapter.loadClusters([priorityColumn])
      : null;

    const dataInfo = adapter.getDatasetInfo();

    await onProgress?.(90, "Loading expression matrix...");
    const matrix = await adapter.fetchFullMatrix();

    await onProgress?.(95, "Finalizing dataset...");

    // Create StandardizedDataset
    const dataset = new StandardizedDataset({
      id: dataInfo.id || "custom",
      name: dataInfo.name || "Custom S3 Dataset",
      type: dataInfo.type || "custom",
      spatial,
      embeddings,
      genes,
      clusters,
      metadata: {
        numCells: dataInfo.numCells,
        numGenes: dataInfo.numGenes,
        spatialDimensions: dataInfo.spatialDimensions,
        availableEmbeddings: dataInfo.availableEmbeddings,
        clusterCount: dataInfo.clusterCount,
        loadedFrom: "custom_s3",
        customS3BaseUrl,
      },
      rawData: null,
      adapter,
    });

    // Assign matrix (pre-loaded for custom S3)
    dataset.matrix = matrix;

    // Set deferred cluster column info
    dataset.allClusterColumnNames = columnInfo.names;
    dataset.allClusterColumnTypes = columnInfo.types;
    dataset.clustersFullyLoaded = columnInfo.names.length <= 1;

    // Set deferred embedding loading info
    dataset.allEmbeddingNames = dataInfo.availableEmbeddings || [];
    dataset.embeddingsFullyLoaded = false;

    await onProgress?.(100, "Dataset loaded successfully");

    return dataset;
  }

  /**
   * Load dataset from local chunked files (created by Python script)
   */
  static async fromLocalChunked(
    files: File[],
    onProgress?: (progress: number, message: string) => Promise<void> | void,
    priorityColumnHint?: string,
  ): Promise<StandardizedDataset> {
    console.log("[StandardizedDataset] Loading from local chunked files...");

    // Convert File[] to Map<fileKey, File>
    // File keys should match the structure: manifest.json, coords/spatial.bin.gz, etc.
    const fileMap = new Map<string, File>();

    for (const file of files) {
      // Extract relative path from webkitRelativePath
      const relativePath = file.webkitRelativePath;

      if (!relativePath) {
        console.warn("File missing webkitRelativePath:", file.name);
        continue;
      }

      // Remove the root folder name to get the file key
      // e.g., "my_dataset/manifest.json" -> "manifest.json"
      const parts = relativePath.split("/");
      const fileKey = parts.slice(1).join("/"); // Remove first part (root folder)

      fileMap.set(fileKey, file);
      console.log(`Mapped file: ${fileKey}`);
    }

    console.log(`Total files mapped: ${fileMap.size}`);

    // Generate a temporary dataset ID
    const datasetId = `local_${Date.now()}`;

    await onProgress?.(10, "Initializing local chunked adapter...");

    // Create ChunkedDataAdapter in local mode
    const { ChunkedDataAdapter } = await import(
      "./adapters/ChunkedDataAdapter"
    );
    const adapter = new ChunkedDataAdapter(datasetId, fileMap);

    await adapter.initialize();

    await onProgress?.(30, "Loading spatial coordinates...");
    const spatial = await adapter.loadSpatialCoordinates();

    // Skip eager embedding loading — embeddings are loaded on demand
    const embeddings: Record<string, number[][]> = {};

    await onProgress?.(70, "Loading genes...");
    const genes = await adapter.loadGenes();

    // Deferred cluster loading: use URL hint if valid, otherwise auto-detect
    const columnInfo = adapter.getClusterColumnInfo();
    let priorityColumn: string | null = null;

    if (
      priorityColumnHint &&
      columnInfo.names.includes(priorityColumnHint)
    ) {
      priorityColumn = priorityColumnHint;
    } else {
      priorityColumn = selectBestClusterColumnByName(
        columnInfo.names,
        columnInfo.types,
      );
    }

    await onProgress?.(85, "Loading priority cluster column...");
    const clusters = priorityColumn
      ? await adapter.loadClusters([priorityColumn])
      : null;

    const dataInfo = adapter.getDatasetInfo();

    await onProgress?.(95, "Finalizing dataset...");

    // Expression matrix is NOT pre-loaded - will be fetched on-demand via adapter
    const matrix = null;

    // Create StandardizedDataset
    const dataset = new StandardizedDataset({
      id: datasetId,
      name: dataInfo.name,
      type: dataInfo.type,
      spatial: {
        coordinates: spatial.coordinates,
        dimensions: spatial.dimensions,
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
      },
      adapter: adapter,
      rawData: null,
    });

    // Attach matrix (even though it's null for chunked datasets)
    dataset.matrix = matrix;

    // Set deferred cluster column info
    dataset.allClusterColumnNames = columnInfo.names;
    dataset.allClusterColumnTypes = columnInfo.types;
    dataset.clustersFullyLoaded = columnInfo.names.length <= 1;

    // Set deferred embedding loading info
    dataset.allEmbeddingNames = dataInfo.availableEmbeddings || [];
    dataset.embeddingsFullyLoaded = false;

    await onProgress?.(100, "Dataset loaded successfully!");

    return dataset;
  }
}
