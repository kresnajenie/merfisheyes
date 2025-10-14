import { H5adAdapter } from "./adapters/H5adAdapter";
import { XeniumAdapter } from "./adapters/XeniumAdapter";
import { MerscopeAdapter } from "./adapters/MerscopeAdapter";
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
  static async fromH5ad(file: File, onProgress?: (progress: number, message: string) => Promise<void> | void): Promise<StandardizedDataset> {
    const adapter = new H5adAdapter();
    await adapter.initialize(file, onProgress);

    // Load all data through adapter
    await onProgress?.(92, "Loading spatial coordinates...");
    const spatial = adapter.loadSpatialCoordinates();
    console.log("Spatial data:", spatial);
    await onProgress?.(94, "Loading embeddings...");
    const embeddings = adapter.loadEmbeddings();
    console.log("Embeddings:", embeddings);
    await onProgress?.(96, "Loading genes...");
    const genes = await adapter.loadGenes();
    console.log("Genes:", genes.length, "genes loaded");
    await onProgress?.(98, "Loading clusters...");
    const clusters = await adapter.loadClusters();
    console.log("Clusters:", clusters);
    const dataInfo = adapter.getDatasetInfo();
    console.log("Dataset info:", dataInfo);

    // Generate dataset ID
    const timestamp = Date.now();
    const random = Math.random().toString(36).substring(2, 11);
    const id = `h5ad_${file.name.replace(".h5ad", "")}_${timestamp}_${random}`;

    return new StandardizedDataset({
      id: id,
      name: file.name.replace(".h5ad", ""),
      type: "h5ad",
      spatial: spatial,
      embeddings: embeddings,
      genes: genes,
      clusters: clusters,
      metadata: {
        ...adapter.metadata,
        originalFileName: file.name,
        numCells: dataInfo.numCells,
        numGenes: dataInfo.numGenes,
      },
      rawData: adapter.h5File,
      adapter: adapter,
    });
  }

  /**
   * Create StandardizedDataset from Xenium files
   */
  static async fromXenium(files: File[], onProgress?: (progress: number, message: string) => Promise<void> | void): Promise<StandardizedDataset> {
    const adapter = new XeniumAdapter();
    await adapter.initialize(files, onProgress);

    // Load all data through adapter
    await onProgress?.(92, "Loading spatial coordinates...");
    const spatial = adapter.loadSpatialCoordinates();
    console.log("Spatial data:", spatial);
    await onProgress?.(94, "Loading embeddings...");
    const embeddings = adapter.loadEmbeddings();
    console.log("Embeddings:", embeddings);
    await onProgress?.(96, "Loading genes...");
    const genes = await adapter.loadGenes();
    console.log("Genes:", genes.length, "genes loaded");
    await onProgress?.(98, "Loading clusters...");
    const clusterData = await adapter.loadClusters();
    console.log("Clusters:", clusterData);

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
    console.log("Dataset info:", dataInfo);

    // Generate dataset ID and name from folder
    const timestamp = Date.now();
    const random = Math.random().toString(36).substring(2, 11);
    const folderName =
      files[0]?.webkitRelativePath?.split("/")[0] || "xenium_data";
    const id = `xenium_${folderName}_${timestamp}_${random}`;

    return new StandardizedDataset({
      id: id,
      name: folderName,
      type: "xenium",
      spatial: spatial,
      embeddings: embeddings,
      genes: genes,
      clusters: clusters,
      metadata: {
        ...((adapter as any).metadata || {}),
        fileCount: files.length,
        numCells: dataInfo.numCells,
        numGenes: dataInfo.numGenes,
      },
      rawData: files,
      adapter: adapter,
    });
  }

  /**
   * Create StandardizedDataset from MERSCOPE files
   */
  static async fromMerscope(files: File[], onProgress?: (progress: number, message: string) => Promise<void> | void): Promise<StandardizedDataset> {
    const adapter = new MerscopeAdapter();
    await adapter.initialize(files, onProgress);

    // Load all data through adapter
    await onProgress?.(92, "Loading spatial coordinates...");
    const spatial = adapter.loadSpatialCoordinates();
    console.log("Spatial data:", spatial);
    await onProgress?.(94, "Loading embeddings...");
    const rawEmbeddings = adapter.loadEmbeddings();
    console.log("Raw embeddings:", rawEmbeddings);

    // Convert embeddings to proper format (filter out undefined values)
    const embeddings: Record<string, number[][]> = {};
    if (rawEmbeddings && typeof rawEmbeddings === "object") {
      Object.entries(rawEmbeddings).forEach(([key, value]) => {
        if (value && Array.isArray(value)) {
          embeddings[key] = value;
        }
      });
    }

    await onProgress?.(96, "Loading genes...");
    const genes = await adapter.loadGenes();
    console.log("Genes:", genes.length, "genes loaded");
    await onProgress?.(98, "Loading clusters...");
    const clusterData = await adapter.loadClusters();
    console.log("Clusters:", clusterData);

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
    console.log("Dataset info:", dataInfo);

    // Generate dataset ID and name from folder
    const timestamp = Date.now();
    const random = Math.random().toString(36).substring(2, 11);
    const folderName =
      files[0]?.webkitRelativePath?.split("/")[0] || "merscope_data";
    const id = `merscope_${folderName}_${timestamp}_${random}`;

    return new StandardizedDataset({
      id: id,
      name: folderName,
      type: "merscope",
      spatial: spatial,
      embeddings: embeddings,
      genes: genes,
      clusters: clusters,
      metadata: {
        ...((adapter as any).metadata || {}),
        fileCount: files.length,
        numCells: dataInfo.numCells,
        numGenes: dataInfo.numGenes,
      },
      rawData: files,
      adapter: adapter,
    });
  }
}
