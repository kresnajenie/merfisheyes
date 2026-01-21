import Papa from "papaparse";
import { ungzip } from "pako";

import { hyparquetService } from "./services/hyparquetService";
import { normalizeCoordinates } from "./utils/coordinates";
import {
  MOLECULE_COLUMN_MAPPINGS,
  MoleculeDatasetType,
} from "./config/moleculeColumnMappings";
import { shouldFilterGene } from "./utils/gene-filters";

/**
 * Generate a random bright color for dark background
 */
function generateBrightColor(): string {
  const hue = Math.random() * 360;
  const saturation = 70 + Math.random() * 30; // 70-100%
  const lightness = 50 + Math.random() * 20; // 50-70%

  return `hsl(${hue}, ${saturation}%, ${lightness}%)`;
}

/**
 * Gene visualization properties
 */
export interface GeneProperties {
  color: string;
  size: number;
}

/**
 * Format elapsed time in a human-readable format
 */
function formatElapsedTime(ms: number): string {
  if (ms < 1000) {
    return `${ms.toFixed(0)}ms`;
  } else if (ms < 60000) {
    return `${(ms / 1000).toFixed(2)}s`;
  } else {
    const minutes = Math.floor(ms / 60000);
    const seconds = ((ms % 60000) / 1000).toFixed(1);

    return `${minutes}m ${seconds}s`;
  }
}

/**
 * Standardized dataset format for single molecule data
 * Stores molecule coordinates with fast gene-based lookup
 */
export class SingleMoleculeDataset {
  id: string;
  name: string;
  type: string;

  // Core data
  uniqueGenes: string[];

  // Fast gene lookup with pre-computed normalized coordinates
  private geneIndex: Map<string, number[]>; // gene -> normalized [x1,y1,z1, x2,y2,z2, ...]

  // Gene visualization properties (color and size for each gene)
  geneColors: Record<string, GeneProperties>;

  dimensions: 2 | 3;
  scalingFactor: number;
  metadata: Record<string, any>;
  rawData: any;

  constructor({
    id,
    name,
    type,
    uniqueGenes,
    geneIndex,
    dimensions,
    scalingFactor,
    metadata = {},
    rawData = null,
  }: {
    id: string;
    name: string;
    type: string;
    uniqueGenes: string[];
    geneIndex: Map<string, number[]>;
    dimensions: 2 | 3;
    scalingFactor: number;
    metadata?: Record<string, any>;
    rawData?: any;
  }) {
    this.id = id;
    this.name = name;
    this.type = type;
    this.uniqueGenes = uniqueGenes;
    this.geneIndex = geneIndex;
    this.dimensions = dimensions;
    this.scalingFactor = scalingFactor;
    this.metadata = {
      ...metadata,
      spatialScalingFactor: scalingFactor,
    };
    this.rawData = rawData;

    this.validateStructure();

    // Initialize gene colors from localStorage or generate new ones
    this.geneColors = this.initializeGeneColors();
  }

  /**
   * Initialize gene colors from localStorage or generate new ones
   * Persists across page reloads for the same dataset ID
   */
  private initializeGeneColors(): Record<string, GeneProperties> {
    const storageKey = `sm_gene_colors_${this.id}`;

    // Try to load from localStorage
    if (typeof window !== "undefined") {
      try {
        const stored = localStorage.getItem(storageKey);

        if (stored) {
          const parsed = JSON.parse(stored);

          console.log(
            `[SingleMoleculeDataset] Loaded gene colors from localStorage for dataset: ${this.id}`,
          );

          return parsed;
        }
      } catch (error) {
        console.warn(
          `[SingleMoleculeDataset] Failed to load gene colors from localStorage:`,
          error,
        );
      }
    }

    // Generate new colors for all genes
    console.log(
      `[SingleMoleculeDataset] Generating new gene colors for ${this.uniqueGenes.length} genes`,
    );
    const geneColors: Record<string, GeneProperties> = {};

    for (const gene of this.uniqueGenes) {
      geneColors[gene] = {
        color: generateBrightColor(),
        size: 1.0, // Default local size multiplier
      };
    }

    // Save to localStorage
    if (typeof window !== "undefined") {
      try {
        localStorage.setItem(storageKey, JSON.stringify(geneColors));
        console.log(
          `[SingleMoleculeDataset] Saved gene colors to localStorage for dataset: ${this.id}`,
        );
      } catch (error) {
        console.warn(
          `[SingleMoleculeDataset] Failed to save gene colors to localStorage:`,
          error,
        );
      }
    }

    return geneColors;
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

    if (!Array.isArray(this.uniqueGenes)) {
      throw new Error("uniqueGenes must be an array");
    }

    if (!(this.geneIndex instanceof Map)) {
      throw new Error("geneIndex must be a Map");
    }

    if (![2, 3].includes(this.dimensions)) {
      throw new Error("Dimensions must be 2 or 3");
    }
  }

  /**
   * Get the number of molecules
   */
  getMoleculeCount(): number {
    let total = 0;

    for (const coords of this.geneIndex.values()) {
      total += coords.length / 3; // Each molecule has x, y, z
    }

    return total;
  }

  /**
   * Get dataset summary
   */
  getSummary() {
    return {
      id: this.id,
      name: this.name,
      type: this.type,
      moleculeCount: this.getMoleculeCount(),
      spatialDimensions: this.dimensions,
      uniqueGenes: this.uniqueGenes.length,
    };
  }

  /**
   * Get gene index entries for serialization (used by web worker)
   */
  getGeneIndexEntries(): [string, number[]][] {
    return Array.from(this.geneIndex.entries());
  }

  /**
   * Reconstruct dataset from serializable data (used after web worker processing)
   */
  static fromSerializedData(data: {
    id: string;
    name: string;
    type: string;
    uniqueGenes: string[];
    geneIndexEntries: [string, number[]][];
    dimensions: 2 | 3;
    scalingFactor: number;
    metadata: Record<string, any>;
  }): SingleMoleculeDataset {
    return new SingleMoleculeDataset({
      id: data.id,
      name: data.name,
      type: data.type,
      uniqueGenes: data.uniqueGenes,
      geneIndex: new Map(data.geneIndexEntries),
      dimensions: data.dimensions,
      scalingFactor: data.scalingFactor,
      metadata: data.metadata,
      rawData: null,
    });
  }

  /**
   * Get coordinates for a specific gene
   * Returns flat array: [x1,y1,z1, x2,y2,z2, ...]
   * Throws error if gene not found
   */
  getCoordinatesByGene(geneName: string): number[] {
    // Check if gene exists
    if (!this.geneIndex.has(geneName)) {
      throw new Error(
        `Gene '${geneName}' not found. Available genes: ${this.uniqueGenes.length}`,
      );
    }

    // Direct lookup - already pre-computed!
    return this.geneIndex.get(geneName)!;
  }

  /**
   * Get all genes (for backward compatibility)
   */
  get genes(): string[] {
    return this.uniqueGenes;
  }

  /**
   * Create SingleMoleculeDataset from Parquet file
   * Automatically normalizes coordinates to [-1, 1] range
   */
  static async fromParquet(
    file: File,
    datasetType: MoleculeDatasetType = "xenium",
    onProgress?: (progress: number, message: string) => Promise<void> | void,
  ): Promise<SingleMoleculeDataset> {
    const startTime = performance.now();

    console.log(
      `[SingleMoleculeDataset] Starting parquet parsing: ${file.name}`,
    );

    // Get column mapping for this dataset type
    const columnMapping = MOLECULE_COLUMN_MAPPINGS[datasetType];

    // Determine which columns to read
    const columnsToRead = [
      columnMapping.gene,
      columnMapping.x,
      columnMapping.y,
    ];

    // Add z column if it exists (for 3D data)
    if (columnMapping.z) {
      columnsToRead.push(columnMapping.z);
    }

    // Read parquet file using hyparquet
    const columnData = await hyparquetService.readParquetColumns(
      file,
      columnsToRead,
      onProgress,
    );

    await onProgress?.(30, "Extracting columns...");

    // Extract columns from the returned Map
    const moleculeGenes = columnData.get(columnMapping.gene);
    const xData = columnData.get(columnMapping.x);
    const yData = columnData.get(columnMapping.y);
    const zData = columnData.get(columnMapping.z || "");

    if (!moleculeGenes || !xData || !yData) {
      throw new Error(
        `Missing required columns. Expected: ${columnMapping.gene}, ${columnMapping.x}, ${columnMapping.y}`,
      );
    }

    await onProgress?.(50, "Converting to typed arrays...");

    // Convert to typed arrays for efficiency
    const xCoords = new Float32Array(xData);
    const yCoords = new Float32Array(yData);
    const zCoords = zData
      ? new Float32Array(zData)
      : new Float32Array(xCoords.length); // Fill with 0s if 2D

    const dimensions: 2 | 3 = zData ? 3 : 2;

    await onProgress?.(60, "Normalizing coordinates...");

    // Normalize coordinates to [-1, 1]
    const coords2D: number[][] = [];

    for (let i = 0; i < xCoords.length; i++) {
      coords2D.push([xCoords[i], yCoords[i], zCoords[i]]);
    }

    const normalized = normalizeCoordinates(coords2D);
    let scalingFactor = 1;
    let normalizedCoords = coords2D;

    if (normalized) {
      scalingFactor = normalized.scalingFactor;
      normalizedCoords = normalized.normalized;
    }

    await onProgress?.(70, "Building gene index...");

    // Build gene index with pre-computed normalized coordinates
    const geneIndex = new Map<string, number[]>();
    const uniqueGenesSet = new Set<string>();

    const totalMolecules = moleculeGenes.length;
    const progressInterval = Math.max(1, Math.floor(totalMolecules / 20)); // Report every 5%

    for (let i = 0; i < moleculeGenes.length; i++) {
      const gene = moleculeGenes[i];

      uniqueGenesSet.add(gene);

      if (!geneIndex.has(gene)) {
        geneIndex.set(gene, []);
      }

      // Store normalized coordinates [x, y, z] for this molecule
      const coords = geneIndex.get(gene)!;

      coords.push(
        normalizedCoords[i][0],
        normalizedCoords[i][1],
        normalizedCoords[i][2],
      );

      // Report progress every 5% and yield to browser
      if (i > 0 && i % progressInterval === 0) {
        const elapsed = performance.now() - startTime;
        const progress = 70 + Math.floor((i / totalMolecules) * 20); // 70-90%

        await onProgress?.(
          progress,
          `Indexing molecules: ${((i / totalMolecules) * 100).toFixed(1)}% (${formatElapsedTime(elapsed)})`,
        );
        // Yield to browser to allow UI updates
        await new Promise((resolve) => setTimeout(resolve, 0));
      }
    }

    await onProgress?.(85, "Filtering control genes...");

    // Filter out control/unassigned genes
    const filteredGenes: string[] = [];
    let filteredCount = 0;

    for (const gene of Array.from(uniqueGenesSet)) {
      if (shouldFilterGene(gene)) {
        geneIndex.delete(gene);
        filteredCount++;
      } else {
        filteredGenes.push(gene);
      }
    }

    if (filteredCount > 0) {
      console.log(
        `[SingleMoleculeDataset] Filtered ${filteredCount} control/unassigned genes from parquet dataset`,
      );
    }

    const uniqueGenes = filteredGenes;

    await onProgress?.(90, "Creating dataset...");

    // Generate dataset ID
    const timestamp = Date.now();
    const random = Math.random().toString(36).substring(2, 11);
    const id = `parquet_${file.name.replace(/\.(parquet|csv)$/, "")}_${timestamp}_${random}`;

    const dataset = new SingleMoleculeDataset({
      id,
      name: file.name.replace(/\.(parquet|csv)$/, ""),
      type: datasetType,
      uniqueGenes,
      geneIndex,
      dimensions,
      scalingFactor,
      metadata: {
        originalFileName: file.name,
        moleculeCount: xCoords.length,
        uniqueGeneCount: uniqueGenes.length,
        datasetType,
        columnMapping,
      },
      rawData: file,
    });

    const elapsedTime = performance.now() - startTime;

    console.log(
      `[SingleMoleculeDataset] ✅ Parquet parsing complete: ${formatElapsedTime(elapsedTime)} | ` +
        `${dataset.getMoleculeCount().toLocaleString()} molecules | ` +
        `${dataset.genes.length.toLocaleString()} genes`,
    );
    await onProgress?.(
      100,
      `Dataset loaded successfully in ${formatElapsedTime(elapsedTime)}`,
    );

    return dataset;
  }

  /**
   * Create SingleMoleculeDataset from CSV file
   * Automatically normalizes coordinates to [-1, 1] range
   */
  static async fromCSV(
    file: File,
    datasetType: MoleculeDatasetType = "xenium",
    onProgress?: (progress: number, message: string) => Promise<void> | void,
  ): Promise<SingleMoleculeDataset> {
    const startTime = performance.now();

    console.log(`[SingleMoleculeDataset] Starting CSV parsing: ${file.name}`);

    await onProgress?.(10, "Reading CSV file...");

    // Get column mapping for this dataset type
    const columnMapping = MOLECULE_COLUMN_MAPPINGS[datasetType];

    // Parse CSV file
    const text = await file.text();

    await onProgress?.(30, "Parsing CSV data...");

    const parseResult = Papa.parse(text, {
      header: true,
      skipEmptyLines: true,
      dynamicTyping: true,
    });

    if (parseResult.errors.length > 0) {
      console.warn("CSV parsing warnings:", parseResult.errors);
    }

    const rows = parseResult.data as any[];

    await onProgress?.(50, "Extracting columns...");

    // Extract columns from parsed data
    const moleculeGenes: string[] = [];
    const xCoords: number[] = [];
    const yCoords: number[] = [];
    const zCoords: number[] = [];

    let hasZ = false;

    for (let i = 0; i < rows.length; i++) {
      const row = rows[i];

      // Skip invalid rows
      if (!row || !row[columnMapping.gene]) continue;

      moleculeGenes.push(String(row[columnMapping.gene]));
      xCoords.push(Number(row[columnMapping.x]) || 0);
      yCoords.push(Number(row[columnMapping.y]) || 0);

      if (
        columnMapping.z &&
        row[columnMapping.z] !== undefined &&
        row[columnMapping.z] !== null
      ) {
        zCoords.push(Number(row[columnMapping.z]) || 0);
        hasZ = true;
      } else {
        zCoords.push(0);
      }

      // Report progress every 5% and yield to browser
      if (i > 0 && i % Math.max(1, Math.floor(rows.length / 20)) === 0) {
        const elapsed = performance.now() - startTime;
        const progress = 50 + Math.floor((i / rows.length) * 10); // 50-60%

        await onProgress?.(
          progress,
          `Extracting data: ${((i / rows.length) * 100).toFixed(1)}% (${formatElapsedTime(elapsed)})`,
        );
        // Yield to browser to allow UI updates
        await new Promise((resolve) => setTimeout(resolve, 0));
      }
    }

    const dimensions: 2 | 3 = hasZ ? 3 : 2;

    await onProgress?.(60, "Normalizing coordinates...");

    // Normalize coordinates to [-1, 1]
    const coords2D: number[][] = [];

    for (let i = 0; i < xCoords.length; i++) {
      coords2D.push([xCoords[i], yCoords[i], zCoords[i]]);
    }

    const normalized = normalizeCoordinates(coords2D);
    let scalingFactor = 1;
    let normalizedCoords = coords2D;

    if (normalized) {
      scalingFactor = normalized.scalingFactor;
      normalizedCoords = normalized.normalized;
    }

    await onProgress?.(70, "Building gene index...");

    // Build gene index with pre-computed normalized coordinates
    const geneIndex = new Map<string, number[]>();
    const uniqueGenesSet = new Set<string>();

    const totalMolecules = moleculeGenes.length;
    const progressInterval = Math.max(1, Math.floor(totalMolecules / 20)); // Report every 5%

    for (let i = 0; i < moleculeGenes.length; i++) {
      const gene = moleculeGenes[i];

      uniqueGenesSet.add(gene);

      if (!geneIndex.has(gene)) {
        geneIndex.set(gene, []);
      }

      // Store normalized coordinates [x, y, z] for this molecule
      const coords = geneIndex.get(gene)!;

      coords.push(
        normalizedCoords[i][0],
        normalizedCoords[i][1],
        normalizedCoords[i][2],
      );

      // Report progress every 5% and yield to browser
      if (i > 0 && i % progressInterval === 0) {
        const elapsed = performance.now() - startTime;
        const progress = 70 + Math.floor((i / totalMolecules) * 20); // 70-90%

        await onProgress?.(
          progress,
          `Indexing molecules: ${((i / totalMolecules) * 100).toFixed(1)}% (${formatElapsedTime(elapsed)})`,
        );
        // Yield to browser to allow UI updates
        await new Promise((resolve) => setTimeout(resolve, 0));
      }
    }

    await onProgress?.(85, "Filtering control genes...");

    // Filter out control/unassigned genes
    const filteredGenes: string[] = [];
    let filteredCount = 0;

    for (const gene of Array.from(uniqueGenesSet)) {
      if (shouldFilterGene(gene)) {
        geneIndex.delete(gene);
        filteredCount++;
      } else {
        filteredGenes.push(gene);
      }
    }

    if (filteredCount > 0) {
      console.log(
        `[SingleMoleculeDataset] Filtered ${filteredCount} control/unassigned genes from CSV dataset`,
      );
    }

    const uniqueGenes = filteredGenes;

    await onProgress?.(90, "Creating dataset...");

    // Generate dataset ID
    const timestamp = Date.now();
    const random = Math.random().toString(36).substring(2, 11);
    const id = `csv_${file.name.replace(/\.(parquet|csv)$/, "")}_${timestamp}_${random}`;

    const dataset = new SingleMoleculeDataset({
      id,
      name: file.name.replace(/\.(parquet|csv)$/, ""),
      type: datasetType,
      uniqueGenes,
      geneIndex,
      dimensions,
      scalingFactor,
      metadata: {
        originalFileName: file.name,
        moleculeCount: xCoords.length,
        uniqueGeneCount: uniqueGenes.length,
        datasetType,
        columnMapping,
      },
      rawData: file,
    });

    const elapsedTime = performance.now() - startTime;

    console.log(
      `[SingleMoleculeDataset] ✅ CSV parsing complete: ${formatElapsedTime(elapsedTime)} | ` +
        `${dataset.getMoleculeCount().toLocaleString()} molecules | ` +
        `${dataset.genes.length.toLocaleString()} genes`,
    );
    await onProgress?.(
      100,
      `Dataset loaded successfully in ${formatElapsedTime(elapsedTime)}`,
    );

    return dataset;
  }

  /**
   * Create SingleMoleculeDataset from S3 with lazy loading
   * Only loads manifest initially, gene data loaded on-demand and cached
   */
  static async fromS3(
    datasetId: string,
    onProgress?: (progress: number, message: string) => Promise<void> | void,
  ): Promise<SingleMoleculeDataset> {
    const startTime = performance.now();

    console.log(`[SingleMoleculeDataset] Loading from S3: ${datasetId}`);

    await onProgress?.(10, "Fetching dataset metadata...");

    // Fetch dataset metadata and manifest URL from API
    const response = await fetch(`/api/single-molecule/${datasetId}`);

    if (!response.ok) {
      const error = await response.json();

      throw new Error(error.error || "Failed to fetch dataset metadata");
    }

    const apiData = await response.json();

    await onProgress?.(30, "Downloading manifest...");

    // Download and decompress manifest from S3
    const manifestResponse = await fetch(apiData.manifestUrl);

    if (!manifestResponse.ok) {
      throw new Error("Failed to download manifest from S3");
    }

    const manifestCompressed = await manifestResponse.arrayBuffer();
    const manifestJson = ungzip(new Uint8Array(manifestCompressed), {
      to: "string",
    });
    const manifest = JSON.parse(manifestJson);

    await onProgress?.(60, "Initializing dataset...");

    // Extract metadata from manifest
    const uniqueGenes = manifest.genes.unique_gene_names;
    const dimensions = manifest.statistics.spatial_dimensions;
    const scalingFactor = manifest.processing.scaling_factor;

    // Create empty gene index - genes will be loaded on-demand
    const geneIndex = new Map<string, number[]>();

    // Create dataset with lazy-loading capability
    const dataset = new SingleMoleculeDataset({
      id: datasetId,
      name: manifest.name || apiData.title,
      type: manifest.type,
      uniqueGenes,
      geneIndex,
      dimensions,
      scalingFactor,
      metadata: {
        ...manifest,
        loadedFrom: "s3",
        manifestUrl: apiData.manifestUrl,
        moleculeCount: manifest.statistics.total_molecules,
        uniqueGeneCount: manifest.statistics.unique_genes,
      },
      rawData: null,
    });

    // Override getCoordinatesByGene to support lazy loading from S3
    const originalGetCoordinates = dataset.getCoordinatesByGene.bind(dataset);

    dataset.getCoordinatesByGene = async function (
      geneName: string,
    ): Promise<number[]> {
      // Check if already cached
      if (geneIndex.has(geneName)) {
        return geneIndex.get(geneName)!;
      }

      // Check if gene exists in manifest
      if (!uniqueGenes.includes(geneName)) {
        throw new Error(
          `Gene '${geneName}' not found. Available genes: ${uniqueGenes.length}`,
        );
      }

      console.log(
        `[SingleMoleculeDataset] Lazy-loading gene '${geneName}' from S3...`,
      );

      // Get presigned URL from API
      const urlResponse = await fetch(
        `/api/single-molecule/${datasetId}/gene/${encodeURIComponent(geneName)}`,
      );

      if (!urlResponse.ok) {
        const error = await urlResponse.json();

        throw new Error(
          error.error || `Failed to get URL for gene '${geneName}'`,
        );
      }

      const { url: geneUrl } = await urlResponse.json();

      // Download and decompress gene file
      const geneResponse = await fetch(geneUrl);

      if (!geneResponse.ok) {
        throw new Error(`Failed to download gene file for '${geneName}'`);
      }

      const geneCompressed = await geneResponse.arrayBuffer();
      const geneBuffer = ungzip(new Uint8Array(geneCompressed));

      // Convert to Float32Array, then to regular array
      const float32Array = new Float32Array(geneBuffer.buffer);
      const coordinates = Array.from(float32Array);

      // Cache for future use
      geneIndex.set(geneName, coordinates);

      console.log(
        `[SingleMoleculeDataset] ✅ Loaded gene '${geneName}': ${coordinates.length / dimensions} molecules`,
      );

      return coordinates;
    } as any; // Type override for lazy loading

    const elapsedTime = performance.now() - startTime;

    console.log(
      `[SingleMoleculeDataset] ✅ S3 dataset initialized: ${formatElapsedTime(elapsedTime)} | ` +
        `${manifest.statistics.total_molecules.toLocaleString()} molecules | ` +
        `${uniqueGenes.length.toLocaleString()} genes (lazy-loaded)`,
    );
    await onProgress?.(
      100,
      `Dataset ready in ${formatElapsedTime(elapsedTime)}`,
    );

    return dataset;
  }

  /**
   * Create SingleMoleculeDataset from locally uploaded chunked folder
   * Uses ProcessedSingleMoleculeAdapter for lazy loading from local files
   */
  static async fromLocalChunked(
    files: File[],
    onProgress?: (progress: number, message: string) => Promise<void> | void,
  ): Promise<SingleMoleculeDataset> {
    const startTime = performance.now();

    console.log(
      "[SingleMoleculeDataset] Loading from local chunked files...",
      files.length,
      "files",
    );

    await onProgress?.(10, "Preparing file map...");

    // Convert File[] to Map<fileKey, File>
    // File keys should match the structure: manifest.json.gz, genes/GENE1.bin.gz, etc.
    const fileMap = new Map<string, File>();

    for (const file of files) {
      // Extract relative path from webkitRelativePath
      const relativePath = file.webkitRelativePath;

      if (!relativePath) {
        console.warn("File missing webkitRelativePath:", file.name);
        continue;
      }

      // Remove the root folder name to get the file key
      // e.g., "my_dataset/manifest.json.gz" -> "manifest.json.gz"
      const parts = relativePath.split("/");
      const fileKey = parts.slice(1).join("/"); // Remove first part (root folder)

      fileMap.set(fileKey, file);
      console.log(`  Mapped file: ${fileKey}`);
    }

    console.log(`  Total files mapped: ${fileMap.size}`);

    // Generate a temporary dataset ID
    const datasetId = `local_sm_${Date.now()}`;

    await onProgress?.(30, "Initializing adapter...");

    // Create ProcessedSingleMoleculeAdapter in local mode
    const { ProcessedSingleMoleculeAdapter } = await import(
      "./adapters/ProcessedSingleMoleculeAdapter"
    );
    const adapter = new ProcessedSingleMoleculeAdapter(datasetId, fileMap);

    await adapter.initialize();

    await onProgress?.(60, "Loading manifest...");

    // Get manifest from adapter
    const manifest = adapter.getManifest();

    await onProgress?.(80, "Creating dataset...");

    // Create empty gene index - genes will be loaded on-demand
    const geneIndex = new Map<string, number[]>();

    // Create dataset with lazy-loading capability
    const dataset = new SingleMoleculeDataset({
      id: datasetId,
      name: `local_chunked_${Date.now()}`,
      type: "processed_chunked",
      uniqueGenes: manifest.uniqueGenes,
      geneIndex,
      dimensions: manifest.dimensions,
      scalingFactor: manifest.scalingFactor,
      metadata: {
        loadedFrom: "local_chunked",
        totalMolecules: manifest.totalMolecules,
        geneCount: manifest.geneCount,
        isPreChunked: true, // Mark as pre-chunked
      },
      rawData: null,
    });

    // Override getCoordinatesByGene to support lazy loading from local files
    dataset.getCoordinatesByGene = async function (
      geneName: string,
    ): Promise<number[]> {
      // Delegate to adapter which handles caching
      return adapter.getCoordinatesByGene(geneName);
    };

    // Attach adapter to dataset for future use
    (dataset as any).adapter = adapter;

    // Attach file map for upload
    (dataset as any).fileMap = fileMap;

    const elapsedTime = performance.now() - startTime;

    console.log(
      `[SingleMoleculeDataset] ✅ Local chunked dataset ready: ${formatElapsedTime(elapsedTime)} | ` +
        `${manifest.totalMolecules.toLocaleString()} molecules | ` +
        `${manifest.geneCount.toLocaleString()} genes (lazy-loaded)`,
    );

    await onProgress?.(
      100,
      `Dataset ready in ${formatElapsedTime(elapsedTime)}`,
    );

    return dataset;
  }
}
