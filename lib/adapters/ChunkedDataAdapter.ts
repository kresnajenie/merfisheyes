/**
 * Chunked Data Adapter
 * Reconstructs StandardizedDataset from chunked compressed files
 * Supports three modes:
 * - Remote: S3 storage via presigned URLs (dataset ID)
 * - Local: Local files from folder upload (File objects)
 * - Custom: Direct S3 URLs (custom S3 base URL)
 */
export class ChunkedDataAdapter {
  private datasetId: string;
  private downloadUrls: Record<string, string> = {};
  private localFiles: Map<string, File> | null = null;
  private customS3BaseUrl: string | null = null;
  private manifest: any = null;
  private expressionIndex: any = null;
  private loadedChunks = new Map<number, any>();
  private obsMetadata: any = null;
  private mode: "remote" | "local" | "custom";

  constructor(
    datasetId: string,
    localFiles?: Map<string, File>,
    customS3BaseUrl?: string,
  ) {
    this.datasetId = datasetId;
    this.localFiles = localFiles || null;
    this.customS3BaseUrl = customS3BaseUrl || null;

    // Determine mode based on parameters
    if (customS3BaseUrl) {
      this.mode = "custom";
    } else if (localFiles) {
      this.mode = "local";
    } else {
      this.mode = "remote";
    }
  }

  /**
   * Initialize adapter by fetching presigned URLs (remote), using local files (local), or custom S3 URLs (custom)
   */
  async initialize() {
    console.log("Initializing ChunkedDataAdapter...", {
      datasetId: this.datasetId,
      mode: this.mode,
      customS3BaseUrl: this.customS3BaseUrl,
    });

    try {
      if (this.mode === "remote") {
        // Fetch dataset metadata and presigned URLs from API
        // Use absolute URL for worker compatibility
        // 'self' is available in both workers and main thread
        const baseUrl =
          typeof self !== "undefined" && self.location
            ? self.location.origin
            : "";
        const url = `${baseUrl}/api/datasets/${this.datasetId}`;
        const response = await fetch(url);

        if (!response.ok) {
          throw new Error(
            `Failed to fetch dataset: ${response.status} ${response.statusText}`,
          );
        }

        const data = await response.json();

        console.log("Dataset API response:", data);

        // Check if there's an error in the response (e.g., dataset not ready)
        if (data.error) {
          throw new Error(`${data.error}: ${data.message || ""}`);
        }

        // Check if files object exists
        if (!data.files || typeof data.files !== "object") {
          throw new Error(
            `Invalid response structure: files object missing. Status: ${data.status || "unknown"}`,
          );
        }

        this.downloadUrls = data.files;
        console.log(
          "Available files:",
          Object.keys(this.downloadUrls).length,
          "files",
        );
      } else if (this.mode === "custom") {
        // Custom S3 mode - files will be fetched directly using base URL
        console.log("Custom S3 mode: using base URL", this.customS3BaseUrl);
      } else {
        // Local mode - files already provided
        console.log(
          "Local mode: using provided files",
          this.localFiles?.size,
          "files",
        );
      }

      // Load manifest
      // For custom S3 mode, try both .json and .json.gz formats
      if (this.mode === "custom") {
        try {
          this.manifest = await this.fetchJSON("manifest.json");
          console.log("Loaded manifest (uncompressed):", this.manifest);
        } catch (error) {
          console.log(
            "Failed to load manifest.json, trying manifest.json.gz...",
          );
          // If manifest.json doesn't exist, it might be compressed
          // Try to fetch and decompress it
          try {
            const buffer = await this.fetchBinary("manifest.json.gz");
            const jsonString = new TextDecoder().decode(buffer);

            this.manifest = JSON.parse(jsonString);
            console.log("Loaded manifest (compressed):", this.manifest);
          } catch (gzError) {
            console.error(
              "Failed to load both manifest.json and manifest.json.gz",
            );
            throw error; // Throw original error
          }
        }
      } else {
        this.manifest = await this.fetchJSON("manifest.json");
        console.log("Loaded manifest:", this.manifest);
      }

      // Load expression index
      this.expressionIndex = await this.fetchJSON("expr/index.json");
      console.log("Loaded expression index:", {
        totalGenes: this.expressionIndex.total_genes,
        numChunks: this.expressionIndex.num_chunks,
        chunkSize: this.expressionIndex.chunk_size,
      });

      // Load and cache observation metadata
      this.obsMetadata = await this.fetchJSON("obs/metadata.json");
      console.log(
        "Loaded observation metadata:",
        Object.keys(this.obsMetadata),
      );

      return true;
    } catch (error) {
      console.error("Failed to initialize ChunkedDataAdapter:", error);
      throw new Error(`Adapter initialization failed: ${error}`);
    }
  }

  /**
   * Fetch and parse JSON file (from S3 URL or local File)
   */
  private async fetchJSON(fileKey: string) {
    console.log(`Fetching JSON: ${fileKey} (mode: ${this.mode})`);

    if (this.mode === "local") {
      const file = this.localFiles?.get(fileKey);

      if (!file) {
        throw new Error(`No local file found for ${fileKey}`);
      }

      const text = await file.text();

      return JSON.parse(text);
    } else if (this.mode === "custom") {
      // Custom S3 mode - construct URL directly
      const url = `${this.customS3BaseUrl}/${fileKey}`;

      console.log(`Fetching JSON from custom S3: ${url}`);

      const response = await fetch(url);

      if (!response.ok) {
        throw new Error(
          `Failed to fetch ${fileKey} from custom S3: ${response.status} ${response.statusText}`,
        );
      }

      return await response.json();
    } else {
      // Remote mode - use presigned URLs from API
      const url = this.downloadUrls[fileKey];

      if (!url) {
        throw new Error(`No download URL found for ${fileKey}`);
      }

      const response = await fetch(url);

      if (!response.ok) {
        throw new Error(
          `Failed to fetch ${fileKey}: ${response.status} ${response.statusText}`,
        );
      }

      return await response.json();
    }
  }

  /**
   * Fetch and decompress binary file (from S3 URL or local File)
   */
  private async fetchBinary(fileKey: string): Promise<ArrayBuffer> {
    console.log(`Fetching binary: ${fileKey} (mode: ${this.mode})`);

    if (this.mode === "local") {
      const file = this.localFiles?.get(fileKey);

      if (!file) {
        throw new Error(`No local file found for ${fileKey}`);
      }

      // File is already gzipped, decompress it
      return await this.decompress(file);
    } else if (this.mode === "custom") {
      // Custom S3 mode - construct URL directly
      const url = `${this.customS3BaseUrl}/${fileKey}`;

      console.log(`Fetching from custom S3: ${url}`);

      const response = await fetch(url);

      if (!response.ok) {
        throw new Error(
          `Failed to fetch ${fileKey} from custom S3: ${response.status} ${response.statusText}`,
        );
      }

      const compressedBlob = await response.blob();

      return await this.decompress(compressedBlob);
    } else {
      // Remote mode - use presigned URLs from API
      const url = this.downloadUrls[fileKey];

      if (!url) {
        throw new Error(`No download URL found for ${fileKey}`);
      }

      const response = await fetch(url);

      if (!response.ok) {
        throw new Error(
          `Failed to fetch ${fileKey}: ${response.status} ${response.statusText}`,
        );
      }

      const compressedBlob = await response.blob();

      return await this.decompress(compressedBlob);
    }
  }

  /**
   * Decompress gzip data
   */
  private async decompress(compressedBlob: Blob): Promise<ArrayBuffer> {
    const stream = compressedBlob
      .stream()
      .pipeThrough(new DecompressionStream("gzip"));

    return await new Response(stream).arrayBuffer();
  }

  /**
   * Load spatial coordinates
   */
  async loadSpatialCoordinates() {
    console.log("Loading spatial coordinates...");

    try {
      const spatialBuffer = await this.fetchBinary("coords/spatial.bin.gz");
      const coordinates = this.parseCoordinateBuffer(spatialBuffer);

      console.log("Loaded spatial coordinates:", coordinates.data.length);

      return {
        coordinates: coordinates.data,
        dimensions: coordinates.dimensions,
      };
    } catch (error) {
      console.error("Failed to load spatial coordinates:", error);
      throw error;
    }
  }

  /**
   * Load embeddings (UMAP, etc.)
   */
  async loadEmbeddings() {
    console.log("Loading embeddings...");
    const embeddings: Record<string, number[][]> = {};

    // Use manifest.statistics.available_embeddings as the source of truth.
    const availableEmbeddings = Array.isArray(
      this.manifest?.statistics?.available_embeddings,
    )
      ? this.manifest.statistics.available_embeddings
      : [];

    const embeddingKeys = Array.from(new Set(availableEmbeddings))
      .map((value) =>
        String(value)
          .trim()
          .replace(/\.bin\.gz$/i, ""),
      )
      .filter((value) => value.length > 0);

    for (const coordType of embeddingKeys) {
      try {
        const buffer = await this.fetchBinary(`coords/${coordType}.bin.gz`);
        const coordinates = this.parseCoordinateBuffer(buffer);

        embeddings[coordType] = coordinates.data;
        console.log(
          `Loaded ${coordType} embedding:`,
          coordinates.data.length,
          "points",
        );
      } catch (error) {
        console.warn(`Failed to load ${coordType} embedding:`, error);
      }
    }

    return embeddings;
  }

  /**
   * Load a single embedding by name (on-demand).
   * Reuses fetchBinary + parseCoordinateBuffer.
   */
  async loadEmbedding(
    name: string,
  ): Promise<{ name: string; data: number[][] } | null> {
    const coordType = name
      .trim()
      .replace(/\.bin\.gz$/i, "");

    if (!coordType) return null;

    try {
      const buffer = await this.fetchBinary(`coords/${coordType}.bin.gz`);
      const coordinates = this.parseCoordinateBuffer(buffer);

      console.log(
        `Loaded ${coordType} embedding on demand:`,
        coordinates.data.length,
        "points",
      );

      return { name: coordType, data: coordinates.data };
    } catch (error) {
      console.warn(`Failed to load ${coordType} embedding on demand:`, error);

      return null;
    }
  }

  /**
   * Parse binary coordinate buffer
   */
  private parseCoordinateBuffer(buffer: ArrayBuffer) {
    const view = new DataView(buffer);

    // Read header
    const numPoints = view.getUint32(0, true);
    const dimensions = view.getUint32(4, true);

    // Read coordinates
    const coordinates: number[][] = [];
    let offset = 8;

    for (let i = 0; i < numPoints; i++) {
      const coord: number[] = [];

      for (let d = 0; d < dimensions; d++) {
        coord.push(view.getFloat32(offset, true));
        offset += 4;
      }
      coordinates.push(coord);
    }

    return { data: coordinates, dimensions };
  }

  /**
   * Load gene names from expression index
   */
  async loadGenes(): Promise<string[]> {
    if (!this.expressionIndex) {
      throw new Error("Expression index not loaded");
    }

    return this.expressionIndex.genes.map((gene: any) => gene.name);
  }

  /**
   * Load clusters and color palettes
   */
  /**
   * Returns cluster column names and types from cached obsMetadata.
   * Zero network cost — obsMetadata is loaded during initialize().
   */
  getClusterColumnInfo(): { names: string[]; types: Record<string, string> } {
    if (!this.obsMetadata) {
      return { names: [], types: {} };
    }

    const names = Object.keys(this.obsMetadata);
    const types: Record<string, string> = {};
    for (const name of names) {
      types[name] = this.obsMetadata[name].type || "categorical";
    }

    return { names, types };
  }

  async loadClusters(columns?: string[]): Promise<Array<{
    column: string;
    type: string;
    values: any[];
    palette: Record<string, string> | null;
  }> | null> {
    console.log("Loading clusters...", columns ? `(filtered: ${columns.join(", ")})` : "(all)");

    try {
      // Use cached observation metadata to load all cluster columns
      if (!this.obsMetadata) {
        throw new Error("Observation metadata not loaded");
      }
      console.log(
        "Available observation columns:",
        Object.keys(this.obsMetadata),
      );

      const availableColumns = columns || Object.keys(this.obsMetadata);

      if (availableColumns.length === 0) {
        console.warn("No observation columns found");

        return null;
      }

      // Load all cluster columns
      const clusters = [];

      for (const columnName of availableColumns) {
        try {
          console.log(`Loading cluster column: ${columnName}`);

          // Load cluster values
          const clusterValues = await this.fetchCompressedJSON(
            `obs/${columnName}.json.gz`,
          );

          // Get column type from metadata
          const columnType = this.obsMetadata[columnName].type || "categorical";
          const isNumerical = columnType === "numerical";

          // Only load/generate palette for categorical columns
          let palette: Record<string, string> | null = null;

          if (!isNumerical) {
            try {
              palette = await this.fetchJSON(`palettes/${columnName}.json`);
              console.log(`Loaded palette from: palettes/${columnName}.json`);
            } catch (error) {
              // Fall back to default colors if palette not found
              console.log(
                `Palette palettes/${columnName}.json not found, generating default colors`,
              );
              palette = this.generateDefaultPalette(clusterValues);
            }
          } else {
            console.log(`Skipping palette for numerical column: ${columnName}`);
          }

          clusters.push({
            column: columnName,
            type: columnType,
            values: clusterValues,
            palette: palette,
          });
        } catch (error) {
          console.warn(`Failed to load cluster column ${columnName}:`, error);
          // Continue loading other columns even if one fails
        }
      }

      console.log(`Successfully loaded ${clusters.length} cluster columns`);

      return clusters.length > 0 ? clusters : null;
    } catch (error) {
      console.error("Failed to load clusters:", error);

      return null;
    }
  }

  /**
   * Fetch and decompress JSON data
   */
  private async fetchCompressedJSON(fileKey: string) {
    const buffer = await this.fetchBinary(fileKey);
    const jsonString = new TextDecoder().decode(buffer);

    return JSON.parse(jsonString);
  }

  /**
   * Generate default color palette for clusters
   */
  private generateDefaultPalette(values: any[]): Record<string, string> {
    const uniqueValues = [...new Set(values)];
    const colors = [
      "#e6194b",
      "#3cb44b",
      "#ffe119",
      "#4363d8",
      "#f58231",
      "#911eb4",
      "#42d4f4",
      "#f032e6",
      "#bfef45",
      "#fabed4",
      "#469990",
      "#dcbeff",
      "#9a6324",
      "#fffac8",
      "#800000",
      "#aaffc3",
      "#808000",
      "#ffd8b1",
      "#000075",
      "#a9a9a9",
    ];

    const palette: Record<string, string> = {};

    uniqueValues.forEach((value, index) => {
      palette[String(value)] = colors[index % colors.length];
    });

    return palette;
  }

  /**
   * Fetch gene expression data for a specific gene
   */
  async fetchGeneExpression(geneName: string): Promise<number[] | null> {
    console.log(`Fetching gene expression for: ${geneName}`);

    if (!this.expressionIndex) {
      throw new Error("Expression index not loaded");
    }

    // Find gene in index
    const geneInfo = this.expressionIndex.genes.find(
      (g: any) => g.name === geneName,
    );

    if (!geneInfo) {
      console.warn(`Gene not found: ${geneName}`);

      return null;
    }

    const chunkId = geneInfo.chunk_id;
    const positionInChunk = geneInfo.position_in_chunk;

    // Load chunk if not already cached
    let chunk = this.loadedChunks.get(chunkId);

    if (!chunk) {
      chunk = await this.loadExpressionChunk(chunkId);
      this.loadedChunks.set(chunkId, chunk);
    }

    // Extract gene data from chunk
    const geneData = chunk.genes[positionInChunk];

    if (!geneData) {
      throw new Error(
        `Gene data not found in chunk ${chunkId} at position ${positionInChunk}`,
      );
    }

    // Reconstruct dense array from sparse data
    const numCells = this.manifest.statistics.total_cells;
    const denseArray = new Float32Array(numCells);

    // Fill in non-zero values
    for (let i = 0; i < geneData.indices.length; i++) {
      const cellIndex = geneData.indices[i];
      const value = geneData.values[i];

      denseArray[cellIndex] = value;
    }

    console.log(
      `Loaded gene ${geneName}: ${geneData.indices.length} non-zero values out of ${numCells} cells`,
    );

    return Array.from(denseArray);
  }

  /**
   * Load and parse expression chunk
   */
  private async loadExpressionChunk(chunkId: number) {
    const chunkFilename = `chunk_${chunkId.toString().padStart(5, "0")}.bin.gz`;

    console.log(`Loading expression chunk: ${chunkFilename}`);

    const buffer = await this.fetchBinary(`expr/${chunkFilename}`);

    return this.parseExpressionChunk(buffer);
  }

  /**
   * Parse binary expression chunk
   */
  private parseExpressionChunk(buffer: ArrayBuffer) {
    const view = new DataView(buffer);

    // Read header
    const version = view.getUint32(0, true);
    const numGenes = view.getUint32(4, true);
    const chunkId = view.getUint32(8, true);
    const totalCells = view.getUint32(12, true);

    console.log(
      `Parsing chunk ${chunkId}: ${numGenes} genes, ${totalCells} total cells`,
    );

    // Read gene table
    const genes: any[] = [];
    let offset = 16; // After header

    for (let i = 0; i < numGenes; i++) {
      const geneTableOffset = offset + i * 24;
      const geneIndex = view.getUint32(geneTableOffset, true);
      const dataOffset = view.getUint32(geneTableOffset + 4, true);
      const dataSize = view.getUint32(geneTableOffset + 8, true);
      const uncompressedSize = view.getUint32(geneTableOffset + 12, true);
      const numNonZero = view.getUint32(geneTableOffset + 16, true);

      // Parse sparse gene data
      const geneData = this.parseSparseGeneData(buffer, dataOffset, numNonZero);

      genes.push({
        index: geneIndex,
        indices: geneData.indices,
        values: geneData.values,
        numNonZero: numNonZero,
      });
    }

    return {
      chunkId,
      numGenes,
      genes,
    };
  }

  /**
   * Parse sparse gene data from buffer
   */
  private parseSparseGeneData(
    buffer: ArrayBuffer,
    offset: number,
    numNonZero: number,
  ) {
    const view = new DataView(buffer);

    // Read sparse data header
    const numCells = view.getUint32(offset, true);
    const actualNonZero = view.getUint32(offset + 4, true);

    let dataOffset = offset + 8;

    // Read indices
    const indices: number[] = [];

    for (let i = 0; i < actualNonZero; i++) {
      indices.push(view.getUint32(dataOffset, true));
      dataOffset += 4;
    }

    // Read values
    const values: number[] = [];

    for (let i = 0; i < actualNonZero; i++) {
      values.push(view.getFloat32(dataOffset, true));
      dataOffset += 4;
    }

    return { indices, values };
  }

  /**
   * Get dataset info from manifest
   */
  getDatasetInfo() {
    if (!this.manifest) {
      throw new Error("Manifest not loaded");
    }

    return {
      id: this.manifest.dataset_id,
      name: this.manifest.name,
      type: this.manifest.type,
      numCells: this.manifest.statistics.total_cells,
      numGenes: this.manifest.statistics.total_genes,
      spatialDimensions: this.manifest.statistics.spatial_dimensions,
      availableEmbeddings: this.manifest.statistics.available_embeddings,
      clusterCount: this.manifest.statistics.cluster_count,
    };
  }

  /**
   * Fetch full expression matrix (placeholder for interface compatibility)
   * For chunked data, we don't load the full matrix at once
   * Returns null - actual data is fetched per-gene via fetchColumn
   */
  fetchFullMatrix(): null {
    // Don't actually load full matrix for S3 chunked data
    // This is just for interface compatibility with H5adAdapter
    return null;
  }

  /**
   * Fetch column (gene expression) from matrix by gene index
   * Uses chunk caching - if the gene's chunk is already loaded, it's instant
   * @param matrix - Ignored for chunked adapter (always null)
   * @param geneIndex - Index of the gene in the genes array
   * @returns Promise with array of expression values for all cells
   */
  async fetchColumn(matrix: any, geneIndex: number): Promise<number[]> {
    if (!this.expressionIndex) {
      throw new Error("Expression index not loaded");
    }

    // Get gene info by index
    const geneInfo = this.expressionIndex.genes[geneIndex];

    if (!geneInfo) {
      throw new Error(`Gene at index ${geneIndex} not found`);
    }

    const geneName = geneInfo.name;

    console.log(`Fetching column for gene index ${geneIndex}: ${geneName}`);

    // Use the existing fetchGeneExpression which handles chunk caching
    const result = await this.fetchGeneExpression(geneName);

    if (result === null) {
      throw new Error(`Failed to fetch expression data for gene: ${geneName}`);
    }

    return result;
  }
}
