/**
 * Processed Single Molecule Adapter
 * Loads pre-processed single molecule data with lazy loading support
 *
 * Supports both:
 * - Remote: S3 storage via presigned URLs (dataset ID)
 * - Local: Local files from folder upload (File objects)
 *
 * Expected folder structure (from Python script):
 * output_folder/
 * ├── manifest.json.gz
 * └── genes/
 *     ├── GENE1.bin.gz
 *     ├── GENE2.bin.gz
 *     └── ...
 */

import pako from "pako";

interface Manifest {
  version: string;
  created_at: string;
  dataset_id: string;
  name: string;
  type: string;
  statistics: {
    total_molecules: number;
    unique_genes: number;
    spatial_dimensions: 2 | 3;
  };
  genes: {
    unique_gene_names: string[];
  };
  processing: {
    compression: string;
    coordinate_format: string;
    coordinate_range: string;
    scaling_factor: number;
    created_by: string;
    source_file: string;
  };
}

export class ProcessedSingleMoleculeAdapter {
  private datasetId: string;
  private downloadUrls: Record<string, string> = {};
  private localFiles: Map<string, File> | null = null;
  private manifest: Manifest | null = null;
  private mode: "remote" | "local";
  private geneCache: Map<string, number[]> = new Map(); // Cache loaded genes

  constructor(datasetId: string, localFiles?: Map<string, File>) {
    this.datasetId = datasetId;
    this.localFiles = localFiles || null;
    this.mode = localFiles ? "local" : "remote";
  }

  /**
   * Initialize adapter by loading manifest
   */
  async initialize() {
    console.log("[ProcessedSingleMoleculeAdapter] Initializing...", {
      datasetId: this.datasetId,
      mode: this.mode,
    });

    try {
      if (this.mode === "remote") {
        // Fetch dataset metadata and presigned URLs from API
        const baseUrl =
          typeof self !== "undefined" && self.location
            ? self.location.origin
            : "";
        const url = `${baseUrl}/api/single-molecule/${this.datasetId}`;
        const response = await fetch(url);

        if (!response.ok) {
          throw new Error(
            `Failed to fetch dataset: ${response.status} ${response.statusText}`,
          );
        }

        const data = await response.json();

        console.log(
          "[ProcessedSingleMoleculeAdapter] Dataset API response:",
          data,
        );

        // Check if there's an error in the response
        if (data.error) {
          throw new Error(`${data.error}: ${data.message || ""}`);
        }

        // For remote mode, we get the manifest URL
        if (data.manifestUrl) {
          this.downloadUrls["manifest.json.gz"] = data.manifestUrl;
        }
      } else {
        // Local mode - files already provided
        console.log(
          "[ProcessedSingleMoleculeAdapter] Local mode: using provided files",
          this.localFiles?.size,
          "files",
        );
      }

      // Load manifest
      this.manifest = await this.fetchManifest();
      console.log(
        "[ProcessedSingleMoleculeAdapter] Loaded manifest:",
        this.manifest,
      );
    } catch (error) {
      console.error(
        "[ProcessedSingleMoleculeAdapter] Initialization failed:",
        error,
      );
      throw error;
    }
  }

  /**
   * Fetch and parse manifest - supports both manifest.json.gz and manifest.json
   */
  private async fetchManifest(): Promise<Manifest> {
    // In local mode, check which manifest format is available
    if (this.mode === "local") {
      const hasGzipped = this.localFiles?.has("manifest.json.gz");
      const hasPlain = this.localFiles?.has("manifest.json");

      if (hasGzipped) {
        return (await this.fetchGzippedJSON("manifest.json.gz")) as Manifest;
      } else if (hasPlain) {
        return (await this.fetchPlainJSON("manifest.json")) as Manifest;
      }

      throw new Error(
        "Manifest not found in local files. Expected manifest.json.gz or manifest.json",
      );
    }

    // Remote mode always uses gzipped format
    return (await this.fetchGzippedJSON("manifest.json.gz")) as Manifest;
  }

  /**
   * Fetch and parse a plain (non-gzipped) JSON file
   */
  private async fetchPlainJSON(fileKey: string): Promise<any> {
    console.log(`[ProcessedSingleMoleculeAdapter] Fetching plain JSON ${fileKey}...`);

    if (this.mode === "local") {
      const file = this.localFiles?.get(fileKey);

      if (!file) {
        throw new Error(`File not found in local files: ${fileKey}`);
      }

      const text = await file.text();

      return JSON.parse(text);
    }

    // Remote mode
    const url = this.downloadUrls[fileKey];

    if (!url) {
      throw new Error(`No download URL for file: ${fileKey}`);
    }

    const response = await fetch(url);

    if (!response.ok) {
      throw new Error(`Failed to fetch ${fileKey}: ${response.statusText}`);
    }

    return response.json();
  }

  /**
   * Fetch a gzipped JSON file and parse it
   */
  private async fetchGzippedJSON(fileKey: string): Promise<any> {
    console.log(`[ProcessedSingleMoleculeAdapter] Fetching ${fileKey}...`);

    let arrayBuffer: ArrayBuffer;

    if (this.mode === "local") {
      // Local mode: read from File object
      const file = this.localFiles?.get(fileKey);

      if (!file) {
        throw new Error(`File not found in local files: ${fileKey}`);
      }

      arrayBuffer = await file.arrayBuffer();
    } else {
      // Remote mode: fetch from URL
      const url = this.downloadUrls[fileKey];

      if (!url) {
        throw new Error(`No download URL for file: ${fileKey}`);
      }

      const response = await fetch(url);

      if (!response.ok) {
        throw new Error(`Failed to fetch ${fileKey}: ${response.statusText}`);
      }

      arrayBuffer = await response.arrayBuffer();
    }

    // Decompress gzip
    const decompressed = pako.ungzip(new Uint8Array(arrayBuffer), {
      to: "string",
    });

    return JSON.parse(decompressed);
  }

  /**
   * Fetch a gzipped binary file and decompress it
   */
  private async fetchGzippedBinary(fileKey: string): Promise<Uint8Array> {
    console.log(`[ProcessedSingleMoleculeAdapter] Fetching ${fileKey}...`);

    let arrayBuffer: ArrayBuffer;

    if (this.mode === "local") {
      // Local mode: read from File object
      const file = this.localFiles?.get(fileKey);

      if (!file) {
        throw new Error(`File not found in local files: ${fileKey}`);
      }

      arrayBuffer = await file.arrayBuffer();
    } else {
      // Remote mode: fetch from presigned URL via API
      // For remote mode, we need to fetch the presigned URL for the gene file
      const geneName = fileKey.replace("genes/", "").replace(".bin.gz", "");
      const baseUrl =
        typeof self !== "undefined" && self.location
          ? self.location.origin
          : "";
      const url = `${baseUrl}/api/single-molecule/${this.datasetId}/gene/${encodeURIComponent(geneName)}`;

      const response = await fetch(url);

      if (!response.ok) {
        throw new Error(`Failed to fetch gene URL: ${response.statusText}`);
      }

      const data = await response.json();
      const presignedUrl = data.url;

      // Now fetch the actual gene file from S3
      const geneResponse = await fetch(presignedUrl);

      if (!geneResponse.ok) {
        throw new Error(
          `Failed to fetch gene file: ${geneResponse.statusText}`,
        );
      }

      arrayBuffer = await geneResponse.arrayBuffer();
    }

    // Decompress gzip
    const decompressed = pako.ungzip(new Uint8Array(arrayBuffer));

    return decompressed;
  }

  /**
   * Get coordinates for a specific gene (lazy loaded)
   * Returns flat array: [x1,y1,z1, x2,y2,z2, ...]
   */
  async getCoordinatesByGene(geneName: string): Promise<number[]> {
    // Check cache first
    if (this.geneCache.has(geneName)) {
      console.log(
        `[ProcessedSingleMoleculeAdapter] Gene ${geneName} found in cache`,
      );

      return this.geneCache.get(geneName)!;
    }

    if (!this.manifest) {
      throw new Error("Manifest not loaded. Call initialize() first.");
    }

    // Check if gene exists
    if (!this.manifest.genes.unique_gene_names.includes(geneName)) {
      throw new Error(
        `Gene '${geneName}' not found. Available genes: ${this.manifest.statistics.unique_genes}`,
      );
    }

    console.log(
      `[ProcessedSingleMoleculeAdapter] Loading gene ${geneName} from storage...`,
    );

    // Sanitize gene name for filename (matches Python script)
    const sanitizedName = geneName
      .replace(/[^a-zA-Z0-9]/g, "_")
      .replace(/_+/g, "_")
      .replace(/^_|_$/g, "");

    const fileKey = `genes/${sanitizedName}.bin.gz`;

    // Fetch and decompress gene file
    const decompressed = await this.fetchGzippedBinary(fileKey);

    // Convert Uint8Array to Float32Array
    const float32Array = new Float32Array(
      decompressed.buffer,
      decompressed.byteOffset,
      decompressed.byteLength / 4,
    );

    // Convert to regular number array for compatibility
    const coords = Array.from(float32Array);

    // Cache for future use
    this.geneCache.set(geneName, coords);

    console.log(
      `[ProcessedSingleMoleculeAdapter] Loaded gene ${geneName}: ${coords.length / 3} molecules`,
    );

    return coords;
  }

  /**
   * Get manifest data
   */
  getManifest(): Manifest {
    if (!this.manifest) {
      throw new Error("Manifest not loaded. Call initialize() first.");
    }

    return this.manifest;
  }

  /**
   * Get unique genes list
   */
  getUniqueGenes(): string[] {
    if (!this.manifest) {
      throw new Error("Manifest not loaded. Call initialize() first.");
    }

    return this.manifest.genes.unique_gene_names;
  }

  /**
   * Get dimensions (2 or 3)
   */
  getDimensions(): 2 | 3 {
    if (!this.manifest) {
      throw new Error("Manifest not loaded. Call initialize() first.");
    }

    return this.manifest.statistics.spatial_dimensions;
  }

  /**
   * Get scaling factor
   */
  getScalingFactor(): number {
    if (!this.manifest) {
      throw new Error("Manifest not loaded. Call initialize() first.");
    }

    return this.manifest.processing.scaling_factor;
  }

  /**
   * Get total molecule count
   */
  getTotalMolecules(): number {
    if (!this.manifest) {
      throw new Error("Manifest not loaded. Call initialize() first.");
    }

    return this.manifest.statistics.total_molecules;
  }

  /**
   * Get gene count
   */
  getGeneCount(): number {
    if (!this.manifest) {
      throw new Error("Manifest not loaded. Call initialize() first.");
    }

    return this.manifest.statistics.unique_genes;
  }

  /**
   * Clear gene cache (for memory management)
   */
  clearCache() {
    console.log("[ProcessedSingleMoleculeAdapter] Clearing gene cache");
    this.geneCache.clear();
  }

  /**
   * Get cache stats
   */
  getCacheStats() {
    return {
      cachedGenes: this.geneCache.size,
      totalGenes: this.manifest?.statistics.unique_genes || 0,
    };
  }
}
