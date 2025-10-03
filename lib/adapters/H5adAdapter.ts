/**
 * H5AD Data Adapter
 * Handles all H5AD-specific data operations and h5wasm interactions
 */

// We'll load h5wasm dynamically from CDN
import h5wasm from "h5wasm";

// async function loadH5wasm() {
//   if (h5wasm) return h5wasm;
//   const module = await import(
//     "https://cdn.jsdelivr.net/npm/h5wasm@0.4.9/dist/esm/hdf5_hl.js"
//   );
//   h5wasm = module.default;
//   return h5wasm;
// }

export class H5adAdapter {
  h5File: any = null;
  metadata: any = null;

  /**
   * Initialize adapter with H5AD file
   */
  async initialize(file: File) {
    const { FS } = await h5wasm.ready;

    const arrayBuffer = await file.arrayBuffer();
    FS.writeFile(file.name, new Uint8Array(arrayBuffer));
    this.h5File = new h5wasm.File(file.name, "r");

    // Get available keys from different groups
    const rootKeys = this.h5File.keys();
    const obsKeys = this.h5File.get("obs").keys();
    const obsmKeys = this.h5File.get("obsm")?.keys() || [];
    const unsKeys = this.h5File.get("uns")?.keys() || [];
    const varKeys = this.h5File.get("var")?.keys() || [];

    this.metadata = {
      rootKeys,
      obsKeys,
      obsmKeys,
      unsKeys,
      varKeys,
    };
  }

  /**
   * Get basic dataset information
   */
  getDatasetInfo() {
    const xDataset = this.h5File.get("X");
    const [numCells, numGenes] = xDataset.shape || [0, 0];
    return { numCells, numGenes };
  }

  /**
   * Helper method to decode bytes
   */
  decodeBytes(value: any): string {
    const decoder = new TextDecoder();
    if (value instanceof Uint8Array)
      return decoder.decode(value).replace(/\0/g, "");
    if (Array.isArray(value))
      return decoder.decode(new Uint8Array(value)).replace(/\0/g, "");
    return String(value);
  }

  /**
   * Fetch observation data (cell metadata)
   */
  async fetchObs(group: string): Promise<any[]> {
    const obs = this.h5File.get("obs");
    if (group === "index") return Array.from(obs.attrs._index || []);
    const field = obs.get(group);
    const keys = field.keys?.();
    if (keys?.includes("categories") && keys.includes("codes")) {
      const categories = Array.from(
        field.get("categories").value,
        this.decodeBytes.bind(this)
      );
      const codes = Array.from(field.get("codes").value) as number[];
      return codes.map((i) => categories[i] ?? "Unknown");
    }
    return Array.from(field.value, this.decodeBytes.bind(this));
  }

  /**
   * Fetch variable data (gene metadata)
   */
  async fetchVar(group: string): Promise<any[]> {
    const v = this.h5File.get("var");
    if (group === "index") {
      const varIndex = v.get("_index");
      return Array.from(varIndex.value, this.decodeBytes.bind(this));
    }
    const d = v.get(group);
    return Array.from(d.value, this.decodeBytes.bind(this));
  }

  /**
   * Fetch unstructured data
   */
  fetchUns(group: string): any {
    try {
      return this.h5File.get("uns/" + group).value;
    } catch (error) {
      return null;
    }
  }

  /**
   * Fetch observation matrix data
   */
  fetchObsm(group: string): any {
    const m = this.h5File.get("obsm/" + group);
    return m.toArray ? m.toArray() : m.value;
  }

  /**
   * Load and normalize coordinates from dataset
   */
  loadDatasetCoordinates(dataset: any): any {
    if (!dataset) return null;

    let coordinates;
    const flat = dataset.value;
    if (!flat || !flat.length) return null;

    // Determine if we have 2D or 3D coordinates based on array length
    const dimensions = dataset.metadata.shape[1];
    const rows = Math.floor(flat.length / dimensions);

    coordinates = Array.from({ length: rows }, (_, i) => {
      const coords = [flat[i * dimensions], flat[i * dimensions + 1]];
      if (dimensions === 3) {
        coords.push(flat[i * dimensions + 2]);
      }
      return coords;
    });

    return {
      coordinates: coordinates,
      dimensions: dimensions,
    };
  }

  /**
   * Load specific coordinates by key
   */
  loadCoordinates(coordKey: string = "X_spatial"): number[][] | null {
    try {
      const dataset = this.h5File.get("obsm/" + coordKey);
      const result = this.loadDatasetCoordinates(dataset);
      return result ? result.coordinates : null;
    } catch (error) {
      console.warn(`Failed to load coordinates for ${coordKey}:`, error);
      return null;
    }
  }

  /**
   * Load spatial coordinates
   */
  loadSpatialCoordinates(): { coordinates: number[][]; dimensions: number } {
    const result = this.loadCoordinates("X_spatial");
    if (!result) {
      throw new Error("No spatial coordinates found in obsm/X_spatial");
    }

    // Get dimensions from the dataset
    const dataset = this.h5File.get("obsm/X_spatial");
    const dimensions = dataset.metadata.shape[1];

    return {
      coordinates: result,
      dimensions: dimensions,
    };
  }

  /**
   * Load all available embeddings
   */
  loadEmbeddings(): Record<string, number[][]> {
    const embeddings: Record<string, number[][]> = {};

    for (const key of this.metadata.obsmKeys) {
      if (key.startsWith("X_") && key !== "X_spatial") {
        try {
          const coords = this.loadCoordinates(key);
          if (coords) {
            const embeddingName = key.replace("X_", "");
            embeddings[embeddingName] = coords;
          }
        } catch (error) {
          console.warn(`Failed to load embedding ${key}:`, error);
        }
      }
    }

    return embeddings;
  }

  /**
   * Load gene list
   */
  async loadGenes(): Promise<string[]> {
    try {
      const varIndex = await this.fetchVar("_index");
      return varIndex || [];
    } catch (error) {
      console.warn("Failed to load genes:", error);
      return [];
    }
  }

  /**
   * Determine if data is categorical or numerical
   */
  isCategoricalData(values: any[]): boolean {
    if (!values || values.length === 0) return false;

    // Check if values are strings
    if (typeof values[0] === "string") return true;

    // Check if numerical values have limited unique values (threshold: 50)
    const uniqueValues = new Set(values);
    if (uniqueValues.size <= 50) return true;

    // If more than 50 unique numerical values, treat as continuous
    return false;
  }

  /**
   * Parse a single obs column to determine type and palette
   */
  async parseObsColumn(columnName: string): Promise<any> {
    try {
      const values = await this.fetchObs(columnName);

      // Determine if categorical or numerical
      const isCategorical = this.isCategoricalData(values);

      if (isCategorical) {
        const palette = await this.loadClusterPalette(columnName);
        return {
          column: columnName,
          type: "categorical",
          values: values,
          palette: palette,
        };
      } else {
        return {
          column: columnName,
          type: "numerical",
          values: values,
          palette: null,
        };
      }
    } catch (error) {
      console.warn(`Failed to parse column ${columnName}:`, error);
      return null;
    }
  }

  /**
   * Load all cluster/obs columns
   */
  async loadClusters(): Promise<any[]> {
    const allClusters: any[] = [];

    // Get all non-index columns
    const nonIndexColumns = this.metadata.obsKeys.filter(
      (key: string) => !key.startsWith("_")
    );

    // Parse all columns
    for (const columnName of nonIndexColumns) {
      const clusterData = await this.parseObsColumn(columnName);
      if (clusterData) {
        allClusters.push(clusterData);
      }
    }

    return allClusters.length > 0 ? allClusters : [];
  }

  /**
   * Normalize hex color to 7 characters (# + 6 hex)
   */
  normalizeHexColor(hex: string): string {
    if (typeof hex !== "string") return "#808080";

    // If not the expected length, cut off last 2 characters if longer than 7
    if (hex.length !== 7 && hex.length > 7) {
      return hex.substring(0, 7);
    }

    return hex;
  }

  /**
   * Load color palette for clusters
   */
  async loadClusterPalette(
    clusterColumn: string
  ): Promise<Record<string, string>> {
    try {
      // Try to get colors from uns
      const colors = this.fetchUns(clusterColumn + "_colors");
      const clusters = await this.fetchObs(clusterColumn);
      const uniqueClusters = Array.from(new Set(clusters)).sort();

      if (colors && colors.length === uniqueClusters.length) {
        const palette: Record<string, string> = {};
        uniqueClusters.forEach((cluster, index) => {
          palette[cluster] = this.normalizeHexColor(colors[index]);
        });
        return palette;
      }
    } catch (error) {
      console.log("No colors found in uns, using default colors");
    }

    // Use default matplotlib colors
    const defaultColors = [
      "#1f77b4",
      "#ff7f0e",
      "#2ca02c",
      "#d62728",
      "#9467bd",
      "#8c564b",
      "#e377c2",
      "#7f7f7f",
      "#bcbd22",
      "#17becf",
    ];

    const clusters = await this.fetchObs(clusterColumn);
    const uniqueClusters = Array.from(new Set(clusters)).sort();
    const palette: Record<string, string> = {};

    uniqueClusters.forEach((cluster, index) => {
      palette[cluster] = this.normalizeHexColor(
        defaultColors[index % defaultColors.length]
      );
    });

    return palette;
  }

  /**
   * Get available observation columns
   */
  getObsColumns(): string[] {
    return this.metadata.obsKeys;
  }

  /**
   * Get available embeddings
   */
  getObsmEmbeddings(): string[] {
    return this.metadata.obsmKeys;
  }

  /**
   * Fetch gene expression data
   */
  async fetchGeneExpression(gene: string): Promise<number[] | null> {
    try {
      const matrix = this.fetchX();
      const genes = await this.fetchVar("_index");
      const index = genes.indexOf(gene);
      if (index === -1) throw new Error("Gene not found");
      return this.fetchColumn(matrix, index);
    } catch (error) {
      console.warn(`Failed to fetch gene expression for ${gene}:`, error);
      return null;
    }
  }

  /**
   * Fetch expression matrix
   */
  fetchX(): any {
    const m = this.h5File.get("X");
    return m.value;
  }

  /**
   * Fetch column from matrix
   */
  fetchColumn(matrix: any, column: number): number[] {
    if (Array.isArray(matrix[0])) return matrix.map((row: any) => row[column]);

    const xDataset = this.h5File.get("X");
    const [n_obs, n_genes] = xDataset.shape || [0, 0];

    if (ArrayBuffer.isView(matrix)) {
      if (column >= n_genes) throw new Error("Column index out of bounds");
      const typedArray = matrix as any; // Type assertion for indexing
      return Array.from(
        { length: n_obs },
        (_, i) => typedArray[i * n_genes + column]
      );
    }

    throw new Error("Unsupported matrix format");
  }
}
