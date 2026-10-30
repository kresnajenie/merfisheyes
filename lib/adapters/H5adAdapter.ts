/**
 * H5AD Data Adapter
 * Handles all H5AD-specific data operations and h5wasm interactions
 */

// We'll load h5wasm dynamically from CDN
import h5wasm from "h5wasm";

import { columnPalette, compareLabels } from "../utils/color-palette";
import { isCategorical as detectCategorical } from "../utils/column-type-detection";
import { CONVENTIONS, POSITION_COLUMNS } from "../spatialdata/conventions";

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
  async initialize(
    file: File,
    onProgress?: (progress: number, message: string) => Promise<void> | void,
  ) {
    await onProgress?.(10, "Initializing h5wasm...");
    const { FS } = await h5wasm.ready;

    await onProgress?.(30, "Reading file into memory...");
    const arrayBuffer = await file.arrayBuffer();

    await onProgress?.(50, "Loading HDF5 structure...");
    FS.writeFile(file.name, new Uint8Array(arrayBuffer));
    this.h5File = new h5wasm.File(file.name, "r");

    // Get available keys from different groups
    await onProgress?.(70, "Reading dataset metadata...");
    const rootKeys = this.h5File.keys();
    const obsKeys = this.h5File.get("obs").keys();
    const obsmKeys = this.h5File.get("obsm")?.keys() || [];
    const unsKeys = this.h5File.get("uns")?.keys() || [];
    const varKeys = this.h5File.get("var")?.keys() || [];

    await onProgress?.(90, "Finalizing...");
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
    // A sparse X is a group: its shape lives in the attrs.
    const [numCells, numGenes] = xDataset.shape || this.sparseShape(xDataset);

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
   * A column anndata writes as a group of `values` and a boolean `mask`
   * (true where missing): its nullable-string-array, nullable-integer and
   * nullable-boolean encodings. anndata 0.13 with pandas 3 writes every
   * string column this way, the index included. A missing string is "" (the
   * label a missing category gets too), a missing number NaN, and booleans
   * are labelled as pandas prints them.
   */
  private readMasked(node: any): any[] | null {
    const encoding = String(this.attrValue(node, "encoding-type") ?? "");

    if (!encoding.startsWith("nullable-")) return null;
    const mask = node.get("mask").value as ArrayLike<unknown>;
    const values = node.get("values").value;

    if (encoding === "nullable-string-array") {
      return Array.from(values, (v, i) =>
        mask[i] ? "" : this.decodeBytes(v),
      );
    }
    if (encoding === "nullable-boolean") {
      return Array.from(values as ArrayLike<unknown>, (v, i) =>
        mask[i] ? "" : v ? "True" : "False",
      );
    }

    return Array.from(values as ArrayLike<number>, (v, i) =>
      mask[i] ? NaN : Number(v),
    );
  }

  /** A string column in either encoding: a plain array, or a masked group. */
  private readStrings(node: any): string[] {
    return (
      this.readMasked(node) ??
      Array.from(node.value, this.decodeBytes.bind(this))
    );
  }

  /** The column `obs` or `var` names as its index (`_index` by default). */
  private indexName(group: "obs" | "var"): string {
    const named = this.attrValue(this.h5File.get(group), "_index");

    return typeof named === "string" && named ? named : "_index";
  }

  /**
   * Fetch observation data (cell metadata)
   */
  async fetchObs(group: string): Promise<any[]> {
    const obs = this.h5File.get("obs");

    if (group === "index") {
      return this.readStrings(obs.get(this.indexName("obs")));
    }
    const field = obs.get(group);
    const masked = this.readMasked(field);

    if (masked) return masked;
    const keys = field.keys?.();

    if (keys?.includes("categories") && keys.includes("codes")) {
      const categories = Array.from(
        field.get("categories").value,
        this.decodeBytes.bind(this),
      );
      const codes = Array.from(field.get("codes").value) as number[];

      // Missing category (code -1) → "" — the same label the server pipeline writes.
      return codes.map((i) => categories[i] ?? "");
    }

    // Booleans are stored as an HDF5 enum {FALSE, TRUE}; label them the way
    // pandas stringifies booleans so both pipelines write the same dictionary.
    const members = field.metadata?.enum_type?.members;

    if (members && "FALSE" in members && "TRUE" in members) {
      return Array.from(field.value as ArrayLike<number>, (v) => (v ? "True" : "False"));
    }

    return Array.from(field.value, this.decodeBytes.bind(this));
  }

  /**
   * Fetch variable data (gene metadata)
   */
  async fetchVar(group: string): Promise<any[]> {
    const v = this.h5File.get("var");

    return this.readStrings(
      v.get(group === "index" ? this.indexName("var") : group),
    );
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

    // Keep at most 3 dimensions (the viewer renders 2D/3D; a 50-component PCA
    // is truncated to its first 3, same as the server pipeline).
    const dimensions = dataset.metadata.shape[1];
    const keep = Math.min(dimensions, 3);
    const rows = Math.floor(flat.length / dimensions);

    coordinates = Array.from({ length: rows }, (_, i) =>
      Array.from({ length: keep }, (_, k) => Number(flat[i * dimensions + k])),
    );

    return {
      coordinates: coordinates,
      dimensions: keep,
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
  /**
   * Try to load spatial coordinates from obs columns
   * Checks alternative naming conventions in order
   */
  private loadSpatialFromObs(): {
    coordinates: number[][];
    dimensions: number;
  } | null {
    const obs = this.h5File.get("obs");

    if (!obs) return null;

    const columnNames = obs.keys();

    // Define alternative naming patterns to check in order
    const namingPatterns = [
      { x: "center_x", y: "center_y", z: "center_z" },
      { x: "centerX", y: "centerY", z: "centerZ" },
      { x: "x", y: "y", z: "z" },
      { x: "centroid_x", y: "centroid_y", z: "centroid_z" },
    ];

    // Find first matching pattern
    for (const pattern of namingPatterns) {
      const hasX = columnNames.includes(pattern.x);
      const hasY = columnNames.includes(pattern.y);

      if (hasX && hasY) {
        return this.extractCoordinatesFromObs(pattern.x, pattern.y, pattern.z);
      }
    }

    return null;
  }

  /**
   * Extract and validate coordinates from obs columns
   */
  private extractCoordinatesFromObs(
    xCol: string,
    yCol: string,
    zCol: string,
  ): { coordinates: number[][]; dimensions: number } {
    const obs = this.h5File.get("obs");
    const columnNames = obs.keys();

    // Read x and y columns
    const xData = obs.get(xCol).value as string[] | number[];
    const yData = obs.get(yCol).value as string[] | number[];

    // Check if z column exists
    const hasZ = columnNames.includes(zCol);
    const zData = hasZ ? (obs.get(zCol).value as string[] | number[]) : null;

    const numRows = xData.length;

    // Convert to numbers and track validity
    const xNumbers: (number | null)[] = [];
    const yNumbers: (number | null)[] = [];
    const zNumbers: (number | null)[] = [];

    let validXCount = 0;
    let validYCount = 0;
    let validZCount = 0;
    let nonZeroZCount = 0;

    for (let i = 0; i < numRows; i++) {
      // Parse X
      const xVal = parseFloat(String(xData[i]));
      const xValid = !isNaN(xVal) && isFinite(xVal);

      xNumbers.push(xValid ? xVal : null);
      if (xValid) validXCount++;

      // Parse Y
      const yVal = parseFloat(String(yData[i]));
      const yValid = !isNaN(yVal) && isFinite(yVal);

      yNumbers.push(yValid ? yVal : null);
      if (yValid) validYCount++;

      // Parse Z if exists
      if (zData) {
        const zVal = parseFloat(String(zData[i]));
        const zValid = !isNaN(zVal) && isFinite(zVal);

        zNumbers.push(zValid ? zVal : null);
        if (zValid) {
          validZCount++;
          if (zVal !== 0) nonZeroZCount++;
        }
      }
    }

    // Validate: at least 90% of x and y must be valid
    const xValidPercent = (validXCount / numRows) * 100;
    const yValidPercent = (validYCount / numRows) * 100;

    if (xValidPercent < 90) {
      throw new Error(
        `Invalid x coordinates in obs/${xCol}: only ${xValidPercent.toFixed(1)}% are valid (need ≥90%)`,
      );
    }

    if (yValidPercent < 90) {
      throw new Error(
        `Invalid y coordinates in obs/${yCol}: only ${yValidPercent.toFixed(1)}% are valid (need ≥90%)`,
      );
    }

    // Determine dimensions based on z validity
    let dimensions = 2;
    let useZ = false;

    if (zData) {
      const zValidPercent = (validZCount / numRows) * 100;
      const nonZeroPercent =
        validZCount > 0 ? (nonZeroZCount / validZCount) * 100 : 0;

      // Use z if ≥90% are valid and ≥90% of valid values are non-zero
      if (zValidPercent >= 90 && nonZeroPercent >= 90) {
        dimensions = 3;
        useZ = true;
      }
    }

    // One position per cell, in obs order. A cell without a valid position
    // keeps its row as NaN (it is not drawn): dropping the row would shift
    // every later position onto the wrong cell's clusters and expression,
    // and the store saved from this keeps every row too.
    const coordinates: number[][] = [];

    for (let i = 0; i < numRows; i++) {
      const position = useZ
        ? [xNumbers[i], yNumbers[i], zNumbers[i]]
        : [xNumbers[i], yNumbers[i]];

      // Any coordinate missing means no position at all, not half of one
      coordinates.push(
        position.some((v) => v === null)
          ? position.map(() => NaN)
          : (position as number[]),
      );
    }

    return { coordinates, dimensions };
  }

  loadSpatialCoordinates(): { coordinates: number[][]; dimensions: number } {
    // Try obsm/X_spatial first
    let result = this.loadCoordinates("X_spatial");

    if (result) {
      const dataset = this.h5File.get("obsm/X_spatial");
      const dimensions = dataset.metadata.shape[1];

      return { coordinates: result, dimensions };
    }

    // Try obsm/spatial second
    result = this.loadCoordinates("spatial");
    if (result) {
      const dataset = this.h5File.get("obsm/spatial");
      const dimensions = dataset.metadata.shape[1];

      return { coordinates: result, dimensions };
    }

    // Fallback to obs columns
    const obsResult = this.loadSpatialFromObs();

    if (obsResult) {
      return obsResult;
    }

    // If nothing found, throw error
    throw new Error(
      "No spatial coordinates found. Checked: obsm/X_spatial, obsm/spatial, " +
        "obs columns (center_x/y/z, centerX/Y/Z, x/y/z, centroid_x/y/z)",
    );
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
   * Load gene list: the column `var`'s `_index` attribute names (its index,
   * `gene_ids` in a table indexed by Ensembl id), else var/_index, var/gene,
   * var/genes — the same rule as SpatialDataTableAdapter.loadGenes, so a
   * preview and the store saved from it list the same genes.
   */
  async loadGenes(): Promise<string[]> {
    const possibleColumns = [
      ...new Set([
        this.indexName("var"),
        ...CONVENTIONS.tables.featureIndexColumns,
      ]),
    ];

    for (const columnName of possibleColumns) {
      try {
        console.log(
          `[H5adAdapter] Trying to load genes from var['${columnName}']`,
        );
        const varIndex = await this.fetchVar(columnName);

        if (varIndex && varIndex.length > 0) {
          console.log(
            `[H5adAdapter] ✅ Successfully loaded ${varIndex.length} genes from var['${columnName}']`,
          );

          return varIndex;
        }
      } catch (error) {
        console.warn(
          `[H5adAdapter] Column '${columnName}' not found in var, trying next...`,
        );
      }
    }

    console.error("[H5adAdapter] Failed to load genes from any known column");

    return [];
  }

  /**
   * Determine if data is categorical or numerical
   * Uses centralized detection logic from column-type-detection utility
   */
  isCategoricalData(values: any[], columnName?: string): boolean {
    return detectCategorical(values, columnName);
  }

  /**
   * Parse a single obs column to determine type and palette
   */
  async parseObsColumn(columnName: string): Promise<any> {
    try {
      const values = await this.fetchObs(columnName);

      // All-empty columns are not exported (same rule as the server pipeline).
      const isMissing = (v: any) =>
        v == null || (typeof v === "number" && Number.isNaN(v)) || String(v).trim() === "" || String(v) === "NaN";

      if (values.every(isMissing)) return null;

      // Determine if categorical or numerical (pass columnName for special cases)
      const isCategorical = this.isCategoricalData(values, columnName);

      // Build indexed representation
      const valueToIndex = new Map<string, number>();
      const uniqueValuesList: string[] = [];

      for (let i = 0; i < values.length; i++) {
        const str = String(values[i]);

        if (!valueToIndex.has(str)) {
          valueToIndex.set(str, uniqueValuesList.length);
          uniqueValuesList.push(str);
        }
      }

      const uniqueValues = uniqueValuesList.sort(compareLabels);

      const sortedIndexMap = new Map<string, number>();

      for (let i = 0; i < uniqueValues.length; i++) {
        sortedIndexMap.set(uniqueValues[i], i);
      }

      const IndexArray = uniqueValues.length <= 65535 ? Uint16Array : Uint32Array;
      const valueIndices = new IndexArray(values.length);

      for (let i = 0; i < values.length; i++) {
        valueIndices[i] = sortedIndexMap.get(String(values[i]))!;
      }

      if (isCategorical) {
        const palette = await this.loadClusterPalette(columnName);

        return {
          column: columnName,
          type: "categorical",
          values: [],
          valueIndices: valueIndices,
          palette: palette,
          uniqueValues: uniqueValues,
        };
      } else {
        return {
          column: columnName,
          type: "numerical",
          values: [],
          valueIndices: valueIndices,
          palette: null,
          uniqueValues: uniqueValues,
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

    // Get all non-index, non-coordinate columns (the same fixed list of
    // coordinate names the server pipeline leaves out of obs/).
    const nonIndexColumns = this.metadata.obsKeys.filter(
      (key: string) => !key.startsWith("_") && !POSITION_COLUMNS.has(key),
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
   * Stored category order of a categorical obs column (null when the column
   * isn't categorical). scanpy's `uns[col + "_colors"]` follows this order.
   */
  private fetchObsCategories(group: string): string[] | null {
    const field = this.h5File.get("obs").get(group);
    const keys = field?.keys?.();

    if (keys?.includes("categories") && keys.includes("codes")) {
      return Array.from(field.get("categories").value, this.decodeBytes.bind(this)) as string[];
    }

    return null;
  }

  /**
   * Load color palette for clusters: the shared rule (columnPalette) — the
   * default colours dealt in compareLabels order, replaced by
   * `uns[col + "_colors"]` when the file stores one colour per category —
   * the same rule the SpatialData viewer applies to the store saved from
   * this file, so its colours survive the save.
   */
  async loadClusterPalette(
    clusterColumn: string,
  ): Promise<Record<string, string>> {
    const labels = Array.from(
      new Set((await this.fetchObs(clusterColumn)).map(String)),
    ).sort(compareLabels);
    let colors: unknown[] | null = null;

    try {
      const stored = this.fetchUns(clusterColumn + "_colors");

      colors = stored ? Array.from(stored as ArrayLike<unknown>) : null;
    } catch {
      // No colours stored: the defaults stand
    }

    return columnPalette(labels, {
      colors,
      categories: this.fetchObsCategories(clusterColumn),
    });
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
   * Note: This method is deprecated in favor of using cached matrix in StandardizedDataset
   * Kept for backward compatibility
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
   * Fetch full expression matrix (for caching/bulk processing)
   */
  fetchFullMatrix(): any {
    const m = this.h5File.get("X");

    // Sparse X (anndata csr_matrix / csc_matrix group): expand to the same
    // row-major Float32Array a dense X yields, so gene extraction is uniform.
    if (m && typeof m.keys === "function" && m.keys().includes("indptr")) {
      return this.denseFromSparseGroup(m);
    }

    // Work in float32 like the server pipeline (X is read straight into
    // float32 there): DE means then accumulate the same values on both sides,
    // and a float64 X costs half the memory.
    return this.toFloat32(m.value);
  }

  /**
   * Coerce any h5wasm-returned numeric array to Float32Array.
   *
   * An integer-count X stored as int64 comes back as a BigInt64Array; left
   * as-is, its BigInt values reach downstream numeric code and throw
   * "Cannot mix BigInt and other types". The Float32Array constructor also
   * refuses a BigInt typed array, so those are mapped through Number first.
   */
  private toFloat32(value: any): Float32Array {
    if (value instanceof Float32Array) return value;
    if (value instanceof BigInt64Array || value instanceof BigUint64Array) {
      return Float32Array.from(value, Number);
    }

    return new Float32Array(value);
  }

  /** h5wasm exposes attributes as { value, shape, dtype } records. */
  private attrValue(obj: any, name: string): any {
    const a = obj?.attrs?.[name];

    return a && typeof a === "object" && "value" in a ? a.value : a;
  }

  private sparseShape(group: any): [number, number] {
    const shape = this.attrValue(group, "shape"); // BigInt64Array or number[]

    return shape ? [Number(shape[0]), Number(shape[1])] : [0, 0];
  }

  private denseFromSparseGroup(group: any): Float32Array {
    const [nObs, nGenes] = this.sparseShape(group);
    const data = group.get("data").value as ArrayLike<number>;
    const indices = group.get("indices").value as ArrayLike<number | bigint>;
    const indptr = group.get("indptr").value as ArrayLike<number | bigint>;
    const encoding = String(this.attrValue(group, "encoding-type") ?? "csr_matrix");
    const dense = new Float32Array(nObs * nGenes);

    if (encoding.startsWith("csc")) {
      for (let j = 0; j < nGenes; j++) {
        for (let p = Number(indptr[j]); p < Number(indptr[j + 1]); p++) {
          dense[Number(indices[p]) * nGenes + j] = Number(data[p]);
        }
      }
    } else {
      for (let i = 0; i < nObs; i++) {
        for (let p = Number(indptr[i]); p < Number(indptr[i + 1]); p++) {
          dense[i * nGenes + Number(indices[p])] = Number(data[p]);
        }
      }
    }

    return dense;
  }

  /**
   * Fetch column from matrix
   */
  fetchColumn(matrix: any, column: number): number[] {
    if (Array.isArray(matrix[0])) return matrix.map((row: any) => row[column]);

    const { numCells: n_obs, numGenes: n_genes } = this.getDatasetInfo();

    if (ArrayBuffer.isView(matrix)) {
      if (column >= n_genes) throw new Error("Column index out of bounds");
      const typedArray = matrix as any; // Type assertion for indexing

      return Array.from(
        { length: n_obs },
        (_, i) => typedArray[i * n_genes + column],
      );
    }

    throw new Error("Unsupported matrix format");
  }
}
