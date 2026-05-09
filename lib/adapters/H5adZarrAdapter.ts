/**
 * H5AD-Zarr Data Adapter
 *
 * Reads AnnData stored in zarr format (h5ad → zarr conversion) using anndata.js.
 *
 * Lazy gene expression access works for **dense** and **CSC** X matrices
 * (column slicing on CSC is O(nnz_gene)). For **CSR** the caller is expected
 * to densify the matrix once up front via `lib/workers/zarr-densify.worker.ts`,
 * cache the result on `dataset.matrix`, and route gene queries through that.
 */
import { readZarr, get, AnnData } from "anndata.js";
import * as zarr from "zarrita";

import { FileMapStore } from "../storage/FileMapStore";
import { DEFAULT_COLOR_PALETTE } from "../utils/color-palette";
import { isCategorical as detectCategorical } from "../utils/column-type-detection";

export type XFormat = "dense" | "csr" | "csc" | "missing";

interface ClusterColumn {
  column: string;
  type: "categorical" | "numerical";
  values: any[];
  valueIndices?: Uint16Array | Uint32Array;
  palette: Record<string, string> | null;
  uniqueValues?: string[];
}

export class H5adZarrAdapter {
  fileMap: Map<string, File>;
  store: FileMapStore;
  adata: AnnData<FileMapStore, any, any> | null = null;

  numCells = 0;
  numGenes = 0;
  spatialDimensions = 2;
  xFormat: XFormat = "missing";

  // Enumerated up front from the file map so we don't need group.keys()
  obsColumns: string[] = [];
  obsmKeys: string[] = [];
  varColumns: string[] = [];

  // Names + types for cluster columns (filtered subset of obsColumns)
  clusterColumnNames: string[] = [];
  clusterColumnTypes: Record<string, string> = {};

  // Cached gene list
  private genes: string[] | null = null;

  // Cache for lazy gene queries (dense/CSC path)
  private geneExprCache = new Map<string, number[]>();

  constructor(fileMap: Map<string, File>) {
    this.fileMap = fileMap;
    this.store = new FileMapStore(fileMap);
  }

  /**
   * Open the zarr root, detect X format, enumerate top-level obs/obsm/var columns.
   * Does NOT eagerly read column values — caller picks what to load.
   */
  async initialize(
    onProgress?: (progress: number, message: string) => Promise<void> | void,
  ) {
    await onProgress?.(15, "Opening zarr store...");

    this.adata = await readZarr(this.store);

    await onProgress?.(25, "Reading dataset shape...");

    // Detect X format and shape
    if (this.adata.X === undefined) {
      this.xFormat = "missing";
    } else if ("format" in this.adata.X) {
      // SparseArray
      this.xFormat = (this.adata.X as any).format as "csr" | "csc";
      const shape = (this.adata.X as any).shape as number[];

      this.numCells = shape[0];
      this.numGenes = shape[1];
    } else {
      // dense zarr.Array
      this.xFormat = "dense";
      const shape = (this.adata.X as any).shape as number[];

      this.numCells = shape[0];
      this.numGenes = shape[1];
    }

    await onProgress?.(35, "Enumerating obs / obsm / var columns...");

    this.obsColumns = enumerateChildArrays(this.fileMap, "obs");
    this.obsmKeys = enumerateChildArrays(this.fileMap, "obsm");
    this.varColumns = enumerateChildArrays(this.fileMap, "var");

    // Spatial dimensions: derive from obsm/spatial or obsm/X_spatial shape if present
    const spatialKey = this.obsmKeys.includes("X_spatial")
      ? "X_spatial"
      : this.obsmKeys.includes("spatial")
        ? "spatial"
        : null;

    if (spatialKey) {
      try {
        const arr = await this.adata.obsm.get(spatialKey);
        const shape = (arr as any).shape as number[];

        if (shape && shape.length === 2) {
          this.spatialDimensions = shape[1] >= 3 ? 3 : 2;
        }
      } catch {
        // fall back to default 2
      }
    }

    // Build cluster column metadata (filter out the cell index)
    this.clusterColumnNames = this.obsColumns.filter(
      (k) => !k.startsWith("_") && k !== "_index",
    );
    // Type detection happens on-demand in loadClusters; default everything to
    // "categorical" pre-load so the picker can show them.
    for (const c of this.clusterColumnNames) {
      this.clusterColumnTypes[c] = "categorical";
    }

    await onProgress?.(45, "Initialization complete");
  }

  /**
   * Read spatial coordinates as a flat Float32Array (numCells × dims, row-major).
   */
  async loadSpatialCoordinates(): Promise<{
    coordinates: Float32Array;
    dimensions: number;
  }> {
    if (!this.adata) throw new Error("Adapter not initialized");

    const key = this.obsmKeys.includes("X_spatial")
      ? "X_spatial"
      : this.obsmKeys.includes("spatial")
        ? "spatial"
        : null;

    if (!key) {
      throw new Error(
        "No spatial coordinates found in zarr (looked for obsm/X_spatial and obsm/spatial)",
      );
    }

    const arr = await this.adata.obsm.get(key);
    const chunk = await get(arr as any, [null, null]);
    const shape = (chunk as any).shape as number[];
    const numRows = shape[0];
    const dims = shape[1];
    const out = new Float32Array(numRows * Math.min(dims, 3));
    const data = (chunk as any).data as ArrayLike<number>;
    const outDims = Math.min(dims, 3);

    for (let i = 0; i < numRows; i++) {
      for (let d = 0; d < outDims; d++) {
        out[i * outDims + d] = Number(data[i * dims + d]);
      }
    }

    this.spatialDimensions = outDims >= 3 ? 3 : 2;

    return { coordinates: out, dimensions: this.spatialDimensions };
  }

  /**
   * Read the gene list from var/_index (or var/gene, var/genes).
   */
  async loadGenes(): Promise<string[]> {
    if (this.genes) return this.genes;
    if (!this.adata) throw new Error("Adapter not initialized");

    const candidates = ["_index", "gene", "genes"];

    for (const candidate of candidates) {
      if (!this.varColumns.includes(candidate)) continue;
      try {
        const arr = await this.adata.var.get(candidate);
        const values = await get(arr as any, [null]);
        const data = (values as any).data;
        const list = arrayToStringArray(data);

        if (list.length > 0) {
          this.genes = list;

          return list;
        }
      } catch (e) {
        console.warn(
          `[H5adZarrAdapter] Failed to read var/${candidate}:`,
          e,
        );
      }
    }

    throw new Error("No gene list found in var/_index, var/gene, or var/genes");
  }

  /**
   * Read one or more obs columns and build ClusterColumn entries with
   * indexed values + palette.
   */
  async loadClusters(columns: string[]): Promise<ClusterColumn[]> {
    if (!this.adata) throw new Error("Adapter not initialized");
    const out: ClusterColumn[] = [];

    for (const columnName of columns) {
      if (!this.obsColumns.includes(columnName)) continue;
      try {
        const arr = await this.adata.obs.get(columnName);
        const chunk = await get(arr as any, [null]);
        const data = (chunk as any).data;
        const values = arrayToStringArray(data);

        const isCategorical = detectCategorical(values, columnName);

        const valueToIndex = new Map<string, number>();
        const uniqueValuesList: string[] = [];

        for (let i = 0; i < values.length; i++) {
          const s = values[i];

          if (!valueToIndex.has(s)) {
            valueToIndex.set(s, uniqueValuesList.length);
            uniqueValuesList.push(s);
          }
        }

        const uniqueValues = uniqueValuesList.sort((a, b) =>
          a.localeCompare(b, undefined, { numeric: true }),
        );
        const sortedMap = new Map<string, number>();

        for (let i = 0; i < uniqueValues.length; i++) {
          sortedMap.set(uniqueValues[i], i);
        }

        const IndexArray =
          uniqueValues.length <= 65535 ? Uint16Array : Uint32Array;
        const valueIndices = new IndexArray(values.length);

        for (let i = 0; i < values.length; i++) {
          valueIndices[i] = sortedMap.get(values[i])!;
        }

        const type: "categorical" | "numerical" = isCategorical
          ? "categorical"
          : "numerical";

        this.clusterColumnTypes[columnName] = type;

        out.push({
          column: columnName,
          type,
          values: [],
          valueIndices,
          palette: isCategorical ? buildPalette(uniqueValues) : null,
          uniqueValues,
        });
      } catch (e) {
        console.warn(
          `[H5adZarrAdapter] Failed to load cluster column "${columnName}":`,
          e,
        );
      }
    }

    return out;
  }

  /**
   * Read one embedding (e.g. "X_umap") as nested number[][].
   */
  async loadEmbedding(name: string): Promise<number[][] | null> {
    if (!this.adata) throw new Error("Adapter not initialized");
    const key = this.obsmKeys.includes(name)
      ? name
      : this.obsmKeys.includes(`X_${name}`)
        ? `X_${name}`
        : null;

    if (!key) return null;

    const arr = await this.adata.obsm.get(key);
    const chunk = await get(arr as any, [null, null]);
    const shape = (chunk as any).shape as number[];
    const data = (chunk as any).data as ArrayLike<number>;
    const rows = shape[0];
    const cols = Math.min(shape[1], 3);
    const out: number[][] = new Array(rows);

    for (let i = 0; i < rows; i++) {
      const row = new Array(cols);

      for (let j = 0; j < cols; j++) {
        row[j] = Number(data[i * shape[1] + j]);
      }
      out[i] = row;
    }

    return out;
  }

  /**
   * Lazy per-gene fetch. Only safe for **dense** and **CSC** X.
   * For CSR, the caller should densify up-front and feed `dataset.matrix`.
   */
  async fetchGeneExpression(geneName: string): Promise<number[] | null> {
    if (!this.adata || !this.adata.X) return null;

    const cached = this.geneExprCache.get(geneName);

    if (cached) return cached;

    const genes = await this.loadGenes();
    const geneIndex = genes.indexOf(geneName);

    if (geneIndex < 0) return null;

    if (this.xFormat === "csr") {
      console.warn(
        "[H5adZarrAdapter] fetchGeneExpression called on CSR X — this will scan the full matrix. Densify CSR up-front instead.",
      );
    }

    const chunk = await get(this.adata.X as any, [null, geneIndex]);
    const data = (chunk as any).data as ArrayLike<number>;
    const result: number[] = new Array(data.length);

    for (let i = 0; i < data.length; i++) {
      result[i] = Number(data[i]);
    }

    this.geneExprCache.set(geneName, result);

    return result;
  }

  fetchFullMatrix(): null {
    return null;
  }

  /**
   * For interface parity with ChunkedDataAdapter.
   * `matrix` arg is unused; we go through fetchGeneExpression's cache path.
   */
  async fetchColumn(_matrix: any, geneIndex: number): Promise<number[]> {
    const genes = await this.loadGenes();
    const gene = genes[geneIndex];

    if (!gene) throw new Error(`Gene index ${geneIndex} out of range`);
    const result = await this.fetchGeneExpression(gene);

    if (!result) throw new Error(`Failed to fetch expression for ${gene}`);

    return result;
  }

  getClusterColumnInfo(): { names: string[]; types: Record<string, string> } {
    return {
      names: this.clusterColumnNames,
      types: { ...this.clusterColumnTypes },
    };
  }

  getDatasetInfo() {
    const availableEmbeddings = this.obsmKeys
      .filter((k) => k !== "X_spatial" && k !== "spatial")
      .map((k) => (k.startsWith("X_") ? k.slice(2) : k));

    return {
      id: undefined as string | undefined,
      name: undefined as string | undefined,
      type: "h5ad-zarr",
      numCells: this.numCells,
      numGenes: this.numGenes,
      spatialDimensions: this.spatialDimensions,
      availableEmbeddings,
      clusterCount: this.clusterColumnNames.length,
      normalized: false,
      xFormat: this.xFormat,
    };
  }
}

/**
 * Enumerate immediate child names under a top-level zarr group prefix
 * by inspecting which paths in the file map start with `${prefix}/...`.
 *
 * A child is considered an array/group if any of its keys looks like
 * `<prefix>/<child>/.zarray`, `.zgroup`, or `zarr.json`.
 */
function enumerateChildArrays(
  fileMap: Map<string, File>,
  prefix: string,
): string[] {
  const found = new Set<string>();
  const markers = [".zarray", ".zgroup", "zarr.json"];

  for (const key of fileMap.keys()) {
    if (!key.startsWith(`${prefix}/`)) continue;
    const rest = key.slice(prefix.length + 1);
    const parts = rest.split("/");

    if (parts.length < 2) continue;
    const child = parts[0];
    const last = parts[parts.length - 1];

    if (markers.includes(last)) {
      // We only mark when the marker is at depth-1 (direct child) or deeper
      // (the deeper case still implies the child group exists).
      found.add(child);
    }
  }

  return Array.from(found).sort();
}

function arrayToStringArray(data: any): string[] {
  // anndata.js's LazyCategoricalArray returns string typed arrays;
  // raw zarrita arrays return TypedArray or string arrays.
  if (data && typeof data.length === "number" && typeof data.get === "function") {
    // zarr.UnicodeStringArray / ByteStringArray
    const out: string[] = new Array(data.length);

    for (let i = 0; i < data.length; i++) out[i] = String(data.get(i));

    return out;
  }
  if (Array.isArray(data) || ArrayBuffer.isView(data)) {
    const arr = data as ArrayLike<unknown>;
    const out: string[] = new Array(arr.length);

    for (let i = 0; i < arr.length; i++) out[i] = String(arr[i]);

    return out;
  }

  return [];
}

function buildPalette(uniqueValues: string[]): Record<string, string> {
  const palette: Record<string, string> = {};

  for (let i = 0; i < uniqueValues.length; i++) {
    palette[uniqueValues[i]] =
      DEFAULT_COLOR_PALETTE[i % DEFAULT_COLOR_PALETTE.length];
  }

  return palette;
}

// silence unused-import warnings if tree-shake removes them
void zarr;
