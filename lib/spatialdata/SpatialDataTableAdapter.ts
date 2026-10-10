/**
 * SpatialData table adapter
 *
 * Reads the AnnData table of a SpatialData zarr store (`tables/<name>`) with
 * zarrita, straight from the store as written — nothing is converted.
 *
 * Lazy gene expression access works for **dense** and **CSC** X matrices
 * (column slicing on CSC is O(nnz_gene)). A **CSR** X with a CSC copy in
 * `layers/` (e.g. `layers/X_csc`) reads genes from that copy; a CSR X with no
 * such copy is scanned in full for every gene.
 *
 * Store-agnostic: takes any zarrita `AsyncReadable` + a list of all keys in
 * the store. The key list is used to enumerate `obs`/`obsm`/`var` children
 * since `AsyncReadable` has no LIST operation.
 */
import type { SpatialDataStore } from "./store";

// Zarr v3 string arrays need a newer zarrita than the one the h5ad-zarr
// path is pinned to.
import * as zarr from "zarrita-v3";

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

export class SpatialDataTableAdapter {
  // Routing signal used by umap-panel / visualization-panel / load-cluster-column /
  // sync hooks: "local" means "call this adapter on the main thread, don't ship
  // to the S3 worker." Both FileMap-backed and S3-backed zarr adapters use
  // this signal — the call-site routing only needs to know "is this S3
  // chunked or not?", and zarr (whether local or S3) is "not S3 chunked."
  readonly mode = "local" as const;

  store: SpatialDataStore["store"];
  // Keys under the AnnData root, relative to it
  storeKeys: string[];
  // Path of the table inside the store: "tables/<name>/"
  tablePrefix = "";
  private root: zarr.Group<zarr.Readable> | null = null;
  // Where gene expression is read from: X, or a CSC layer when X is CSR
  private expr: zarr.Array<zarr.DataType, zarr.Readable> | SparseMatrix | null =
    null;
  private indptr: ArrayLike<number | bigint> | null = null;

  numCells = 0;
  numGenes = 0;
  spatialDimensions = 2;
  xFormat: XFormat = "missing";

  // Enumerated up front from `storeKeys` so we don't need group.keys()
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

  constructor(store: SpatialDataStore["store"], storeKeys: string[]) {
    this.store = store;
    this.storeKeys = storeKeys;
  }

  /**
   * Open the table, detect X format, enumerate top-level obs/obsm/var columns.
   * Does NOT eagerly read column values — caller picks what to load.
   */
  async initialize(
    onProgress?: (progress: number, message: string) => Promise<void> | void,
  ) {
    await onProgress?.(15, "Opening zarr store...");

    const storeRoot = await zarr.open.v3(this.store as zarr.Readable, {
      kind: "group",
    });

    this.tablePrefix = spatialDataTablePrefix(storeRoot.attrs, this.storeKeys);
    this.storeKeys = this.storeKeys
      .filter((k) => k.startsWith(this.tablePrefix))
      .map((k) => k.slice(this.tablePrefix.length));
    this.root = await zarr.open.v3(storeRoot.resolve(this.tablePrefix), {
      kind: "group",
    });

    await onProgress?.(25, "Reading dataset shape...");

    // Detect X format and shape
    const x = hasNode(this.storeKeys, "X") ? await this.openMatrix("X") : null;

    if (!x) {
      this.xFormat = "missing";
    } else {
      this.expr = x;
      this.xFormat = "format" in x ? x.format : "dense";
      this.numCells = x.shape[0];
      this.numGenes = x.shape[1];
    }

    // A CSR X can't be sliced by gene; use a CSC copy from layers/ if present.
    // The tiling manifest names one that the key list may not include (it is
    // written after the store's metadata is consolidated).
    if (this.xFormat === "csr") {
      const hinted = (storeRoot.attrs.spatial_tiling as any)?.spatialdata
        ?.expression_index?.csc?.layer;
      const layers = new Set(enumerateChildArrays(this.storeKeys, "layers"));

      if (hinted) layers.add(hinted);

      for (const layer of layers) {
        const m = await this.openMatrix(`layers/${layer}`).catch(() => null);

        if (
          m &&
          "format" in m &&
          m.format === "csc" &&
          m.shape[0] === this.numCells &&
          m.shape[1] === this.numGenes
        ) {
          this.expr = m;
          this.xFormat = "csc";
          break;
        }
      }
    }

    await onProgress?.(35, "Enumerating obs / obsm / var columns...");

    this.obsColumns = enumerateChildArrays(this.storeKeys, "obs");
    this.obsmKeys = enumerateChildArrays(this.storeKeys, "obsm");
    this.varColumns = enumerateChildArrays(this.storeKeys, "var");

    // Spatial dimensions: derive from obsm/spatial or obsm/X_spatial shape if present
    const spatialKey = this.obsmKeys.includes("X_spatial")
      ? "X_spatial"
      : this.obsmKeys.includes("spatial")
        ? "spatial"
        : null;

    if (spatialKey) {
      try {
        const arr = await this.openArray(`obsm/${spatialKey}`);
        const shape = arr.shape;

        if (shape && shape.length === 2) {
          this.spatialDimensions = shape[1] >= 3 ? 3 : 2;
        }
      } catch {
        // fall back to default 2
      }
    }

    // Build cluster column metadata (filter out the cell index, and the
    // SpatialData join keys, which identify cells rather than annotate them)
    const joinKeys = [this.root.attrs.instance_key, this.root.attrs.region_key];

    this.clusterColumnNames = this.obsColumns.filter(
      (k) => !k.startsWith("_") && k !== "_index" && !joinKeys.includes(k),
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
    if (!this.root) throw new Error("Adapter not initialized");

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

    const arr = await this.openArray(`obsm/${key}`);
    const chunk = await zarr.get(arr, [null, null]);
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
    if (!this.root) throw new Error("Adapter not initialized");

    const candidates = ["_index", "gene", "genes"];

    for (const candidate of candidates) {
      if (!this.varColumns.includes(candidate)) continue;
      try {
        const list = await this.readColumn(`var/${candidate}`);

        if (list.length > 0) {
          this.genes = list;

          return list;
        }
      } catch (e) {
        console.warn(
          `[SpatialDataTableAdapter] Failed to read var/${candidate}:`,
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
    if (!this.root) throw new Error("Adapter not initialized");
    const out: ClusterColumn[] = [];

    for (const columnName of columns) {
      if (!this.obsColumns.includes(columnName)) continue;
      try {
        const values = await this.readColumn(`obs/${columnName}`);

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
          `[SpatialDataTableAdapter] Failed to load cluster column "${columnName}":`,
          e,
        );
      }
    }

    return out;
  }

  /**
   * Read one embedding (e.g. "umap" or "X_umap") on demand.
   * Return shape matches ChunkedDataAdapter.loadEmbedding so callers
   * (umap-panel, etc.) can destructure `{ name, data }`.
   */
  async loadEmbedding(
    name: string,
  ): Promise<{ name: string; data: number[][] } | null> {
    if (!this.root) throw new Error("Adapter not initialized");

    const exact = this.obsmKeys.includes(name) ? name : null;
    const prefixed = this.obsmKeys.includes(`X_${name}`) ? `X_${name}` : null;
    const key = exact ?? prefixed;

    if (!key) return null;

    // Strip the "X_" prefix in the returned name so it matches the convention
    // the rest of the app uses (dataset.embeddings keyed by "umap", not "X_umap").
    const returnedName = key.startsWith("X_") ? key.slice(2) : key;

    const arr = await this.openArray(`obsm/${key}`);
    const chunk = await zarr.get(arr, [null, null]);
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

    return { name: returnedName, data: out };
  }

  /**
   * Lazy per-gene fetch. Cheap for **dense** and **CSC**; a **CSR** matrix
   * with no CSC layer is scanned in full, so densify it up-front instead.
   */
  async fetchGeneExpression(geneName: string): Promise<number[] | null> {
    if (!this.expr) return null;

    const cached = this.geneExprCache.get(geneName);

    if (cached) return cached;

    const genes = await this.loadGenes();
    const geneIndex = genes.indexOf(geneName);

    if (geneIndex < 0) return null;

    const expr = this.expr;
    const result: number[] = new Array(this.numCells).fill(0);

    if (!("format" in expr)) {
      const chunk = await zarr.get(expr, [null, geneIndex]);
      const data = (chunk as any).data as ArrayLike<number>;

      for (let i = 0; i < data.length; i++) result[i] = Number(data[i]);
    } else {
      if (!this.indptr) {
        this.indptr = (await zarr.get(expr.indptr, [null])).data as ArrayLike<
          number | bigint
        >;
      }
      const indptr = this.indptr;

      if (expr.format === "csc") {
        const start = Number(indptr[geneIndex]);
        const stop = Number(indptr[geneIndex + 1]);

        if (stop > start) {
          const sel = [zarr.slice(start, stop)];
          const [rows, values] = await Promise.all([
            zarr.get(expr.indices, sel),
            zarr.get(expr.data, sel),
          ]);
          const r = rows.data as ArrayLike<number | bigint>;
          const v = values.data as ArrayLike<number | bigint>;

          for (let k = 0; k < r.length; k++)
            result[Number(r[k])] = Number(v[k]);
        }
      } else {
        console.warn(
          "[SpatialDataTableAdapter] fetchGeneExpression called on CSR X — this will scan the full matrix. Densify CSR up-front instead.",
        );
        const [cols, values] = await Promise.all([
          zarr.get(expr.indices, [null]),
          zarr.get(expr.data, [null]),
        ]);
        const c = cols.data as ArrayLike<number | bigint>;
        const v = values.data as ArrayLike<number | bigint>;

        for (let row = 0; row < this.numCells; row++) {
          const stop = Number(indptr[row + 1]);

          for (let k = Number(indptr[row]); k < stop; k++) {
            if (Number(c[k]) === geneIndex) result[row] = Number(v[k]);
          }
        }
      }
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

  private openArray(path: string) {
    return zarr.open.v3(this.root!.resolve(path), { kind: "array" });
  }

  /** Open X or a layer: a dense array, or a csr/csc group. */
  private async openMatrix(
    path: string,
  ): Promise<zarr.Array<zarr.DataType, zarr.Readable> | SparseMatrix> {
    const node = await zarr.open.v3(this.root!.resolve(path));

    if (node instanceof zarr.Array) return node;

    const [indptr, indices, data] = await Promise.all(
      ["indptr", "indices", "data"].map((k) =>
        zarr.open.v3(node.resolve(k), { kind: "array" }),
      ),
    );

    return {
      format: String(node.attrs["encoding-type"]).startsWith("csc")
        ? "csc"
        : "csr",
      shape: (node.attrs.shape ?? node.attrs.h5sparse_shape) as number[],
      indptr,
      indices,
      data,
    };
  }

  /**
   * Read a dataframe column (`obs/<name>`, `var/<name>`) as strings,
   * resolving categoricals to their labels.
   */
  private async readColumn(path: string): Promise<string[]> {
    const node = await zarr.open.v3(this.root!.resolve(path));

    if (node instanceof zarr.Array) {
      return arrayToStringArray((await zarr.get(node, [null])).data);
    }

    const [codes, categories] = await Promise.all(
      ["codes", "categories"].map((k) =>
        zarr.open
          .v3(node.resolve(k), { kind: "array" })
          .then((arr) => zarr.get(arr, [null])),
      ),
    );
    const labels = arrayToStringArray(categories.data);
    const c = codes.data as ArrayLike<number | bigint>;
    const out: string[] = new Array(c.length);

    // Code -1 marks a missing value
    for (let i = 0; i < c.length; i++) out[i] = labels[Number(c[i])] ?? "";

    return out;
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

interface SparseMatrix {
  format: "csr" | "csc";
  shape: number[];
  indptr: zarr.Array<zarr.DataType, zarr.Readable>;
  indices: zarr.Array<zarr.DataType, zarr.Readable>;
  data: zarr.Array<zarr.DataType, zarr.Readable>;
}

/** The path of the store's table (`tables/<name>/`). */
function spatialDataTablePrefix(
  rootAttrs: Record<string, any>,
  storeKeys: string[],
): string {
  if (!rootAttrs.spatialdata_attrs) {
    throw new Error(
      "Not a SpatialData store (no spatialdata_attrs at the root)",
    );
  }

  const tables = enumerateChildArrays(storeKeys, "tables");

  if (tables.length === 0) {
    throw new Error("SpatialData store has no table under tables/");
  }
  const table = tables.includes("table") ? "table" : tables[0];

  return `tables/${table}/`;
}

function hasNode(storeKeys: string[], path: string): boolean {
  return storeKeys.some((k) => k.startsWith(`${path}/`));
}

/**
 * Enumerate immediate child names under a top-level zarr group prefix by
 * inspecting which keys in the store start with `${prefix}/...`.
 *
 * A child is considered an array/group if any of its keys looks like
 * `<prefix>/<child>/.zarray`, `.zgroup`, or `zarr.json`.
 */
function enumerateChildArrays(storeKeys: string[], prefix: string): string[] {
  const found = new Set<string>();
  const markers = [".zarray", ".zgroup", "zarr.json"];

  for (const key of storeKeys) {
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
  // zarrita returns TypedArray, string[] or a fixed-width string array.
  if (
    data &&
    typeof data.length === "number" &&
    typeof data.get === "function"
  ) {
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
