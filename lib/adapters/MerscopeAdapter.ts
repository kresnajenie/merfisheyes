// src/data/adapters/MerscopeAdapter.ts
import { fileToTextMaybeGz } from "@/lib/utils/gzip";
import Papa, { ParseResult } from "papaparse";

interface MerscopeMetadata {
  obsKeys: string[];
  clusterColumn: string | null;
  hasUMAP: boolean;
  hasGenes: boolean;
}

interface SpatialCoordinates {
  coordinates: number[][];
  dimensions: number;
}

interface ClusterData {
  column: string;
  values: string[];
  palette: Record<string, string>;
}

type RowData = Record<string, any>;

export class MerscopeAdapter {
  files: File[];
  metadata: MerscopeMetadata;
  _rows: RowData[];
  _obsKeys: string[];
  _clusterColumn: string | null;
  _umap: SpatialCoordinates | null;
  _genes: string[];
  _exprByGene: Map<string, Float32Array>;
  _cellIndex: Map<string, number>;
  _onProgress?: (progress: number, message: string) => Promise<void> | void;

  constructor() {
    this.files = [];
    this.metadata = {
      obsKeys: [],
      clusterColumn: null,
      hasUMAP: false,
      hasGenes: false,
    };

    this._rows = [];
    this._obsKeys = [];
    this._clusterColumn = null;
    this._umap = null;

    this._genes = [];
    this._exprByGene = new Map();
    this._cellIndex = new Map();
  }

  async initialize(files: File[], onProgress?: (progress: number, message: string) => Promise<void> | void): Promise<void> {
    this._onProgress = onProgress;
    this.files = files;

    // 1) cell_metadata.csv (positions live here)
    await this._onProgress?.(10, "Loading cell metadata...");
    const metaRows = await this._readTableOneOf(["cell_metadata.csv"]);
    if (!metaRows.length) {
      throw new Error("MERSCOPE: cell_metadata.csv not found or empty");
    }

    // 2) cell_categories.csv (cluster labels like "leiden")
    await this._onProgress?.(25, "Loading cluster categories...");
    const catRows = await this._readTableOneOf(["cell_categories.csv"]);
    let clusterCol = null;
    let catIndex = null;
    if (catRows.length) {
      // choose the ID key that is present in categories
      const catIdKey = firstPresent(Object.keys(catRows[0]), [
        "EntityID",
        "cell",
        "cell_id",
        "id",
      ]);
      // choose a good label column (e.g., "leiden")
      clusterCol = firstPresent(Object.keys(catRows[0]), [
        "leiden",
        "cluster",
        "clusters",
        "cell_type",
        "celltype",
        "annotation",
        "annotations",
        "class",
        "subclass",
        "label",
        "labels",
      ]);
      if (catIdKey && clusterCol) {
        catIndex = new Map();
        for (const r of catRows) {
          const id = String(r[catIdKey] ?? "");
          const lab = r[clusterCol] ?? "";
          if (id) catIndex.set(id, lab);
        }
      }
    }

    // 3) cell_numeric_categories.csv (UMAP)
    await this._onProgress?.(40, "Loading UMAP embeddings...");
    const numRows = await this._readTableOneOf([
      "cell_numeric_categories.csv",
    ]);
    let umapIndex = null;
    if (numRows.length) {
      const idKey = firstPresent(Object.keys(numRows[0]), [
        "EntityID",
        "cell",
        "cell_id",
        "id",
      ]);
      const ux = firstPresent(Object.keys(numRows[0]), [
        "umap_X",
        "umap_x",
        "UMAP_X",
        "umap1",
        "x_umap",
      ]);
      const uy = firstPresent(Object.keys(numRows[0]), [
        "umap_Y",
        "umap_y",
        "UMAP_Y",
        "umap2",
        "y_umap",
      ]);
      if (idKey && ux && uy) {
        umapIndex = new Map();
        for (const r of numRows) {
          const id = String(r[idKey] ?? "");
          const x = toNum(r[ux]);
          const y = toNum(r[uy]);
          if (id && Number.isFinite(x) && Number.isFinite(y)) {
            umapIndex.set(id, [x, y]);
          }
        }
      }
    }

    // 4) Merge rows — choose a stable ID present in metadata
    await this._onProgress?.(55, "Merging cell data...");
    const metaIdKey = firstPresent(Object.keys(metaRows[0]), [
      "EntityID",
      "cell",
      "cell_id",
      "id",
    ]);
    if (!metaIdKey) {
      throw new Error(
        "MERSCOPE: could not find an ID column in cell_metadata.csv"
      );
    }

    this._rows = metaRows.map((r) => {
      const id = String(r[metaIdKey] ?? "");
      const out = { ...r };

      // attach cluster label if available
      if (catIndex && clusterCol) {
        out[clusterCol] = catIndex.get(id) ?? "";
      }

      // attach umap for convenience (we’ll also expose as an embedding)
      if (umapIndex) {
        const u = umapIndex.get(id);
        if (u) {
          out._umap_x = u[0];
          out._umap_y = u[1];
        }
      }

      return out;
    });

    // Build cell index (id -> row idx) for expression mapping
    this._cellIndex.clear();
    for (let i = 0; i < this._rows.length; i++) {
      const r = this._rows[i];
      const id = String(r[metaIdKey] ?? "");
      if (id) this._cellIndex.set(id, i);
    }

    // 5) genes + expression — from cell_by_gene.csv (supports wide or long)
    await this._onProgress?.(70, "Loading gene expression...");
    try {
      const cbgRows = await this._readTableOneOf(["cell_by_gene.csv"]);
      if (cbgRows.length) {
        await this._ingestCellByGene(cbgRows, metaIdKey);
      } else {
        console.warn(
          "[MerscopeAdapter] cell_by_gene.csv not found; gene coloring disabled."
        );
      }
    } catch (e) {
      console.warn("[MerscopeAdapter] failed reading cell_by_gene.csv:", e);
    }

    // 6) remember obs keys and cluster column
    await this._onProgress?.(85, "Processing clusters and embeddings...");
    this._obsKeys = Array.from(new Set(Object.keys(this._rows[0] || {})));
    const clusterCandidates = [
      "leiden",
      "cluster",
      "clusters",
      "cluster_id",
      "clusterid",
      "cluster label",
      "cluster_label",
      "cluster id",
      "cell_type",
      "cell type",
      "cell-type",
      "celltype",
      "celltypes",
      "cell types",
      "cell_class",
      "cell class",
      "cell_subclass",
      "cell subclass",
      "predicted_celltype",
      "predicted celltype",
      "predicted cell type",
      "annotation",
      "annotations",
      "cell_annotation",
      "cell annotation",
      "cell_annotations",
      "cell annotations",
      "class",
      "class_label",
      "class label",
      "subclass",
      "subclass_label",
      "subclass label",
      "label",
      "labels",
      "group",
      "group_label",
      "group label",
    ];
    this._clusterColumn =
      clusterCandidates.find((k) => this._obsKeys.includes(k)) || null;

    if (!this._clusterColumn) {
      const heuristicCluster = detectLikelyClusterColumn(this._rows);
      if (heuristicCluster) {
        this._clusterColumn = heuristicCluster;
        if (!this._obsKeys.includes(heuristicCluster)) {
          this._obsKeys.push(heuristicCluster);
        }
        console.log(
          `[MerscopeAdapter] Auto-selected "${heuristicCluster}" as cluster column via heuristic detection.`
        );
      }
    }

    if (!this._clusterColumn) {
      console.warn(
        '[MerscopeAdapter] No celltype-like column found; will fall back to single "All" group.'
      );
    }

    // 7) build UMAP embedding map if present on rows
    const haveU =
      this._rows.length &&
      "_umap_x" in this._rows[0] &&
      "_umap_y" in this._rows[0];
    if (haveU) {
      const coords = [];
      for (const r of this._rows) {
        const x = toNum(r._umap_x),
          y = toNum(r._umap_y);
        if (Number.isFinite(x) && Number.isFinite(y)) coords.push([x, y]);
      }
      this._umap = { coordinates: coords, dimensions: 2 };
    }

    // metadata summary
    this.metadata = {
      obsKeys: this._obsKeys,
      clusterColumn: this._clusterColumn,
      hasUMAP: !!this._umap,
      hasGenes: !!this._genes.length,
    };

    console.log(
      `[MerscopeAdapter] Loaded ${this._rows.length} cells, ${this._genes.length} genes.`
    );
  }

  getDatasetInfo(): { numCells: number; numGenes: number } {
    return {
      numCells: this._rows.length,
      numGenes: this._genes.length,
    };
  }

  // --- spatial coords from cell_metadata: center_x/center_y (robust)
  loadSpatialCoordinates(): SpatialCoordinates {
    const rows = this._rows;
    if (!rows.length) return { coordinates: [], dimensions: 2 };

    const xKey = firstPresent(Object.keys(rows[0]), [
      "center_x",
      "centroid_x",
      "x",
      "x_centroid",
      "x_px",
      "x_position",
    ]);
    const yKey = firstPresent(Object.keys(rows[0]), [
      "center_y",
      "centroid_y",
      "y",
      "y_centroid",
      "y_px",
      "y_position",
    ]);

    if (!xKey || !yKey) {
      console.warn(
        "[MerscopeAdapter] Could not detect centroid columns; returning 0 points"
      );
      return { coordinates: [], dimensions: 2 };
    }

    const coords = [];
    let dropped = 0;
    for (const r of rows) {
      const x = toNum(r[xKey]);
      const y = toNum(r[yKey]);
      if (Number.isFinite(x) && Number.isFinite(y)) coords.push([x, y]);
      else dropped++;
    }
    if (dropped)
      console.warn(
        `[MerscopeAdapter] Dropped ${dropped} rows with invalid centroids`
      );

    return { coordinates: coords, dimensions: 2 };
  }

  loadEmbeddings(): Record<string, number[][]> {
    if (this._umap) return { umap: this._umap.coordinates };
    return {};
  }

  async loadGenes(): Promise<string[]> {
    return this._genes;
  }

  async loadClusters(): Promise<ClusterData> {
    if (!this._clusterColumn) {
      // fabricate single group "All"
      return {
        column: "cluster",
        values: new Array(this._rows.length).fill("All"),
        palette: { All: "#1f77b4" },
      };
    }
    const clusterColumn = this._clusterColumn; // Store in non-null variable
    const vals = this._rows.map((r) => String(r[clusterColumn] ?? ""));
    const uniq = Array.from(new Set(vals)).sort();
    const palette: Record<string, string> = {};
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
    uniq.forEach((u, i) => {
      palette[u] = defaultColors[i % defaultColors.length];
    });
    return { column: clusterColumn, values: vals, palette };
  }

  // ========= StandardizedDataset adapter surface =========
  async fetchObs(column: string): Promise<any[]> {
    if (!column || !this._rows.length) return [];
    return this._rows.map((r) => r[column]);
  }

  async fetchVar(_column: string): Promise<any[]> {
    return [];
  }

  fetchUns(_key: string): any {
    return null;
  }

  // NEW: return expression vector for a gene (array length == number of cells)
  async fetchGeneExpression(gene: string): Promise<number[] | null> {
    if (!gene || !this._exprByGene.has(gene)) return null;
    return Array.from(this._exprByGene.get(gene)!);
  }

  /**
   * Fetch full expression matrix (for caching/bulk processing)
   * Returns gene expression data organized by gene
   */
  fetchFullMatrix(): Map<string, Float32Array> {
    return this._exprByGene;
  }

  /**
   * Fetch column from matrix by gene index
   */
  fetchColumn(matrix: Map<string, Float32Array>, geneIndex: number): number[] {
    const gene = this._genes[geneIndex];
    if (!gene || !matrix.has(gene)) return [];
    return Array.from(matrix.get(gene)!);
  }

  getObsColumns(): string[] {
    return this._obsKeys;
  }

  getObsmEmbeddings(): string[] {
    const keys = ["X_spatial"];
    if (this._umap) keys.push("X_umap");
    return keys;
  }

  private async _parseCsvFile(file: File): Promise<ParsedTable> {
    const name = (file.name || "").toLowerCase();
    if (name.endsWith(".gz")) {
      const text = await fileToTextMaybeGz(file);
      if (!text || !text.trim()) return { headers: [], rows: [] };
      return parseCsvWithPapa(text);
    }
    return parseCsvWithPapa(file);
  }

  private async _readTableOneOf(names: string[]): Promise<RowData[]> {
    const map = new Map(
      this.files.map((f) => [(f.webkitRelativePath || f.name).toLowerCase(), f])
    );
    let target: File | null = null;
    for (const n of names) {
      const needle = n.toLowerCase();
      const candidate = Array.from(map.values()).find((file) => {
        const base = file.name.toLowerCase();
        const rel = (file.webkitRelativePath || "").toLowerCase();
        return base === needle || rel.endsWith("/" + needle);
      });
      if (candidate) {
        target = candidate;
        break;
      }
    }
    if (!target) return [];
    const { rows } = await this._parseCsvFile(target);
    return rows;
  }

  /* ----------------- gene ingestion helpers ----------------- */

  async _ingestCellByGene(rows: RowData[], _metaIdKey: string): Promise<void> {
    if (!rows.length) return;

    const first = rows[0];
    const keys = Object.keys(first);
    const lower = keys.map((k) => k.toLowerCase());

    // Detect LONG format: explicit cell + gene + count
    const hasLong =
      lower.includes("cell") &&
      lower.includes("gene") &&
      (lower.includes("count") ||
        lower.includes("expression") ||
        lower.includes("value"));

    if (hasLong) {
      const cellKey = keys[lower.indexOf("cell")];
      const geneKey = keys[lower.indexOf("gene")];
      const countKey = lower.includes("count")
        ? keys[lower.indexOf("count")]
        : lower.includes("expression")
          ? keys[lower.indexOf("expression")]
          : lower.includes("value")
            ? keys[lower.indexOf("value")]
            : null;

      if (!countKey) {
        console.warn(
          "[MerscopeAdapter] long cell_by_gene.csv has no count-like column; skipping."
        );
        return;
      }

      // gather genes
      const gset = new Set<string>();
      for (const r of rows) {
        const g = String(r[geneKey] ?? "").trim();
        if (g) gset.add(g);
      }
      this._genes = Array.from(gset).sort();

      // init vectors per gene
      const N = this._rows.length;
      for (const g of this._genes) this._exprByGene.set(g, new Float32Array(N));

      // fill
      for (const r of rows) {
        const id = String(r[cellKey] ?? "");
        const g = String(r[geneKey] ?? "");
        const v = toNumBool(r[countKey]); // handles 0/1/"TRUE"/"FALSE"
        const idx = this._cellIndex.get(id);
        if (idx == null || !this._exprByGene.has(g)) continue;
        const vec = this._exprByGene.get(g);
        if (vec) vec[idx] = Number.isFinite(v) ? v : 0;
      }

      console.log(
        `[MerscopeAdapter] Ingested LONG cell_by_gene: genes=${this._genes.length}`
      );
      return;
    }

    // Otherwise: assume WIDE format (first column is cell id, rest are genes)
    const idColGuess =
      firstPresent(keys, ["cell", "EntityID", "cell_id", "id"]) || keys[0];
    const geneCols = keys.filter((k) => k !== idColGuess && !k.startsWith("_"));
    this._genes = geneCols.slice(); // keep order as in file

    const N = this._rows.length;
    for (const g of this._genes) this._exprByGene.set(g, new Float32Array(N));

    for (const r of rows) {
      const id = String(r[idColGuess] ?? "");
      const idx = this._cellIndex.get(id);
      if (idx == null) continue;
      for (const g of this._genes) {
        const v = toNumBool(r[g]);
        const vec = this._exprByGene.get(g);
        if (vec) vec[idx] = Number.isFinite(v) ? v : 0;
      }
    }

    console.log(
      `[MerscopeAdapter] Ingested WIDE cell_by_gene: genes=${this._genes.length}`
    );
  }
}

/* ----------------- small helpers ----------------- */
function toNum(v: any): number {
  if (v === "" || v == null) return NaN;
  const n = Number(v);
  return Number.isFinite(n) ? n : NaN;
}

// Accepts 0/1, numeric strings, TRUE/FALSE, Yes/No as 1/0
function toNumBool(v: any): number {
  if (v == null) return 0;
  if (typeof v === "number") return Number.isFinite(v) ? v : 0;
  const s = String(v).trim();
  if (s === "") return 0;
  const n = Number(s);
  if (Number.isFinite(n)) return n;
  const l = s.toLowerCase();
  if (["true", "yes", "on", "present", "detected"].includes(l)) return 1;
  if (["false", "no", "off", "absent", "undetected"].includes(l)) return 0;
  return 0;
}

function firstPresent(keys: string[], list: string[]): string | null {
  const set = new Set(keys.map((k) => String(k)));
  for (const cand of list) if (set.has(cand)) return cand;
  // also try case-insensitive
  const lowerMap = new Map(keys.map((k) => [String(k).toLowerCase(), k]));
  for (const cand of list) {
    const k = lowerMap.get(String(cand).toLowerCase());
    if (k) return k;
  }
  return null;
}

type ParsedTable = { headers: string[]; rows: RowData[] };

async function parseCsvWithPapa(input: File | string): Promise<ParsedTable> {
  if (typeof input === "string" && !input.trim()) {
    return { headers: [], rows: [] };
  }

  return new Promise<ParsedTable>((resolve, reject) => {
    const rows: RowData[] = [];
    let headers: string[] = [];

    const pushRows = (data: RowData[] | undefined) => {
      if (!data?.length) return;
      for (const row of data) {
        if (!row) continue;
        if (!headers.length) {
          headers = extractHeadersFromRow(row);
        }
        const normalized = normalizeRow(row, headers);
        if (!hasNonEmptyValue(normalized, headers)) continue;
        rows.push(normalized);
      }
    };

    Papa.parse<RowData>(input as any, {
      header: true,
      skipEmptyLines: "greedy",
      worker: false,
      dynamicTyping: false,
      transformHeader: (header: string | undefined) =>
        (header ?? "").trim(),
      chunk: (results: ParseResult<RowData>) => {
        if (!results) return;
        headers = ensureHeaders(headers, results.meta?.fields);
        pushRows(results.data as RowData[]);
      },
      complete: (results: ParseResult<RowData>) => {
        if (!results) {
          resolve({ headers, rows });
          return;
        }
        headers = ensureHeaders(headers, results.meta?.fields);
        if (!rows.length && Array.isArray(results.data)) {
          pushRows(results.data as RowData[]);
        }
        if (results.errors?.length) {
          console.warn(
            "[MerscopeAdapter] PapaParse completed with errors:",
            results.errors
          );
        }
        resolve({ headers, rows });
      },
      error: (error: unknown) => reject(error),
    });
  });
}

function normalizeRow(row: RowData, headers: string[]): RowData {
  if (!row) return {};
  const keys = headers.length ? headers : Object.keys(row);
  for (const key of keys) {
    if (!key) continue;
    if (!(key in row) || row[key] == null) row[key] = "";
  }
  return row;
}

function hasNonEmptyValue(row: RowData, headers: string[]): boolean {
  if (!row) return false;
  const keys = headers.length ? headers : Object.keys(row);
  for (const key of keys) {
    if (!key) continue;
    const value = row[key];
    if (value == null) continue;
    if (typeof value === "number") {
      if (Number.isFinite(value)) return true;
      continue;
    }
    if (typeof value === "boolean") {
      if (value) return true;
      continue;
    }
    const s = String(value).trim();
    if (s !== "") return true;
  }
  return false;
}

function ensureHeaders(
  existing: string[],
  fields?: (string | undefined)[]
): string[] {
  if (existing.length || !fields?.length) return existing;
  const seen = new Set<string>();
  const trimmed: string[] = [];
  for (const h of fields) {
    const header = (h ?? "").trim();
    if (!header || seen.has(header)) continue;
    seen.add(header);
    trimmed.push(header);
  }
  return trimmed;
}

function extractHeadersFromRow(row: RowData): string[] {
  const seen = new Set<string>();
  const headers: string[] = [];
  for (const key of Object.keys(row || {})) {
    const trimmed = key.trim();
    if (!trimmed || seen.has(trimmed)) continue;
    seen.add(trimmed);
    headers.push(trimmed);
  }
  return headers;
}

function detectLikelyClusterColumn(rows: RowData[]): string | null {
  if (!rows.length) return null;

  const sampleKeys = Object.keys(rows[0] || {});
  if (!sampleKeys.length) return null;

  const skipPrefixes = [
    "x",
    "y",
    "z",
    "row",
    "col",
    "column",
    "coord",
    "n_",
    "sum",
    "total",
    "count",
    "umi",
    "umis",
    "reads",
    "intensity",
    "area",
    "volume",
  ];
  const skipRegex = new RegExp(`^(${skipPrefixes.join("|")})`, "i");

  let bestKey: string | null = null;
  let bestScore = Number.POSITIVE_INFINITY;

  for (const key of sampleKeys) {
    if (!key) continue;
    if (key.startsWith("_")) continue;
    const lower = key.toLowerCase();
    if (skipRegex.test(lower)) continue;
    if (lower.includes("gene") || lower.includes("umi")) continue;

    let nonEmpty = 0;
    let numericish = 0;
    const seen = new Set<string>();
    for (const row of rows) {
      const raw = row?.[key];
      if (raw == null || raw === "") continue;
      const str = String(raw).trim();
      if (!str) continue;
      nonEmpty++;
      if (!Number.isNaN(Number(str))) numericish++;
      seen.add(str);
      if (seen.size > 200) break;
    }

    if (seen.size < 2) continue;
    if (nonEmpty === 0) continue;
    if (seen.size > Math.max(50, 0.2 * rows.length)) continue;
    if (numericish / Math.max(nonEmpty, 1) > 0.25) continue;

    const score = seen.size + numericish / Math.max(nonEmpty, 1);
    if (score < bestScore) {
      bestScore = score;
      bestKey = key;
    }
  }

  return bestKey;
}
