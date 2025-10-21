// /src/data/adapters/XeniumAdapter.ts
import { fileToTextMaybeGz } from "@/lib/utils/gzip";
import Papa, { ParseResult } from "papaparse";

interface XeniumMetadata {
  hasPolygons: boolean;
  hasFeatures: boolean;
  obsKeys: string[];
  clusterColumn: string | null;
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

export class XeniumAdapter { 
  files: File[];
  metadata: XeniumMetadata;
  _rows: RowData[];
  _obsKeys: string[];
  _clusterColumn: string | null;
  _genes: string[];
  _exprByGene: Map<string, Float32Array>;
  _cellIndex: Map<string, number>;
  _onProgress?: (progress: number, message: string) => Promise<void> | void;

  constructor() {
    this.files = [];
    this.metadata = {
      hasPolygons: false,
      hasFeatures: false,
      obsKeys: [],
      clusterColumn: null,
    };

    this._rows = [];
    this._obsKeys = [];
    this._clusterColumn = null;

    this._genes = [];
    this._exprByGene = new Map();
    this._cellIndex = new Map();
  }

  async initialize(files: File[], onProgress?: (progress: number, message: string) => Promise<void> | void): Promise<void> {
    this._onProgress = onProgress;
    this.files = files;

    // 1) Load cells table (csv or csv.gz)
    await this._onProgress?.(10, "Loading cell metadata...");
    this._rows = await this._readTableOneOf(["cells.csv", "cells.csv.gz"]);
    if (!this._rows.length)
      throw new Error("Xenium: cells.csv(.gz) not found or empty");

    // 2) Save obs keys
    this._obsKeys = Object.keys(this._rows[0] || {});

    // Build index: choose a best cell id key and map id -> row idx
    await this._onProgress?.(20, "Building cell index...");
    const cellIdKey =
      pickFirstPresent(this._rows[0], [
        "cell_id",
        "id",
        "barcode",
        "cell",
        "cellid",
        "cell_id_original",
      ]) || this._obsKeys[0];
    this._cellIndex.clear();
    for (let i = 0; i < this._rows.length; i++) {
      const cid = String(this._rows[i]?.[cellIdKey] ?? "");
      if (cid) this._cellIndex.set(cid, i);
    }

    // 3) Try to augment clusters from analysis/clustering CSVs
    await this._onProgress?.(30, "Loading cluster data...");
    await this._augmentClustersFromAnalysis(cellIdKey);

    // 4) If still no cluster column, look for common ones in cells.csv
    if (!this._clusterColumn) {
      this._clusterColumn = pickFirstPresent(this._rows[0], [
        "cell_type",
        "celltype",
        "cell type",
        "cell-type",
        "celltypes",
        "cell types",
        "predicted_celltype",
        "predicted celltype",
        "predicted cell type",
        "cluster",
        "clusters",
        "leiden",
        "louvain",
        "cell.types",
        "cell_types",
        "annotation",
        "annotations",
        "class",
        "class_label",
        "class label",
        "subclass",
        "subclass_label",
        "subclass label",
        "label",
        "labels",
        "labels",
        "annotation",
        "annotations",
        "cell_annotation",
        "cell annotation",
        "cell_annotations",
        "cell annotations",
        "cell_class",
        "cell class",
        "cell_subclass",
        "cell subclass",
        "cluster_label",
        "cluster label",
        "cluster_id",
        "cluster id",
        "group",
        "group_label",
        "group label",
      ]);
    }

    // Guard: if the chosen cluster column is basically unique per cell, ignore it
    if (this._clusterColumn) {
      const clusterColumn = this._clusterColumn;
      const vals = this._rows.map((r) => String(r[clusterColumn] ?? ""));
      const uniqCount = new Set(vals.filter((v) => v !== "")).size;
      const n = this._rows.length;
      if (uniqCount > Math.max(50, 0.1 * n)) {
        // too many unique labels => not a real cluster column, likely an ID
        console.warn(
          '[XeniumAdapter] Cluster column looked per-cell unique; falling back to single "All" group.'
        );
        this._clusterColumn = null;
      }
    }
    if (!this._clusterColumn) {
      console.warn(
        'No celltype-like column found; will use single "All" group.'
      );
    }

    // 5) Genes: features, else transcripts
    //    (we'll still build expression from transcripts or cells-wide if possible)
    await this._onProgress?.(50, "Loading gene features...");
    let feats: RowData[] = [];
    try {
      feats = await this._readTableOneOf([
        "features.tsv",
        "features.tsv.gz",
        "features.csv",
        "features.csv.gz",
      ]);
    } catch {}
    if (feats.length) {
      const key = pickFirstPresent(feats[0], [
        "gene",
        "gene_symbol",
        "feature_name",
        "name",
        "target",
        "id",
      ]);
      if (key)
        this._genes = feats.map((r) => String(r[key] || "")).filter(Boolean);
    }

    // 6) Expression: try transcripts.csv(.gz) (long format) first
    await this._onProgress?.(70, "Loading gene expression...");
    let gotExpr = false;
    try {
      const tx = await this._readTableOneOf([
        "transcripts.csv",
        "transcripts.csv.gz",
      ]);
      if (tx.length) {
        await this._ingestFromTranscripts(tx, cellIdKey);
        gotExpr = this._genes.length > 0 && this._exprByGene.size > 0;
      }
    } catch {}

    // 7) Fallback: derive genes/expr from cells.csv if it contains wide gene columns
    if (!gotExpr) {
      const maybeWide = detectWideGeneColumns(this._rows[0], this._genes);
      if (maybeWide.geneCols.length) {
        await this._ingestFromCellsWide(maybeWide.geneCols);
        gotExpr = this._genes.length > 0 && this._exprByGene.size > 0;
      }
    }

    if (!this._genes.length)
      console.warn("No gene names found (or none usable for expression).");

    await this._onProgress?.(90, "Finalizing dataset...");
    this.metadata = {
      hasPolygons: false,
      hasFeatures: !!this._genes.length,
      obsKeys: this._obsKeys,
      clusterColumn: this._clusterColumn,
    };
  }

  getDatasetInfo(): { numCells: number; numGenes: number } {
    return { numCells: this._rows.length, numGenes: this._genes.length };
  }

  // --- Spatial coordinates (robust centroid detection) ---
  loadSpatialCoordinates(): SpatialCoordinates {
    const rows = this._rows;
    if (!rows.length) return { coordinates: [], dimensions: 2 };

    const { xKey, yKey } = detectCentroidKeys(rows[0]);
    if (!xKey || !yKey) {
      console.warn(
        "[XeniumAdapter] Could not detect centroid columns; returning 0 points"
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
        `[XeniumAdapter] Dropped ${dropped} rows with invalid centroids (from ${rows.length})`
      );

    return { coordinates: coords, dimensions: 2 };
  }

  loadEmbeddings(): Record<string, number[][]> {
    return {};
  }

  async loadGenes(): Promise<string[]> {
    return this._genes;
  }

  async loadClusters(): Promise<ClusterData> {
    if (!this._clusterColumn) {
      // fabricate single cluster "All"
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
      "#393b79",
      "#637939",
      "#8c6d31",
      "#843c39",
      "#7b4173",
      "#3182bd",
      "#31a354",
      "#756bb1",
      "#636363",
      "#e6550d",
    ];
    uniq.forEach((u, i) => {
      palette[u] = defaultColors[i % defaultColors.length];
    });
    return { column: clusterColumn, values: vals, palette };
  }

  // ---- obs/var interface expected by StandardizedDataset ----
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
    return ["X_spatial"];
  }

  /* ========================================================
   * Internal helpers
   * ====================================================== */

  async _augmentClustersFromAnalysis(cellIdKey: string): Promise<void> {
    // find analysis/clustering CSVs
    const candidates = this.files.filter((f) => {
      const p = (f.webkitRelativePath || f.name || "").toLowerCase();
      return /analysis\/clustering\//.test(p) && /\.(csv|tsv)(\.gz)?$/.test(p);
    });
    if (!candidates.length) return;

    for (const f of candidates) {
      try {
        const { headers, rows } = await this._parseCsvFile(f);
        if (!rows.length || !headers.length) continue;

        const cellIdKeyAlt = pickFirstPresent(headersObj(headers), [
          "cell_id",
          "cell",
          "id",
          "barcode",
          "cellid",
          "cell_id_original",
        ]);
        const labelKey = pickFirstPresent(headersObj(headers), [
          "cluster",
          "clusters",
          "cell_type",
          "celltype",
          "annotation",
          "annotations",
          "label",
          "labels",
          "group",
          "subclass",
          "class",
        ]);
        if (!cellIdKeyAlt || !labelKey) continue;

        // build map of id->label
        const lut = new Map();
        for (const r of rows) {
          const cid = String(r[cellIdKeyAlt] ?? "").trim();
          const lab = String(r[labelKey] ?? "").trim();
          if (cid) lut.set(cid, lab);
        }

        let assigned = 0;
        for (const r of this._rows) {
          const cid = String(r[cellIdKey] ?? "").trim();
          if (!cid) continue;
          const lab = lut.get(cid);
          if (lab != null) {
            r[labelKey] = lab;
            assigned++;
          }
        }

        if (assigned > 0) {
          this._clusterColumn = labelKey;
          if (!this._obsKeys.includes(labelKey)) this._obsKeys.push(labelKey);
          console.log(
            `[XeniumAdapter] Joined ${assigned} labels from ${f.name} onto cells (column "${labelKey}")`
          );
          return;
        }
      } catch (e) {
        console.warn(
          "[XeniumAdapter] Failed reading clustering file:",
          f.name,
          e
        );
      }
    }
  }

  async _ingestFromTranscripts(rows: RowData[], cellIdKey: string): Promise<void> {
    // Expect long format with columns like: cell_id, gene (or feature_name/target), and per-transcript lines
    const sample = rows[0];
    const geneKey = pickFirstPresent(sample, [
      "gene",
      "feature_name",
      "target",
      "name",
      "id",
    ]);
    const cellKey =
      pickFirstPresent(sample, [
        cellIdKey,
        "cell_id",
        "cell",
        "id",
        "barcode",
        "cellid",
      ]) || cellIdKey;

    if (!geneKey || !cellKey) {
      console.warn(
        "[XeniumAdapter] transcripts.csv lacks recognizable cell/gene columns; skipping."
      );
      return;
    }

    // Build gene set
    const gset = new Set(this._genes); // keep any from features, but add if needed
    for (const r of rows) {
      const g = String(r[geneKey] ?? "").trim();
      if (g) gset.add(g);
    }
    this._genes = Array.from(gset).sort();

    // Init vectors per gene
    const N = this._rows.length;
    for (const g of this._genes) this._exprByGene.set(g, new Float32Array(N));

    // Fill counts per (cell,gene)
    for (const r of rows) {
      const cid = String(r[cellKey] ?? "");
      const g = String(r[geneKey] ?? "");
      const idx = this._cellIndex.get(cid);
      if (idx == null || !g) continue;
      const vec = this._exprByGene.get(g);
      if (vec) vec[idx] += 1; // increment transcript count
    }

    console.log(
      `[XeniumAdapter] Ingested expression from transcripts.csv: genes=${this._genes.length}`
    );
  }

  async _ingestFromCellsWide(geneCols: string[]): Promise<void> {
    // Wide layout: cells.csv has many gene columns (binary or counts)
    // Use provided geneCols; if _genes empty, adopt them as list in this order.
    if (!geneCols.length) return;

    if (!this._genes.length) {
      this._genes = geneCols.slice();
    } else {
      // keep existing genes order; add any new ones at end
      const gset = new Set(this._genes);
      for (const g of geneCols) if (!gset.has(g)) this._genes.push(g);
    }

    const N = this._rows.length;
    for (const g of this._genes)
      if (!this._exprByGene.has(g))
        this._exprByGene.set(g, new Float32Array(N));

    for (let i = 0; i < N; i++) {
      const row = this._rows[i];
      for (const g of geneCols) {
        const v = toNumBool(row[g]);
        const vec = this._exprByGene.get(g);
        if (vec) vec[i] = Number.isFinite(v) ? v : 0;
      }
    }

    console.log(
      `[XeniumAdapter] Ingested expression from cells.csv wide columns: genes+=${geneCols.length}`
    );
  }

  async _parseCsvFile(file: File): Promise<ParsedTable> {
    const name = (file.name || "").toLowerCase();
    if (name.endsWith(".gz")) {
      const text = await fileToTextMaybeGz(file);
      if (!text || !text.trim()) return { headers: [], rows: [] };
      return parseCsvWithPapa(text);
    }
    return parseCsvWithPapa(file);
  }

  async _readTableOneOf(names: string[]): Promise<RowData[]> {
    const map = new Map(
      this.files.map((f) => [(f.webkitRelativePath || f.name).toLowerCase(), f])
    );
    let f = null;
    for (const n of names) {
      const needle = n.toLowerCase();
      f = Array.from(map.values()).find((file) => {
        const base = file.name.toLowerCase();
        const rel = (file.webkitRelativePath || "").toLowerCase();
        return base === needle || rel.endsWith("/" + needle);
      });
      if (f) break;
    }
    if (!f) return [];
    const { rows } = await this._parseCsvFile(f);
    return rows;
  }
}

/* ----------------- parsing & utility helpers ----------------- */

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
        // Skip rows that are entirely empty after normalization
        if (
          Object.keys(normalized).length === 0 ||
          Object.values(normalized).every((val) => val === "")
        )
          continue;
        rows.push(normalized);
      }
    };

    Papa.parse<RowData>(input as any, {
      header: true,
      skipEmptyLines: "greedy",
      worker: false,
      dynamicTyping: false,
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
        if (results) {
          headers = ensureHeaders(headers, results.meta?.fields);
          if (!rows.length && Array.isArray(results.data)) {
            pushRows(results.data as RowData[]);
          }
          if (results.errors?.length) {
            console.warn(
              "[XeniumAdapter] PapaParse completed with errors:",
              results.errors
            );
          }
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
  const lookup = buildRowKeyLookup(row);
  for (const key of keys) {
    if (!key) continue;
    const sourceKey = lookup.get(key);
    if (sourceKey && sourceKey !== key) {
      row[key] = row[sourceKey];
    }
    if (row[key] == null) row[key] = "";
  }
  return row;
}

function ensureHeaders(
  existing: string[],
  fields?: (string | undefined)[]
): string[] {
  if (existing.length || !fields?.length) return existing;
  const seen = new Set<string>();
  const trimmed = [];
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

function buildRowKeyLookup(row: RowData | null | undefined): Map<string, string> {
  const map = new Map<string, string>();
  if (!row) return map;
  for (const rawKey of Object.keys(row)) {
    if (!map.has(rawKey)) map.set(rawKey, rawKey);
    const trimmed = rawKey.trim();
    if (trimmed && !map.has(trimmed)) map.set(trimmed, rawKey);
  }
  return map;
}

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

function pickFirstPresent(objOrHeaders: any, keys: string[]): string | null {
  if (!objOrHeaders) return null;
  const keysList = Array.isArray(objOrHeaders)
    ? objOrHeaders
    : Object.keys(objOrHeaders);
  const lowerToOrig = new Map(
    keysList.map((k) => [String(k).toLowerCase(), k])
  );
  for (const k of keys) {
    const found = lowerToOrig.get(k.toLowerCase());
    if (found) return found;
  }
  return null;
}

function headersObj(headers: string[]): Record<string, boolean> {
  const o: Record<string, boolean> = {};
  for (const h of headers) o[h] = true;
  return o;
}

// Try centroid key variants; otherwise heuristic fallback
function detectCentroidKeys(sampleRow: RowData): { xKey: string | null; yKey: string | null } {
  const keys = Object.keys(sampleRow || {});
  const lower = keys.map((k) => k.toLowerCase());
  const xCandidates = [
    "cell_centroid_x",
    "x_centroid",
    "centroid_x",
    "cx",
    "x",
    "x_um",
    "x_px",
    "x_position",
    "xpos",
  ];
  const yCandidates = [
    "cell_centroid_y",
    "y_centroid",
    "centroid_y",
    "cy",
    "y",
    "y_um",
    "y_px",
    "y_position",
    "ypos",
  ];
  const xKey = firstMatch(keys, lower, xCandidates);
  const yKey = firstMatch(keys, lower, yCandidates);
  if (xKey && yKey) return { xKey, yKey };

  // Heuristic: numeric-looking x*/y* names
  const numish = keys.filter((k) => {
    const v = sampleRow[k];
    const n = Number(v);
    return typeof v === "string" || typeof v === "number"
      ? Number.isFinite(n)
      : false;
  });
  const xHeu = numish.find((k) =>
    /(^|_)x($|_|[a-z])|centroid.*x|x.*centroid|cx/i.test(k)
  );
  const yHeu = numish.find((k) =>
    /(^|_)y($|_|[a-z])|centroid.*y|y.*centroid|cy/i.test(k)
  );
  return { xKey: xHeu || null, yKey: yHeu || null };
}

function firstMatch(origKeys: string[], lowerKeys: string[], candidates: string[]): string | null {
  for (const cand of candidates) {
    const idx = lowerKeys.indexOf(cand.toLowerCase());
    if (idx !== -1) return origKeys[idx];
  }
  return null;
}

// Heuristic to find gene columns inside cells.csv (wide form)
function detectWideGeneColumns(sampleRow: RowData, knownGenes: string[] = []): { geneCols: string[] } {
  if (!sampleRow) return { geneCols: [] };
  const keys = Object.keys(sampleRow);
  const geneSet = new Set(knownGenes.map((g) => g.toLowerCase()));
  const geneLike = [];

  for (const k of keys) {
    if (k.startsWith("_")) continue; // skip internal
    const v = sampleRow[k];
    const numish = Number.isFinite(Number(v));
    // if column name matches a known gene OR the column is numeric and looks gene-ish (uppercase letters, etc.)
    const looksKnown = geneSet.has(k.toLowerCase());
    const looksGeneish =
      /^[A-Za-z0-9\-\._]+$/.test(k) && numish && k.length <= 30;
    if (looksKnown || looksGeneish) {
      // Avoid obvious non-gene fields:
      if (
        !/^(x|y|cell|id|barcode|area|volume|centroid|cx|cy|nuc|umit|umi|count|reads|cluster|leiden|louvain|class|type|annotation)/i.test(
          k
        )
      ) {
        geneLike.push(k);
      }
    }
  }
  return { geneCols: geneLike };
}
