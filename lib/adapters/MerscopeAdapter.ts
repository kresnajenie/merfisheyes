// src/data/adapters/MerscopeAdapter.ts
import { fileToTextMaybeGz } from "@/lib/utils/gzip";

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

    // --- read helpers
    const readTable = async (needles: string[]): Promise<RowData[]> => {
      const map = new Map(
        this.files.map((f) => [
          (f.webkitRelativePath || f.name).toLowerCase(),
          f,
        ])
      );
      let f = null;
      for (const n of needles) {
        const needle = n.toLowerCase();
        f = Array.from(map.values()).find((file) => {
          const base = file.name.toLowerCase();
          const rel = (file.webkitRelativePath || "").toLowerCase();
          return base === needle || rel.endsWith("/" + needle);
        });
        if (f) break;
      }
      if (!f) return [];
      const text = await fileToTextMaybeGz(f);
      if (!text || !text.trim()) return [];
      const delim = sniffDelimiter(text);
      const { headers, rows } = parseDelimited(text, delim);
      return rows;
    };

    // 1) cell_metadata.csv (positions live here)
    await this._onProgress?.(10, "Loading cell metadata...");
    const metaRows = await readTable(["cell_metadata.csv"]);
    if (!metaRows.length) {
      throw new Error("MERSCOPE: cell_metadata.csv not found or empty");
    }

    // 2) cell_categories.csv (cluster labels like "leiden")
    await this._onProgress?.(25, "Loading cluster categories...");
    const catRows = await readTable(["cell_categories.csv"]);
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
    const numRows = await readTable(["cell_numeric_categories.csv"]);
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
      if (catIndex) {
        out.leiden = catIndex.get(id) ?? "";
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
      const cbgRows = await readTable(["cell_by_gene.csv"]);
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
    this._clusterColumn =
      [
        "leiden",
        "cluster",
        "clusters",
        "cell_type",
        "celltype",
        "annotation",
        "class",
      ].find((k) => this._obsKeys.includes(k)) || null;

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

  getObsColumns(): string[] {
    return this._obsKeys;
  }

  getObsmEmbeddings(): string[] {
    const keys = ["X_spatial"];
    if (this._umap) keys.push("X_umap");
    return keys;
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

function sniffDelimiter(sampleText: string): string {
  const first = sampleText.split(/\r?\n/)[0] || "";
  const tabs = (first.match(/\t/g) || []).length;
  const commas = (first.match(/,/g) || []).length;
  const semis = (first.match(/;/g) || []).length;
  if (tabs > commas && tabs > semis) return "\t";
  if (semis > commas) return ";";
  return ","; // default
}

function parseDelimited(txt: string, delim: string): { headers: string[]; rows: RowData[] } {
  const lines = txt.split(/\r?\n/).filter(Boolean);
  if (!lines.length) return { headers: [], rows: [] };
  const headers = splitLine(lines[0], delim);
  const rows: RowData[] = [];
  for (let i = 1; i < lines.length; i++) {
    const cols = splitLine(lines[i], delim);
    if (!cols.length) continue;
    const obj: RowData = {};
    for (let j = 0; j < headers.length; j++) obj[headers[j]] = cols[j] ?? "";
    rows.push(obj);
  }
  return { headers, rows };
}

function splitLine(line: string, delim: string): string[] {
  const res: string[] = [];
  let cur = "";
  let inQ = false;
  for (let i = 0; i < line.length; i++) {
    const c = line[i];
    if (c === '"') {
      if (inQ && line[i + 1] === '"') {
        cur += '"';
        i++;
      } else inQ = !inQ;
    } else if (c === delim && !inQ) {
      res.push(cur);
      cur = "";
    } else cur += c;
  }
  res.push(cur);
  return res;
}
