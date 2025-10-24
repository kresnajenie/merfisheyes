// /src/data/adapters/XeniumAdapter.ts
import Papa, { ParseResult } from "papaparse";
import { ungzip } from "pako";

import { fileToTextMaybeGz } from "@/lib/utils/gzip";
import { shouldFilterGene } from "@/lib/utils/gene-filters";

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
  _featureTypes: Map<string, string>;

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
    this._featureTypes = new Map();
  }

  async initialize(
    files: File[],
    onProgress?: (progress: number, message: string) => Promise<void> | void,
  ): Promise<void> {
    const expandedArchives = await extractRelevantArchives(files);

    this.files = expandedArchives.length
      ? [...files, ...expandedArchives]
      : files;

    // 1) Load cells table (csv or csv.gz)
    await onProgress?.(10, "Loading cell metadata...");
    this._rows = await this._readTableOneOf(["cells.csv", "cells.csv.gz"]);
    if (!this._rows.length)
      throw new Error("Xenium: cells.csv(.gz) not found or empty");

    // 2) Save obs keys
    this._obsKeys = Object.keys(this._rows[0] || {});

    // Build index: choose a best cell id key and map id -> row idx
    await onProgress?.(20, "Building cell index...");
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
    await onProgress?.(30, "Loading cluster data...");
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
          '[XeniumAdapter] Cluster column looked per-cell unique; falling back to single "All" group.',
        );
        this._clusterColumn = null;
      }
    }
    if (!this._clusterColumn) {
      console.warn(
        'No celltype-like column found; will use single "All" group.',
      );
    }

    // 5) Genes: features, else transcripts
    //    (we'll still build expression from transcripts or cells-wide if possible)
    await onProgress?.(50, "Loading gene features...");
    this._featureTypes.clear();
    const featureFile = this._findFileOneOf([
      "features.tsv",
      "features.tsv.gz",
      "features.csv",
      "features.csv.gz",
      "cell_feature_matrix/features.tsv",
      "cell_feature_matrix/features.tsv.gz",
      "cell_feature_matrix/features.csv",
      "cell_feature_matrix/features.csv.gz",
    ]);

    if (featureFile) {
      try {
        const featureEntries = await parseFeaturesFile(featureFile);

        if (featureEntries.length) {
          const seen = new Set<string>();

          this._genes = [];
          for (const entry of featureEntries) {
            const gene = entry.gene;

            if (!gene || seen.has(gene)) continue;
            seen.add(gene);
            this._genes.push(gene);
            if (entry.type) this._featureTypes.set(gene, entry.type);
          }
        } else {
          console.warn(
            "[XeniumAdapter] features.tsv parsed but yielded no entries.",
          );
        }
      } catch (err) {
        console.warn("[XeniumAdapter] Failed parsing features.tsv:", err);
      }
    } else {
      console.warn(
        "[XeniumAdapter] features.tsv not found; falling back to cells.csv gene inference.",
      );
    }

    // 6) Expression: try transcripts.csv(.gz) (long format) first
    await onProgress?.(70, "Loading gene expression...");
    let gotExpr = false;

    // try {
    //   const tx = await this._readTableOneOf([
    //     "transcripts.csv",
    //     "transcripts.csv.gz",
    //   ]);

    //   if (tx.length) {
    //     await this._ingestFromTranscripts(tx, cellIdKey);
    //     gotExpr = this._genes.length > 0 && this._exprByGene.size > 0;
    //   }
    // } catch {}

    if (!gotExpr) {
      try {
        gotExpr = await this._ingestFromMatrixMarket(cellIdKey);
      } catch (e) {
        console.warn("[XeniumAdapter] matrix.mtx ingestion failed:", e);
      }
    }

    // 7) Fallback: disabled to avoid huge allocations from cells.csv wide columns

    this._applyGeneFilters();

    if (!this._genes.length)
      console.warn("No gene names found (or none usable for expression).");

    await onProgress?.(90, "Finalizing dataset...");
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
        "[XeniumAdapter] Could not detect centroid columns; returning 0 points",
      );

      return { coordinates: [], dimensions: 2 };
    }

    const numRows = rows.length;
    const coords = [];
    let validCount = 0;

    for (const r of rows) {
      const x = toNum(r[xKey]);
      const y = toNum(r[yKey]);

      if (Number.isFinite(x) && Number.isFinite(y)) {
        coords.push([x, y]);
        validCount++;
      }
    }

    // Validate: at least 90% of coordinates must be valid
    const validPercent = (validCount / numRows) * 100;

    if (validPercent < 90) {
      throw new Error(
        `[XeniumAdapter] Invalid coordinates: only ${validPercent.toFixed(1)}% are valid (need ≥90%). ` +
          `Found ${validCount} valid out of ${numRows} rows using keys: ${xKey}, ${yKey}`,
      );
    }

    const droppedCount = numRows - validCount;

    if (droppedCount > 0) {
      console.warn(
        `[XeniumAdapter] Dropped ${droppedCount} rows with invalid centroids (from ${numRows})`,
      );
    }

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

  private _findFileOneOf(names: string[]): File | null {
    const map = new Map(
      this.files.map((f) => [
        (f.webkitRelativePath || f.name).toLowerCase(),
        f,
      ]),
    );

    for (const n of names) {
      const needle = n.toLowerCase();
      const match = Array.from(map.values()).find((file) => {
        const base = file.name.toLowerCase();
        const rel = (file.webkitRelativePath || "").toLowerCase();

        return base === needle || base.endsWith(needle) || rel.endsWith(needle);
      });

      if (match) return match;
    }

    return null;
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

    const preferred = candidates.filter((f) => {
      const p = (f.webkitRelativePath || f.name || "").toLowerCase();

      return p.includes("analysis/clustering/gene_expression_graphclust/");
    });
    const prioritized = preferred.length ? preferred : candidates;

    const nCells = this._rows.length || 1;
    const maxReasonableUnique = Math.max(50, Math.floor(0.1 * nCells));
    let best: {
      file: File;
      labelKey: string;
      lut: Map<string, string>;
      assigned: number;
      uniqCount: number;
    } | null = null;

    for (const f of prioritized) {
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
        const uniqVals = new Set<string>();

        for (const r of rows) {
          const cid = String(r[cellIdKeyAlt] ?? "").trim();
          const lab = String(r[labelKey] ?? "").trim();

          if (cid) {
            lut.set(cid, lab);
            if (lab) uniqVals.add(lab);
          }
        }
        if (!lut.size) continue;

        let assigned = 0;
        const matchedVals: string[] = [];

        for (const r of this._rows) {
          const cid = String(r[cellIdKey] ?? "").trim();

          if (!cid) continue;
          const lab = lut.get(cid);

          if (lab != null) {
            assigned++;
            if (lab) matchedVals.push(lab);
          }
        }

        if (assigned === 0) continue;

        const uniqCount = new Set(matchedVals).size;

        if (uniqCount === 0) continue;
        if (uniqCount > maxReasonableUnique) continue;

        if (
          !best ||
          uniqCount > best.uniqCount ||
          (uniqCount === best.uniqCount && assigned > best.assigned)
        ) {
          best = { file: f, labelKey, lut, assigned, uniqCount };
        }
      } catch (e) {
        console.warn(
          "[XeniumAdapter] Failed reading clustering file:",
          f.name,
          e,
        );
      }
    }

    if (best) {
      const { labelKey, lut, file, assigned, uniqCount } = best;

      for (const r of this._rows) {
        const cid = String(r[cellIdKey] ?? "").trim();

        if (!cid) continue;
        const lab = lut.get(cid);

        if (lab != null) r[labelKey] = lab;
      }
      this._clusterColumn = labelKey;
      if (!this._obsKeys.includes(labelKey)) this._obsKeys.push(labelKey);
      console.log(
        `[XeniumAdapter] Joined ${assigned} labels from ${file.name} onto cells (column "${labelKey}", unique=${uniqCount})`,
      );
    }
  }

  async _ingestFromTranscripts(
    rows: RowData[],
    cellIdKey: string,
  ): Promise<void> {
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
        "[XeniumAdapter] transcripts.csv lacks recognizable cell/gene columns; skipping.",
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
      `[XeniumAdapter] Ingested expression from transcripts.csv: genes=${this._genes.length}`,
    );
  }

  async _ingestFromMatrixMarket(cellIdKey: string): Promise<boolean> {
    const matrixFile = this._findFileOneOf([
      "matrix.mtx",
      "matrix.mtx.gz",
      "cell_feature_matrix/matrix.mtx",
      "cell_feature_matrix/matrix.mtx.gz",
    ]);

    if (!matrixFile) return false;

    const nCells = this._rows.length;

    if (!nCells) return false;

    const barcodesFile = this._findFileOneOf([
      "barcodes.tsv",
      "barcodes.tsv.gz",
      "cell_feature_matrix/barcodes.tsv",
      "cell_feature_matrix/barcodes.tsv.gz",
    ]);
    let barcodes: string[] | null = null;

    if (barcodesFile) {
      try {
        barcodes = await readBarcodesFile(barcodesFile);
      } catch (e) {
        console.warn("[XeniumAdapter] Failed reading barcodes.tsv:", e);
      }
    }

    if (!this._genes.length) {
      console.warn(
        "[XeniumAdapter] matrix.mtx found but gene list empty; cannot align columns.",
      );

      return false;
    }

    this._exprByGene.clear();

    let headerParsed = false;
    let rowCount = 0;
    let colCount = 0;
    let filled = 0;
    let transpose = false;
    let cellsDim = nCells;
    const genesDim = this._genes.length;

    for await (const line of readLinesFromFile(matrixFile)) {
      if (!line || line.startsWith("%")) continue;
      const trimmed = line.trim();

      if (!trimmed) continue;
      const parts = trimmed.split(/\s+/);

      if (!headerParsed) {
        if (parts.length >= 3) {
          rowCount = Number(parts[0]);
          colCount = Number(parts[1]);
          headerParsed = true;
          if (!Number.isFinite(rowCount) || !Number.isFinite(colCount)) {
            console.warn(
              "[XeniumAdapter] matrix.mtx header invalid; aborting ingest.",
            );

            return false;
          }
          const rowsMatchCells = rowCount === nCells;
          const colsMatchGenes = colCount === genesDim;
          const rowsMatchGenes = rowCount === genesDim;
          const colsMatchCells = colCount === nCells;

          if (!(rowsMatchCells && colsMatchGenes)) {
            if (colsMatchCells && rowsMatchGenes) {
              transpose = true;
            } else if (!rowsMatchCells && colsMatchCells) {
              transpose = true;
              console.warn(
                `[XeniumAdapter] matrix.mtx row count (${rowCount}) did not match cell count (${nCells}); treating matrix as transposed.`,
              );
            } else if (!colsMatchGenes && rowsMatchGenes) {
              transpose = true;
              console.warn(
                `[XeniumAdapter] matrix.mtx column count (${colCount}) aligned with cells and rows with genes; transposing interpretation.`,
              );
            } else {
              if (!rowsMatchCells)
                console.warn(
                  `[XeniumAdapter] matrix.mtx row count (${rowCount}) != cell count (${nCells}); using positional fallback.`,
                );
              if (!colsMatchGenes)
                console.warn(
                  `[XeniumAdapter] matrix.mtx column count (${colCount}) != feature count (${genesDim}); extra indices will be ignored.`,
                );
            }
          }
          cellsDim = transpose ? colCount : rowCount;
          if (barcodes && cellsDim !== barcodes.length) {
            console.warn(
              `[XeniumAdapter] matrix.mtx ${transpose ? "column" : "row"} count (${cellsDim}) != barcodes length (${barcodes.length}); using positional fallback.`,
            );
          }
        }
        continue;
      }
      if (parts.length < 3) continue;
      const rowIdx = Number(parts[0]) - 1;
      const colIdx = Number(parts[1]) - 1;
      const value = Number(parts[2]);

      if (
        !Number.isFinite(rowIdx) ||
        !Number.isFinite(colIdx) ||
        !Number.isFinite(value)
      )
        continue;

      const cellIndexRaw = transpose ? colIdx : rowIdx;
      const geneIndexRaw = transpose ? rowIdx : colIdx;

      let cellIdx: number | undefined;

      if (cellIndexRaw >= 0) {
        if (barcodes && cellIndexRaw < barcodes.length) {
          const cid = barcodes[cellIndexRaw]?.trim();

          if (cid) {
            cellIdx = this._cellIndex.get(cid);
            if (cellIdx == null && cid.endsWith("-1"))
              cellIdx = this._cellIndex.get(cid.slice(0, -2));
          }
        }
        if (cellIdx == null && cellIndexRaw < nCells) {
          cellIdx = cellIndexRaw;
        }
      }
      if (cellIdx == null) continue;

      if (geneIndexRaw < 0 || geneIndexRaw >= this._genes.length) continue;
      const gene = this._genes[geneIndexRaw];

      if (!gene) continue;

      let vec = this._exprByGene.get(gene);

      if (!vec) {
        try {
          vec = new Float32Array(nCells);
        } catch (err) {
          console.warn(
            "[XeniumAdapter] Unable to allocate expression vector from matrix.mtx; aborting ingest.",
            err,
          );
          this._exprByGene.clear();

          return false;
        }
        this._exprByGene.set(gene, vec);
      }
      vec[cellIdx] += value;
      filled++;
    }

    if (!filled) {
      console.warn(
        "[XeniumAdapter] matrix.mtx parsed but produced no gene entries.",
      );
      this._exprByGene.clear();

      return false;
    }

    console.log(
      `[XeniumAdapter] Ingested expression from matrix.mtx: genes=${this._genes.length}, nonzeros=${filled}`,
    );

    return true;
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
      `[XeniumAdapter] Ingested expression from cells.csv wide columns: genes+=${geneCols.length}`,
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
    const f = this._findFileOneOf(names);

    if (!f) return [];
    const { rows } = await this._parseCsvFile(f);

    return rows;
  }

  private _applyGeneFilters(): void {
    if (!this._genes.length) return;
    const filtered: string[] = [];
    let removed = 0;
    const seen = new Set<string>();

    for (const gene of this._genes) {
      if (!gene || seen.has(gene)) continue;
      seen.add(gene);
      const type = this._featureTypes.get(gene) || "";

      if (shouldFilterGene(gene, type)) {
        this._exprByGene.delete(gene);
        removed++;
      } else {
        filtered.push(gene);
      }
    }
    if (removed) {
      console.log(
        `[XeniumAdapter] Filtered ${removed} control/unassigned features from gene list.`,
      );
    }
    this._genes = filtered;
  }
}

/* ----------------- parsing & utility helpers ----------------- */

type ParsedTable = { headers: string[]; rows: RowData[] };
interface FeatureEntry {
  gene: string;
  type?: string;
}

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
      transformHeader: (header: string | undefined) => (header ?? "").trim(),
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
              results.errors,
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
  fields?: (string | undefined)[],
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

async function readBarcodesFile(file: File): Promise<string[]> {
  const barcodes: string[] = [];

  for await (const line of readLinesFromFile(file)) {
    if (!line) continue;
    if (line.startsWith("#")) continue;
    const trimmed = line.trim();

    if (!trimmed) continue;
    const first = trimmed.split(/\t|,/)[0]?.trim();

    if (first) barcodes.push(first);
  }

  return barcodes;
}

async function* readLinesFromFile(file: File): AsyncGenerator<string> {
  const name = (file.name || "").toLowerCase();
  const isGz = name.endsWith(".gz");

  if (isGz && typeof DecompressionStream === "undefined") {
    const text = await fileToTextMaybeGz(file);
    const lines = text.split(/\r?\n/);

    for (const line of lines) {
      yield line;
    }

    return;
  }

  let stream: ReadableStream<Uint8Array> = file.stream();

  if (isGz && typeof DecompressionStream !== "undefined") {
    stream = stream.pipeThrough(new DecompressionStream("gzip"));
  }

  const reader = stream.getReader();
  const decoder = new TextDecoder();
  let buffer = "";

  try {
    while (true) {
      const { value, done } = await reader.read();

      if (done) break;
      buffer += decoder.decode(value, { stream: true });
      let newlineIndex = buffer.indexOf("\n");

      while (newlineIndex !== -1) {
        let line = buffer.slice(0, newlineIndex);

        if (line.endsWith("\r")) line = line.slice(0, -1);
        yield line;
        buffer = buffer.slice(newlineIndex + 1);
        newlineIndex = buffer.indexOf("\n");
      }
    }
  } finally {
    reader.releaseLock();
  }

  buffer += decoder.decode();
  if (buffer) {
    if (buffer.endsWith("\r")) buffer = buffer.slice(0, -1);
    yield buffer;
  }
}

async function parseFeaturesFile(file: File): Promise<FeatureEntry[]> {
  const lines: string[] = [];

  for await (const rawLine of readLinesFromFile(file)) {
    if (rawLine == null) continue;
    const trimmed = rawLine.trim();

    if (!trimmed || trimmed.startsWith("#")) continue;
    lines.push(trimmed);
  }
  if (!lines.length) return [];

  const delim = lines[0].includes("\t") ? "\t" : ",";
  const firstParts = lines[0].split(delim).map((s) => s.trim());
  const looksLikeHeader = firstParts.some((p) =>
    /gene|feature|target|type|id/.test(p.toLowerCase()),
  );

  let startIdx = 0;
  let geneIdx = -1;
  let typeIdx = -1;

  if (looksLikeHeader) {
    startIdx = 1;
    firstParts.forEach((col, idx) => {
      const lower = col.toLowerCase();

      if (
        geneIdx === -1 &&
        [
          "gene",
          "gene_symbol",
          "gene_id",
          "gene_name",
          "feature_name",
          "feature",
          "target",
          "name",
          "id",
        ].includes(lower)
      ) {
        geneIdx = idx;
      } else if (
        typeIdx === -1 &&
        (/type/.test(lower) || /class/.test(lower) || /category/.test(lower))
      ) {
        typeIdx = idx;
      }
    });
    if (geneIdx === -1) geneIdx = Math.min(1, firstParts.length - 1);
  } else {
    startIdx = 0;
    geneIdx = firstParts.length >= 2 ? 1 : 0;
    typeIdx = firstParts.length >= 3 ? 2 : -1;
  }

  const entries: FeatureEntry[] = [];

  for (let i = startIdx; i < lines.length; i++) {
    const parts = lines[i].split(delim);

    if (!parts.length) continue;
    const gene = (parts[geneIdx] ?? "").trim();

    if (!gene) continue;
    const type = typeIdx >= 0 ? (parts[typeIdx] ?? "").trim() : "";

    entries.push({ gene, type });
  }

  return entries;
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
    keysList.map((k) => [String(k).toLowerCase(), k]),
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
function detectCentroidKeys(sampleRow: RowData): {
  xKey: string | null;
  yKey: string | null;
} {
  const keys = Object.keys(sampleRow || {});
  const lower = keys.map((k) => k.toLowerCase());
  const xCandidates = [
    "cell_centroid_x",
    "x_centroid",
    "centroid_x",
    "center_x",
    "centerx",
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
    "center_y",
    "centery",
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
    /(^|_)x($|_|[a-z])|centroid.*x|x.*centroid|cx/i.test(k),
  );
  const yHeu = numish.find((k) =>
    /(^|_)y($|_|[a-z])|centroid.*y|y.*centroid|cy/i.test(k),
  );

  return { xKey: xHeu || null, yKey: yHeu || null };
}

function firstMatch(
  origKeys: string[],
  lowerKeys: string[],
  candidates: string[],
): string | null {
  for (const cand of candidates) {
    const idx = lowerKeys.indexOf(cand.toLowerCase());

    if (idx !== -1) return origKeys[idx];
  }

  return null;
}

// Heuristic to find gene columns inside cells.csv (wide form)
function detectWideGeneColumns(
  sampleRow: RowData,
  knownGenes: string[] = [],
): { geneCols: string[] } {
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
          k,
        )
      ) {
        geneLike.push(k);
      }
    }
  }

  return { geneCols: geneLike };
}

type TarEntry = { name: string; data: Uint8Array };

async function extractRelevantArchives(files: File[]): Promise<File[]> {
  const extracted: File[] = [];

  for (const file of files) {
    const fileName = (file.name || "").toLowerCase();

    if (
      !fileName.endsWith(".tar") &&
      !fileName.endsWith(".tar.gz") &&
      !fileName.endsWith(".tgz")
    )
      continue;
    if (!/analysis|cell_feature_matrix/.test(fileName)) continue;
    try {
      const entries = await readTarEntriesFromFile(file);

      for (const entry of entries) {
        const pathLower = entry.name.toLowerCase();
        const isClusterCsv =
          /analysis\/clustering\//.test(pathLower) &&
          /\.(csv|tsv)$/.test(pathLower);
        const isFeatureTsv =
          /cell_feature_matrix\//.test(pathLower) &&
          /features\.tsv(\.gz)?$/.test(pathLower);
        const isBarcodesTsv =
          /cell_feature_matrix\//.test(pathLower) &&
          /barcodes\.tsv(\.gz)?$/.test(pathLower);
        const isMatrixMtx =
          /cell_feature_matrix\//.test(pathLower) &&
          /matrix\.mtx(\.gz)?$/.test(pathLower);

        if (!isClusterCsv && !isFeatureTsv && !isBarcodesTsv && !isMatrixMtx)
          continue;

        const blob = new File([entry.data], entry.name, {
          type: "text/plain",
          lastModified: file.lastModified,
        });

        try {
          Object.defineProperty(blob, "webkitRelativePath", {
            value: entry.name,
            configurable: true,
          });
        } catch {}
        extracted.push(blob);
      }
    } catch (err) {
      console.warn("[XeniumAdapter] Failed to extract archive", file.name, err);
    }
  }

  return extracted;
}

async function readTarEntriesFromFile(file: File): Promise<TarEntry[]> {
  const lower = (file.name || "").toLowerCase();
  let bytes: Uint8Array;

  if (lower.endsWith(".tar.gz") || lower.endsWith(".tgz")) {
    bytes = await decompressGzipToUint8Array(file);
  } else if (lower.endsWith(".tar")) {
    const buf = await file.arrayBuffer();

    bytes = new Uint8Array(buf);
  } else {
    return [];
  }

  return parseTarArchive(bytes);
}

async function decompressGzipToUint8Array(file: File): Promise<Uint8Array> {
  if (typeof DecompressionStream !== "undefined") {
    const ds = new DecompressionStream("gzip");
    const decompressed = file.stream().pipeThrough(ds);
    const arrayBuffer = await new Response(decompressed).arrayBuffer();

    return new Uint8Array(arrayBuffer);
  }
  const buffer = await file.arrayBuffer();

  return ungzip(new Uint8Array(buffer));
}

function parseTarArchive(bytes: Uint8Array): TarEntry[] {
  const entries: TarEntry[] = [];
  const blockSize = 512;
  const decoder = new TextDecoder("utf-8");
  let offset = 0;

  while (offset + blockSize <= bytes.length) {
    const header = bytes.subarray(offset, offset + blockSize);

    offset += blockSize;

    if (isZeroBlock(header)) {
      // end of archive marker (two consecutive zero blocks)
      if (offset + blockSize > bytes.length) break;
      const next = bytes.subarray(offset, offset + blockSize);

      if (isZeroBlock(next)) break;
      else continue;
    }

    const nameRaw = decoder
      .decode(header.subarray(0, 100))
      .replace(/\0.*$/, "");
    const prefixRaw = decoder
      .decode(header.subarray(345, 500))
      .replace(/\0.*$/, "");
    const sizeRaw = decoder
      .decode(header.subarray(124, 136))
      .replace(/\0.*$/, "")
      .trim();

    const fullName = (prefixRaw ? prefixRaw + "/" : "") + nameRaw;

    if (!fullName) {
      const sizeSkip = Math.ceil(0 / blockSize) * blockSize;

      offset += sizeSkip;
      continue;
    }

    const size = sizeRaw ? parseInt(sizeRaw, 8) || 0 : 0;
    const typeflag = header[156];

    const dataStart = offset;
    const dataEnd = dataStart + size;
    const fileData = bytes.subarray(dataStart, dataEnd);
    const paddedSize = Math.ceil(size / blockSize) * blockSize;

    offset += paddedSize;

    if (typeflag === 0 || typeflag === 48 /* '0' */) {
      entries.push({ name: fullName, data: fileData.slice() });
    }
  }

  return entries;
}

function isZeroBlock(block: Uint8Array): boolean {
  for (let i = 0; i < block.length; i++) {
    if (block[i] !== 0) return false;
  }

  return true;
}
