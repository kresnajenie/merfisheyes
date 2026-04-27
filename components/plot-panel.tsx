"use client";

import { useEffect, useMemo, useState } from "react";
import { Rnd } from "react-rnd";

import type { StandardizedDataset } from "@/lib/StandardizedDataset";
import {
  usePanelDatasetStore,
  usePanelVisualizationStore,
} from "@/lib/hooks/usePanelStores";
import { getEffectiveColumnType } from "@/lib/utils/column-type-utils";

import { CelltypeBarplot } from "./plots/celltype-barplot";
import { GeneCelltypeBoxplot } from "./plots/gene-celltype-boxplot";
import { GeneHistogram } from "./plots/gene-histogram";
import { NumericalHistogram } from "./plots/numerical-histogram";
import { SecondaryGroupControls } from "./secondary-group-controls";

// Default size as fractions of the viewport (so the panel scales with the
// window). 0.375 × 0.45 ≈ 480×360 at 1280×800.
const DEFAULT_W_FRAC = 0.375;
const DEFAULT_H_FRAC = 0.45;
const MIN_W = 320;
const MIN_H = 200;
// Left margin clears the sidebar buttons (left-4 + w-14 = 16+56 = 72px) plus a small gap.
const MARGIN_LEFT = 88;
// Extra clearance below so the spatial scale bar (bottom-6 left-6 in the
// scene canvas) stays visible under the panel.
const MARGIN_BOTTOM = 72;

const DEFAULT_TOP_N = 10;
const MAX_TOP_N = 50;
const HEADER_H = 36;

// Parse "20k" → 20000, "1.5m" → 1_500_000, plain numbers as-is. Returns null
// for empty / unparseable input.
function parseNumericShorthand(s: string): number | null {
  const t = s.trim().toLowerCase();
  if (!t) return null;
  const m = t.match(/^([\d,.]+)\s*([kmb])?$/);
  if (!m) return null;
  const n = parseFloat(m[1].replace(/,/g, ""));
  if (!Number.isFinite(n) || n <= 0) return null;
  const mult =
    m[2] === "k" ? 1_000 : m[2] === "m" ? 1_000_000 : m[2] === "b" ? 1_000_000_000 : 1;
  return n * mult;
}

function quantile(sorted: Float32Array, q: number): number {
  const n = sorted.length;
  if (n === 0) return 0;
  const pos = (n - 1) * q;
  const lo = Math.floor(pos);
  const hi = Math.ceil(pos);
  if (lo === hi) return sorted[lo];
  return sorted[lo] + (pos - lo) * (sorted[hi] - sorted[lo]);
}

function csvEscape(v: string | number): string {
  const s = String(v);
  if (/[,"\n]/.test(s)) return `"${s.replace(/"/g, '""')}"`;
  return s;
}

function readClusterPerCell(
  cluster: { valueIndices?: any; uniqueValues?: any; values: any[] },
): (i: number) => string {
  if (cluster.valueIndices && cluster.uniqueValues) {
    return (i) => cluster.uniqueValues[cluster.valueIndices[i]];
  }
  return (i) => String(cluster.values[i]);
}

async function buildBoxplotCsv(
  dataset: StandardizedDataset,
  primary: string,
  gene: string,
  secondary: string | null,
  selectedSecondaryValues: Set<string> | null,
  selectedCelltypes: Set<string> | null,
  topN: number,
): Promise<string> {
  const expr = await dataset.getGeneExpression(gene);
  if (!expr) return "";
  const primaryCluster = dataset.clusters?.find((c) => c.column === primary);
  if (!primaryCluster) return "";
  const cellCount = primaryCluster.valueIndices
    ? primaryCluster.valueIndices.length
    : primaryCluster.values.length;
  const primaryAt = readClusterPerCell(primaryCluster);

  // Bucket gene values per (primary[, secondary]).
  const buckets = new Map<string, Float32Array>();

  if (secondary && selectedSecondaryValues && selectedCelltypes) {
    const secondaryCluster = dataset.clusters?.find(
      (c) => c.column === secondary,
    );
    if (!secondaryCluster) return "";
    const secondaryAt = readClusterPerCell(secondaryCluster);
    const sizes = new Map<string, number>();
    for (let i = 0; i < cellCount; i++) {
      const p = primaryAt(i);
      if (!selectedCelltypes.has(p)) continue;
      const s = secondaryAt(i);
      if (!selectedSecondaryValues.has(s)) continue;
      const k = `${p}\t${s}`;
      sizes.set(k, (sizes.get(k) ?? 0) + 1);
    }
    const cursors = new Map<string, number>();
    for (const [k, n] of sizes) {
      buckets.set(k, new Float32Array(n));
      cursors.set(k, 0);
    }
    for (let i = 0; i < cellCount; i++) {
      const p = primaryAt(i);
      if (!selectedCelltypes.has(p)) continue;
      const s = secondaryAt(i);
      if (!selectedSecondaryValues.has(s)) continue;
      const k = `${p}\t${s}`;
      const arr = buckets.get(k);
      if (!arr) continue;
      const c = cursors.get(k)!;
      arr[c] = expr[i] ?? 0;
      cursors.set(k, c + 1);
    }
  } else {
    const sizes = new Map<string, number>();
    for (let i = 0; i < cellCount; i++) {
      const p = primaryAt(i);
      sizes.set(p, (sizes.get(p) ?? 0) + 1);
    }
    const cursors = new Map<string, number>();
    for (const [k, n] of sizes) {
      buckets.set(k, new Float32Array(n));
      cursors.set(k, 0);
    }
    for (let i = 0; i < cellCount; i++) {
      const p = primaryAt(i);
      const arr = buckets.get(p);
      if (!arr) continue;
      const c = cursors.get(p)!;
      arr[c] = expr[i] ?? 0;
      cursors.set(p, c + 1);
    }
  }

  type Row = {
    primary: string;
    secondary?: string;
    n: number;
    q1: number;
    median: number;
    q3: number;
    max: number;
  };
  const rows: Row[] = [];
  for (const [k, arr] of buckets) {
    const sorted = arr.slice().sort();
    const n = sorted.length;
    const row: Row = {
      primary: secondary ? k.split("\t")[0] : k,
      n,
      q1: quantile(sorted, 0.25),
      median: quantile(sorted, 0.5),
      q3: quantile(sorted, 0.75),
      max: n > 0 ? sorted[n - 1] : 0,
    };
    if (secondary) row.secondary = k.split("\t")[1];
    rows.push(row);
  }

  if (!secondary) {
    rows.sort((a, b) => b.median - a.median);
    const cap = Math.max(0, topN);
    if (rows.length > cap) rows.length = cap;
  }

  const headers = secondary
    ? [primary, secondary, "n", "q1", "median", "q3", "max"]
    : [primary, "n", "q1", "median", "q3", "max"];
  const lines = [headers.map(csvEscape).join(",")];
  for (const r of rows) {
    const cells = secondary
      ? [r.primary, r.secondary ?? "", r.n, r.q1, r.median, r.q3, r.max]
      : [r.primary, r.n, r.q1, r.median, r.q3, r.max];
    lines.push(cells.map(csvEscape).join(","));
  }
  return lines.join("\n") + "\n";
}

function buildBarplotCsv(
  dataset: StandardizedDataset,
  primary: string,
  secondary: string | null,
  selectedSecondaryValues: Set<string> | null,
  selectedCelltypes: Set<string> | null,
  topN: number,
): string {
  const primaryCluster = dataset.clusters?.find((c) => c.column === primary);
  if (!primaryCluster) return "";
  const cellCount = primaryCluster.valueIndices
    ? primaryCluster.valueIndices.length
    : primaryCluster.values.length;
  const primaryAt = readClusterPerCell(primaryCluster);

  if (secondary && selectedSecondaryValues && selectedCelltypes) {
    const secondaryCluster = dataset.clusters?.find(
      (c) => c.column === secondary,
    );
    if (!secondaryCluster) return "";
    const secondaryAt = readClusterPerCell(secondaryCluster);
    const counts = new Map<string, Map<string, number>>();
    for (let i = 0; i < cellCount; i++) {
      const p = primaryAt(i);
      if (!selectedCelltypes.has(p)) continue;
      const s = secondaryAt(i);
      if (!selectedSecondaryValues.has(s)) continue;
      let inner = counts.get(p);
      if (!inner) {
        inner = new Map();
        counts.set(p, inner);
      }
      inner.set(s, (inner.get(s) ?? 0) + 1);
    }
    const lines = [
      [primary, secondary, "cells"].map(csvEscape).join(","),
    ];
    for (const [p, inner] of counts) {
      for (const [s, n] of inner) {
        lines.push([p, s, n].map(csvEscape).join(","));
      }
    }
    return lines.join("\n") + "\n";
  }

  const counts = new Map<string, number>();
  for (let i = 0; i < cellCount; i++) {
    const p = primaryAt(i);
    counts.set(p, (counts.get(p) ?? 0) + 1);
  }
  const rows = Array.from(counts.entries()).sort((a, b) => b[1] - a[1]);
  const cap = Math.max(0, topN);
  const trimmed = rows.length > cap ? rows.slice(0, cap) : rows;
  const lines = [[primary, "cells"].map(csvEscape).join(",")];
  for (const [p, n] of trimmed) {
    lines.push([p, n].map(csvEscape).join(","));
  }
  return lines.join("\n") + "\n";
}

export function PlotPanel() {
  const setPlotPanelOpen = usePanelVisualizationStore(
    (s) => s.setPlotPanelOpen,
  );
  const selectedColumn = usePanelVisualizationStore((s) => s.selectedColumn);
  const selectedCelltypes = usePanelVisualizationStore(
    (s) => s.selectedCelltypes,
  );
  const toggleCelltype = usePanelVisualizationStore((s) => s.toggleCelltype);
  const columnTypeOverrides = usePanelVisualizationStore(
    (s) => s.columnTypeOverrides,
  );
  const clusterVersion = usePanelVisualizationStore((s) => s.clusterVersion);
  const incrementClusterVersion = usePanelVisualizationStore(
    (s) => s.incrementClusterVersion,
  );
  const colormap = usePanelVisualizationStore((s) => s.colormap);
  const selectedGene = usePanelVisualizationStore((s) => s.selectedGene);
  const secondaryColumn = usePanelVisualizationStore((s) => s.secondaryColumn);
  const selectedSecondaryValues = usePanelVisualizationStore(
    (s) => s.selectedSecondaryValues,
  );
  const secondaryPaletteOverrides = usePanelVisualizationStore(
    (s) => s.secondaryPaletteOverrides,
  );
  const setSecondaryColumn = usePanelVisualizationStore(
    (s) => s.setSecondaryColumn,
  );
  const toggleSecondaryValue = usePanelVisualizationStore(
    (s) => s.toggleSecondaryValue,
  );
  const setSecondaryPaletteOverride = usePanelVisualizationStore(
    (s) => s.setSecondaryPaletteOverride,
  );

  const dataset = usePanelDatasetStore((s) => {
    const id = s.currentDatasetId;
    const ds = id ? s.datasets.get(id) : null;
    return ds && "spatial" in ds ? (ds as StandardizedDataset) : null;
  });

  // Size stored as a fraction of the viewport so the panel scales with the
  // window. User manual resize updates these fractions.
  const [sizeFrac, setSizeFrac] = useState({
    width: DEFAULT_W_FRAC,
    height: DEFAULT_H_FRAC,
  });
  // Offsets from bottom-left corner of viewport (panel anchors here).
  const [offsets, setOffsets] = useState({
    left: MARGIN_LEFT,
    bottom: MARGIN_BOTTOM,
  });
  const [minimized, setMinimized] = useState(false);
  const [topN, setTopN] = useState(DEFAULT_TOP_N);
  // When a gene is selected with a categorical column, show the boxplot by
  // default and let the user toggle to a histogram of the gene values.
  const [plotView, setPlotView] = useState<"box" | "histogram">("box");
  // User-controllable y-axis max. Empty string = auto. Accepts shorthand
  // like "20k" → 20000 and "1.5m" → 1500000.
  const [yMaxInput, setYMaxInput] = useState("");
  const yMax = useMemo(() => parseNumericShorthand(yMaxInput), [yMaxInput]);
  const [showOutliers, setShowOutliers] = useState(false);
  const [density, setDensity] = useState(false);
  // Bump on window resize so position/size recompute from fractions.
  // Use ResizeObserver on documentElement (more reliable than window.resize
  // when the page is in a frame, devtools is open, or the OS reports
  // resizes via different events) plus a window.resize fallback.
  const [, setResizeTick] = useState(0);

  useEffect(() => {
    const bump = () => setResizeTick((t) => t + 1);
    window.addEventListener("resize", bump);
    let ro: ResizeObserver | null = null;
    if (typeof ResizeObserver !== "undefined") {
      ro = new ResizeObserver(bump);
      ro.observe(document.documentElement);
    }
    return () => {
      window.removeEventListener("resize", bump);
      if (ro) ro.disconnect();
    };
  }, []);

  // Compute pixel size from fractions × current window size, then clamp to
  // a max that keeps the panel inside the viewport with a small top margin.
  const TOP_MARGIN = 16;
  const viewportH = typeof window === "undefined" ? 0 : window.innerHeight;
  const viewportW = typeof window === "undefined" ? 0 : window.innerWidth;
  const desiredH = sizeFrac.height * viewportH;
  const desiredW = sizeFrac.width * viewportW;
  const maxHeight = Math.max(HEADER_H, viewportH - offsets.bottom - TOP_MARGIN);
  const maxWidth = Math.max(MIN_W, viewportW - offsets.left - TOP_MARGIN);
  const effectiveHeight = minimized
    ? HEADER_H
    : Math.max(MIN_H, Math.min(desiredH, maxHeight));
  const effectiveWidth = Math.max(MIN_W, Math.min(desiredW, maxWidth));
  const position = {
    x: Math.max(0, offsets.left),
    y: Math.max(TOP_MARGIN, viewportH - effectiveHeight - offsets.bottom),
  };

  const columnType =
    dataset && selectedColumn
      ? getEffectiveColumnType(selectedColumn, dataset, columnTypeOverrides)
      : null;
  const isCategorical = columnType === "categorical";
  const isNumerical = columnType === "numerical";
  const hasGene = !!selectedGene;

  // What plot is being shown right now
  // - gene + categorical: boxplot (default) or gene histogram (toggle)
  // - gene + numerical or gene without column: gene histogram only
  // - no gene + categorical: cell count barplot
  // - no gene + numerical: column histogram
  let activePlot:
    | "gene-box"
    | "gene-histogram"
    | "celltype-barplot"
    | "numerical-histogram"
    | "empty" = "empty";
  if (dataset) {
    if (hasGene && isCategorical) {
      activePlot = plotView === "box" ? "gene-box" : "gene-histogram";
    } else if (hasGene) {
      activePlot = "gene-histogram";
    } else if (isCategorical) {
      activePlot = "celltype-barplot";
    } else if (isNumerical) {
      activePlot = "numerical-histogram";
    }
  }
  const showViewToggle = hasGene && isCategorical && !minimized;
  const showTopN =
    !minimized &&
    (activePlot === "celltype-barplot" || activePlot === "gene-box");
  // Secondary grouping controls show when the active plot is something a
  // secondary axis can usefully group: gene-box, gene-histogram, or the
  // no-gene celltype barplot. Numerical column histogram doesn't apply.
  const showSecondaryControls =
    !minimized &&
    !!dataset &&
    isCategorical &&
    (activePlot === "gene-box" ||
      activePlot === "gene-histogram" ||
      activePlot === "celltype-barplot");
  // Secondary picked but no values selected → tell the user instead of
  // silently falling back to the ungrouped plot.
  const secondaryPickedButEmpty =
    !!secondaryColumn &&
    selectedSecondaryValues.size === 0 &&
    selectedCelltypes.size > 0 &&
    (activePlot === "gene-box" ||
      activePlot === "gene-histogram" ||
      activePlot === "celltype-barplot");
  // Per-plot header controls.
  // Boxplot caps the y-axis (expression value); horizontal barplot caps the
  // x-axis (cell count). Same input, different axis depending on plot.
  const showYMax =
    !minimized &&
    (activePlot === "gene-box" || activePlot === "celltype-barplot");
  const valueAxisLabel = activePlot === "celltype-barplot" ? "Xmax" : "Ymax";
  const showOutliersToggle = !minimized && activePlot === "gene-box";
  const showDensityToggle =
    !minimized &&
    activePlot === "gene-histogram" &&
    !!secondaryColumn &&
    selectedSecondaryValues.size > 0 &&
    selectedCelltypes.size > 0;
  const showDownload =
    !minimized &&
    (activePlot === "gene-box" ||
      activePlot === "celltype-barplot");

  const downloadCsv = async () => {
    if (!dataset || !selectedColumn) return;
    const grouping =
      !!secondaryColumn &&
      selectedSecondaryValues.size > 0 &&
      selectedCelltypes.size > 0;
    let csv = "";
    if (activePlot === "gene-box" && selectedGene) {
      csv = await buildBoxplotCsv(
        dataset,
        selectedColumn,
        selectedGene,
        grouping ? secondaryColumn : null,
        grouping ? selectedSecondaryValues : null,
        grouping ? selectedCelltypes : null,
        topN,
      );
    } else if (activePlot === "celltype-barplot") {
      csv = buildBarplotCsv(
        dataset,
        selectedColumn,
        grouping ? secondaryColumn : null,
        grouping ? selectedSecondaryValues : null,
        grouping ? selectedCelltypes : null,
        topN,
      );
    }
    if (!csv) return;
    const blob = new Blob([csv], { type: "text/csv;charset=utf-8" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    const stamp = new Date().toISOString().slice(0, 19).replace(/[:T]/g, "-");
    const base =
      activePlot === "gene-box"
        ? `${selectedGene}_by_${selectedColumn}`
        : `cells_per_${selectedColumn}`;
    a.download = `${base}_${stamp}.csv`;
    document.body.appendChild(a);
    a.click();
    a.remove();
    URL.revokeObjectURL(url);
  };

  // Title reflects what is being shown
  let title = "Plot";
  if (!dataset) {
    title = "Plot — load a dataset";
  } else if (activePlot === "empty") {
    title = "Plot — select a column or gene";
  } else if (activePlot === "gene-box") {
    title = `${selectedGene} by ${selectedColumn}`;
  } else if (activePlot === "gene-histogram") {
    title = `Distribution · ${selectedGene}`;
  } else if (activePlot === "celltype-barplot") {
    title = `Cells per category · ${selectedColumn}`;
  } else if (activePlot === "numerical-histogram") {
    title = `Distribution · ${selectedColumn}`;
  }

  return (
    <div className="fixed inset-0 pointer-events-none z-[10]">
    <Rnd
      bounds="window"
      className="pointer-events-auto"
      dragHandleClassName="plot-panel-drag-handle"
      enableResizing={!minimized}
      minHeight={minimized ? HEADER_H : MIN_H}
      minWidth={MIN_W}
      position={position}
      size={{ width: effectiveWidth, height: effectiveHeight }}
      onDragStop={(_, d) => {
        setOffsets({
          left: d.x,
          bottom: window.innerHeight - d.y - effectiveHeight,
        });
      }}
      onResizeStop={(_, __, ref, ___, pos) => {
        const newW = parseInt(ref.style.width, 10);
        const newH = parseInt(ref.style.height, 10);
        setSizeFrac({
          width: newW / window.innerWidth,
          height: newH / window.innerHeight,
        });
        setOffsets({
          left: pos.x,
          bottom: window.innerHeight - pos.y - newH,
        });
        // Plotly's useResizeHandler only listens to window resize. Tell it to
        // reflow now that the panel has settled at its new size.
        window.dispatchEvent(new Event("resize"));
      }}
    >
      <div className="w-full h-full flex flex-col rounded-xl overflow-hidden border border-white/15 shadow-2xl backdrop-blur-md bg-background/85">
        {/* Header */}
        <div
          className="plot-panel-drag-handle flex items-center gap-2 h-9 px-2 border-b border-white/10 cursor-move select-none"
          style={{ height: HEADER_H }}
        >
          <svg
            className="w-3.5 h-3.5 text-default-400 flex-shrink-0"
            fill="none"
            stroke="currentColor"
            strokeWidth={2}
            viewBox="0 0 24 24"
          >
            <path d="M3 3v18h18" strokeLinecap="round" strokeLinejoin="round" />
            <rect x="6" y="11" width="3" height="7" />
            <rect x="11" y="7" width="3" height="11" />
            <rect x="16" y="13" width="3" height="5" />
          </svg>
          <span
            className="text-xs font-medium text-default-700 truncate flex-1"
            title={title}
          >
            {title}
          </span>

          {showViewToggle && (
            <div
              className="flex items-center text-[10px] rounded overflow-hidden border border-default-300/40"
              onMouseDown={(e) => e.stopPropagation()}
            >
              <button
                aria-pressed={plotView === "box"}
                className={`px-2 py-0.5 ${plotView === "box" ? "bg-primary/30 text-default-800" : "bg-default-100/40 text-default-500 hover:bg-default-100/70"}`}
                type="button"
                onClick={() => setPlotView("box")}
              >
                Box
              </button>
              <button
                aria-pressed={plotView === "histogram"}
                className={`px-2 py-0.5 ${plotView === "histogram" ? "bg-primary/30 text-default-800" : "bg-default-100/40 text-default-500 hover:bg-default-100/70"}`}
                type="button"
                onClick={() => setPlotView("histogram")}
              >
                Hist
              </button>
            </div>
          )}

          {showTopN && (
            <div
              className="flex items-center gap-1 text-[10px] text-default-500"
              onMouseDown={(e) => e.stopPropagation()}
            >
              <span>Top</span>
              <input
                aria-label="Top N"
                className="w-10 px-1 py-0.5 rounded bg-default-100/70 border border-default-300/40 text-xs text-center outline-none focus:border-primary"
                max={MAX_TOP_N}
                min={1}
                type="number"
                value={topN}
                onChange={(e) => {
                  const v = parseInt(e.target.value, 10);
                  if (Number.isFinite(v))
                    setTopN(Math.min(MAX_TOP_N, Math.max(1, v)));
                }}
              />
            </div>
          )}

          {showYMax && (
            <div
              className="flex items-center gap-1 text-[10px] text-default-500"
              onMouseDown={(e) => e.stopPropagation()}
            >
              <span>{valueAxisLabel}</span>
              <input
                aria-label={`${valueAxisLabel} value-axis cap`}
                className="w-20 px-1.5 py-0.5 rounded bg-default-100/70 border border-default-300/40 text-xs text-center outline-none focus:border-primary"
                placeholder="auto"
                title="Type a number, or shorthand like 20k or 1.5m."
                type="text"
                value={yMaxInput}
                onChange={(e) => setYMaxInput(e.target.value)}
              />
            </div>
          )}

          {showOutliersToggle && (
            <button
              aria-pressed={showOutliers}
              className={`text-[10px] px-2 py-0.5 rounded border ${showOutliers ? "border-primary/60 bg-primary/20 text-default-800" : "border-default-300/40 bg-default-100/40 text-default-500 hover:bg-default-100/70"}`}
              type="button"
              onMouseDown={(e) => e.stopPropagation()}
              onClick={() => setShowOutliers((v) => !v)}
            >
              Outliers
            </button>
          )}

          {showDensityToggle && (
            <button
              aria-pressed={density}
              className={`text-[10px] px-2 py-0.5 rounded border ${density ? "border-primary/60 bg-primary/20 text-default-800" : "border-default-300/40 bg-default-100/40 text-default-500 hover:bg-default-100/70"}`}
              type="button"
              onMouseDown={(e) => e.stopPropagation()}
              onClick={() => setDensity((v) => !v)}
            >
              Density
            </button>
          )}

          {showDownload && (
            <button
              aria-label="Download CSV"
              className="p-1 rounded hover:bg-white/10 text-default-500"
              type="button"
              onMouseDown={(e) => e.stopPropagation()}
              onClick={() => {
                void downloadCsv();
              }}
            >
              <svg
                className="w-3.5 h-3.5"
                fill="none"
                stroke="currentColor"
                strokeWidth={2}
                viewBox="0 0 24 24"
              >
                <path
                  d="M12 4v12m0 0l-4-4m4 4l4-4M4 20h16"
                  strokeLinecap="round"
                  strokeLinejoin="round"
                />
              </svg>
            </button>
          )}

          <button
            aria-label={minimized ? "Restore" : "Minimize"}
            className="p-1 rounded hover:bg-white/10 text-default-500"
            type="button"
            onMouseDown={(e) => e.stopPropagation()}
            onClick={() => setMinimized((m) => !m)}
          >
            {minimized ? (
              <svg
                className="w-3 h-3"
                fill="none"
                stroke="currentColor"
                strokeWidth={2.5}
                viewBox="0 0 24 24"
              >
                <path
                  d="M4 14h16M4 14l4-4M4 14l4 4"
                  strokeLinecap="round"
                  strokeLinejoin="round"
                />
              </svg>
            ) : (
              <svg
                className="w-3 h-3"
                fill="none"
                stroke="currentColor"
                strokeWidth={2.5}
                viewBox="0 0 24 24"
              >
                <path d="M5 12h14" strokeLinecap="round" />
              </svg>
            )}
          </button>
          <button
            aria-label="Close"
            className="p-1 rounded hover:bg-white/10 text-default-500"
            type="button"
            onMouseDown={(e) => e.stopPropagation()}
            onClick={() => setPlotPanelOpen(false)}
          >
            <svg
              className="w-3 h-3"
              fill="none"
              stroke="currentColor"
              strokeWidth={2.5}
              viewBox="0 0 24 24"
            >
              <path d="M6 18L18 6M6 6l12 12" strokeLinecap="round" />
            </svg>
          </button>
        </div>

        {showSecondaryControls && dataset && (
          <SecondaryGroupControls
            clusterVersion={clusterVersion}
            columnTypeOverrides={columnTypeOverrides}
            dataset={dataset}
            incrementClusterVersion={incrementClusterVersion}
            primaryColumn={selectedColumn}
            secondaryColumn={secondaryColumn}
            secondaryPaletteOverrides={secondaryPaletteOverrides}
            selectedSecondaryValues={selectedSecondaryValues}
            setSecondaryColumn={setSecondaryColumn}
            setSecondaryPaletteOverride={setSecondaryPaletteOverride}
            toggleSecondaryValue={toggleSecondaryValue}
          />
        )}

        {/* Body */}
        {!minimized && (
          <div className="flex-1 min-h-0 p-2">
            {activePlot === "empty" || !dataset ? (
              <div className="w-full h-full flex items-center justify-center text-xs text-default-400">
                Select a cluster column or gene to plot
              </div>
            ) : secondaryPickedButEmpty ? (
              <div className="w-full h-full flex items-center justify-center text-xs text-default-400 px-4 text-center">
                Select one or more {secondaryColumn} values to compare.
              </div>
            ) : activePlot === "gene-box" ? (
              <GeneCelltypeBoxplot
                clusterVersion={clusterVersion}
                column={selectedColumn!}
                dataset={dataset}
                gene={selectedGene!}
                secondaryColumn={secondaryColumn}
                secondaryPaletteOverrides={secondaryPaletteOverrides}
                selectedCelltypes={selectedCelltypes}
                selectedSecondaryValues={selectedSecondaryValues}
                showOutliers={showOutliers}
                topN={topN}
                yMax={yMax}
                onBarDoubleClick={(name) => toggleCelltype(name)}
              />
            ) : activePlot === "gene-histogram" ? (
              <GeneHistogram
                clusterVersion={clusterVersion}
                colormap={colormap}
                dataset={dataset}
                density={density}
                gene={selectedGene!}
                primaryColumn={selectedColumn}
                secondaryColumn={secondaryColumn}
                secondaryPaletteOverrides={secondaryPaletteOverrides}
                selectedCelltypes={selectedCelltypes}
                selectedSecondaryValues={selectedSecondaryValues}
              />
            ) : activePlot === "celltype-barplot" ? (
              <CelltypeBarplot
                clusterVersion={clusterVersion}
                column={selectedColumn!}
                dataset={dataset}
                secondaryColumn={secondaryColumn}
                secondaryPaletteOverrides={secondaryPaletteOverrides}
                selectedCelltypes={selectedCelltypes}
                selectedSecondaryValues={selectedSecondaryValues}
                topN={topN}
                xMax={yMax}
                onBarDoubleClick={(name) => toggleCelltype(name)}
              />
            ) : activePlot === "numerical-histogram" ? (
              <NumericalHistogram
                clusterVersion={clusterVersion}
                colormap={colormap}
                column={selectedColumn!}
                dataset={dataset}
              />
            ) : (
              <div className="w-full h-full flex items-center justify-center text-xs text-default-400">
                Column not loaded
              </div>
            )}
          </div>
        )}
      </div>
    </Rnd>
    </div>
  );
}
