"use client";

import { useEffect, useState } from "react";
import { Rnd } from "react-rnd";

import type { StandardizedDataset } from "@/lib/StandardizedDataset";
import {
  usePanelDatasetStore,
  usePanelVisualizationStore,
} from "@/lib/hooks/usePanelStores";
import { getEffectiveColumnType } from "@/lib/utils/column-type-utils";

import { CelltypeBarplot } from "./plots/celltype-barplot";
import { NumericalHistogram } from "./plots/numerical-histogram";

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
  const colormap = usePanelVisualizationStore((s) => s.colormap);

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

  // Title reflects what is being shown
  let title = "Plot";
  if (!selectedColumn || !dataset) {
    title = "Plot — select a column";
  } else if (isCategorical) {
    title = `Cells per category · ${selectedColumn}`;
  } else if (isNumerical) {
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

          {!minimized && isCategorical && (
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

        {/* Body */}
        {!minimized && (
          <div className="flex-1 min-h-0 p-2">
            {!dataset || !selectedColumn ? (
              <div className="w-full h-full flex items-center justify-center text-xs text-default-400">
                Select a cluster column to plot
              </div>
            ) : isCategorical ? (
              <CelltypeBarplot
                clusterVersion={clusterVersion}
                column={selectedColumn}
                dataset={dataset}
                selectedCelltypes={selectedCelltypes}
                topN={topN}
                onBarDoubleClick={(name) => toggleCelltype(name)}
              />
            ) : isNumerical ? (
              <NumericalHistogram
                clusterVersion={clusterVersion}
                colormap={colormap}
                column={selectedColumn}
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
