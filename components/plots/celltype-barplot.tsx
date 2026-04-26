"use client";

import { useMemo, useRef } from "react";
import type { Data, Layout, Config } from "plotly.js";

import type { StandardizedDataset } from "@/lib/StandardizedDataset";
import { getColorFromPalette } from "@/lib/utils/color-palette";

import { Plot } from "./plot-loader";

interface CelltypeBarplotProps {
  dataset: StandardizedDataset;
  column: string;
  selectedCelltypes: Set<string>;
  topN: number;
  onBarDoubleClick: (celltype: string) => void;
}

const DOUBLE_CLICK_MS = 350;

export function CelltypeBarplot({
  dataset,
  column,
  selectedCelltypes,
  topN,
  onBarDoubleClick,
}: CelltypeBarplotProps) {
  const lastClickRef = useRef<{ name: string; time: number } | null>(null);
  const counts = useMemo(() => {
    const cluster = dataset.clusters?.find((c) => c.column === column);
    if (!cluster) return null;

    // Count cells per category. Prefer indexed lookup if available.
    const map = new Map<string, number>();
    if (cluster.valueIndices && cluster.uniqueValues) {
      const counts32 = new Uint32Array(cluster.uniqueValues.length);
      for (let i = 0; i < cluster.valueIndices.length; i++) {
        counts32[cluster.valueIndices[i]]++;
      }
      cluster.uniqueValues.forEach((v, idx) => map.set(v, counts32[idx]));
    } else {
      for (const v of cluster.values) {
        const k = String(v);
        map.set(k, (map.get(k) ?? 0) + 1);
      }
    }

    const palette = cluster.palette ?? {};
    const allRows = Array.from(map.entries()).map(([name, count], idx) => ({
      name,
      count,
      color: palette[name] ?? getColorFromPalette(idx),
      isSelected: selectedCelltypes.has(name),
    }));

    return allRows;
  }, [dataset, column, selectedCelltypes]);

  const total = useMemo(
    () => (counts ? counts.reduce((s, r) => s + r.count, 0) : 0),
    [counts],
  );

  const visibleRows = useMemo(() => {
    if (!counts) return [];
    const selected = counts
      .filter((r) => r.isSelected)
      .sort((a, b) => b.count - a.count);
    const nonSelected = counts
      .filter((r) => !r.isSelected)
      .sort((a, b) => b.count - a.count)
      .slice(0, Math.max(0, topN));
    // Selected always shown on top, non-selected truncated to topN.
    return [...selected, ...nonSelected];
  }, [counts, topN]);

  const data = useMemo<Data[]>(() => {
    if (visibleRows.length === 0) return [];
    // Plotly horizontal bar: y = categories (top→bottom in array order),
    // we want highest at top, so reverse the visible order for display.
    const ordered = [...visibleRows].reverse();
    const hasSelection = selectedCelltypes.size > 0;
    return [
      {
        type: "bar",
        orientation: "h",
        x: ordered.map((r) => r.count),
        y: ordered.map((r) => r.name),
        marker: {
          color: ordered.map((r) => r.color),
          opacity: ordered.map((r) =>
            !hasSelection || r.isSelected ? 1.0 : 0.25,
          ),
          line: { width: 0 },
        },
        customdata: ordered.map((r) =>
          total > 0 ? ((r.count / total) * 100).toFixed(2) : "0",
        ),
        hovertemplate:
          "<b>%{y}</b><br>" +
          "Cells: %{x:,}<br>" +
          "Share: %{customdata}%" +
          "<extra></extra>",
      },
    ];
  }, [visibleRows, selectedCelltypes.size, total]);

  const layout = useMemo<Partial<Layout>>(
    () => ({
      autosize: true,
      margin: { l: 120, r: 24, t: 8, b: 36 },
      paper_bgcolor: "rgba(0,0,0,0)",
      plot_bgcolor: "rgba(0,0,0,0)",
      font: { color: "#ddd", size: 11 },
      xaxis: {
        title: { text: "# cells", standoff: 4 },
        gridcolor: "rgba(255,255,255,0.08)",
        zerolinecolor: "rgba(255,255,255,0.2)",
        automargin: true,
      },
      yaxis: {
        gridcolor: "rgba(255,255,255,0)",
        automargin: true,
        tickfont: { size: 10 },
      },
      showlegend: false,
      bargap: 0.15,
      hoverlabel: {
        bgcolor: "rgba(15,15,17,0.92)",
        bordercolor: "rgba(255,255,255,0.2)",
        font: { color: "#f5f5f5", size: 11, family: "inherit" },
        align: "left",
      },
    }),
    [],
  );

  const config = useMemo<Partial<Config>>(
    () => ({
      responsive: true,
      displayModeBar: false,
      doubleClick: false,
    }),
    [],
  );

  if (!counts) {
    return (
      <div className="w-full h-full flex items-center justify-center text-xs text-default-400">
        Column not loaded
      </div>
    );
  }
  if (counts.length === 0) {
    return (
      <div className="w-full h-full flex items-center justify-center text-xs text-default-400">
        No data
      </div>
    );
  }

  return (
    <Plot
      config={config}
      data={data}
      layout={layout}
      style={{ width: "100%", height: "100%" }}
      useResizeHandler
      onClick={(ev) => {
        const pt = ev.points?.[0];
        if (!pt || typeof pt.y !== "string") return;
        const name = pt.y;
        const now = performance.now();
        const last = lastClickRef.current;
        if (last && last.name === name && now - last.time < DOUBLE_CLICK_MS) {
          onBarDoubleClick(name);
          lastClickRef.current = null;
        } else {
          lastClickRef.current = { name, time: now };
        }
      }}
    />
  );
}
