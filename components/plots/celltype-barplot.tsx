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
  // Bumped by the store when cluster columns load; included in memo deps
  // so a freshly-loaded column produces fresh counts.
  clusterVersion: number;
  onBarDoubleClick: (celltype: string) => void;
  // Optional secondary grouping (compare counts across e.g. treatments
  // within each starred celltype). Active when all are set and
  // selectedCelltypes is non-empty.
  secondaryColumn?: string | null;
  selectedSecondaryValues?: Set<string>;
  secondaryPaletteOverrides?: Record<string, string>;
  // User-controllable cap on the count axis (horizontal bars → x-axis).
  // null = auto-range.
  xMax?: number | null;
}

const DOUBLE_CLICK_MS = 350;

export function CelltypeBarplot({
  dataset,
  column,
  selectedCelltypes,
  topN,
  clusterVersion,
  onBarDoubleClick,
  secondaryColumn,
  selectedSecondaryValues,
  secondaryPaletteOverrides,
  xMax,
}: CelltypeBarplotProps) {
  const lastClickRef = useRef<{ name: string; time: number } | null>(null);

  const secondaryActive =
    !!secondaryColumn &&
    !!selectedSecondaryValues &&
    selectedSecondaryValues.size > 0 &&
    selectedCelltypes.size > 0 &&
    !!dataset.clusters?.find((c) => c.column === secondaryColumn);

  // ---------- Ungrouped path: single bar per celltype ----------
  const counts = useMemo(() => {
    if (secondaryActive) return null;
    const cluster = dataset.clusters?.find((c) => c.column === column);
    if (!cluster) return null;

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
    return Array.from(map.entries()).map(([name, count], idx) => ({
      name,
      count,
      color: palette[name] ?? getColorFromPalette(idx),
      isSelected: selectedCelltypes.has(name),
    }));
  }, [dataset, column, selectedCelltypes, clusterVersion, secondaryActive]);

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
    return [...selected, ...nonSelected];
  }, [counts, topN]);

  // ---------- Grouped path: counts per (primary × secondary) pair ----------
  const grouped = useMemo(() => {
    if (!secondaryActive || !secondaryColumn) return null;
    const primaryCluster = dataset.clusters?.find((c) => c.column === column);
    const secondaryCluster = dataset.clusters?.find(
      (c) => c.column === secondaryColumn,
    );
    if (!primaryCluster || !secondaryCluster) return null;

    const cellCount = primaryCluster.valueIndices
      ? primaryCluster.valueIndices.length
      : primaryCluster.values.length;
    const primaryAt = (i: number): string =>
      primaryCluster.valueIndices && primaryCluster.uniqueValues
        ? primaryCluster.uniqueValues[primaryCluster.valueIndices[i]]
        : String(primaryCluster.values[i]);
    const secondaryAt = (i: number): string =>
      secondaryCluster.valueIndices && secondaryCluster.uniqueValues
        ? secondaryCluster.uniqueValues[secondaryCluster.valueIndices[i]]
        : String(secondaryCluster.values[i]);

    // Bucket counts per (primary, secondary) for starred + selected.
    const counts2d = new Map<string, Map<string, number>>();
    for (let i = 0; i < cellCount; i++) {
      const p = primaryAt(i);
      if (!selectedCelltypes.has(p)) continue;
      const s = secondaryAt(i);
      if (!selectedSecondaryValues!.has(s)) continue;
      let inner = counts2d.get(p);
      if (!inner) {
        inner = new Map();
        counts2d.set(p, inner);
      }
      inner.set(s, (inner.get(s) ?? 0) + 1);
    }

    // Total per primary across selected secondaries — used for sort + share %.
    const primaryTotal = new Map<string, number>();
    for (const [p, inner] of counts2d) {
      let sum = 0;
      for (const v of inner.values()) sum += v;
      primaryTotal.set(p, sum);
    }

    // Sort starred primaries by total desc.
    const primaryOrder = Array.from(selectedCelltypes)
      .filter((p) => primaryTotal.has(p))
      .sort((a, b) => (primaryTotal.get(b) ?? 0) - (primaryTotal.get(a) ?? 0));

    const allUnique =
      secondaryCluster.uniqueValues ??
      Array.from(new Set(secondaryCluster.values.map(String)));
    const secondaryOrder = allUnique.filter((v) =>
      selectedSecondaryValues!.has(v),
    );

    const secPalette = secondaryCluster.palette ?? {};
    const colorFor = (v: string, idx: number) =>
      secondaryPaletteOverrides?.[v] ??
      secPalette[v] ??
      getColorFromPalette(idx);

    return {
      primaryOrder,
      secondaryOrder,
      counts2d,
      primaryTotal,
      colorFor,
    };
  }, [
    dataset,
    column,
    secondaryColumn,
    selectedCelltypes,
    selectedSecondaryValues,
    secondaryPaletteOverrides,
    secondaryActive,
    clusterVersion,
  ]);

  const data = useMemo<Data[]>(() => {
    if (grouped) {
      // Horizontal grouped bars. y = primary, x = count, one trace per
      // selected secondary value. Reverse so highest total sits at top.
      const yOrdered = [...grouped.primaryOrder].reverse();
      return grouped.secondaryOrder.map((secondary, i) => {
        const color = grouped.colorFor(secondary, i);
        const xs = yOrdered.map(
          (p) => grouped.counts2d.get(p)?.get(secondary) ?? 0,
        );
        const customdata = yOrdered.map((p) => {
          const tot = grouped.primaryTotal.get(p) ?? 0;
          const c = grouped.counts2d.get(p)?.get(secondary) ?? 0;
          return tot > 0 ? ((c / tot) * 100).toFixed(2) : "0";
        });
        return {
          type: "bar",
          orientation: "h",
          name: secondary,
          x: xs,
          y: yOrdered,
          marker: { color, line: { width: 0 } },
          offsetgroup: secondary,
          legendgroup: secondary,
          showlegend: true,
          customdata,
          hovertemplate:
            `<b>%{y}</b><br>` +
            `${secondaryColumn}: ${secondary}<br>` +
            `Cells: %{x:,}<br>` +
            `Share: %{customdata}%<extra></extra>`,
        } as unknown as Data;
      });
    }

    // Ungrouped data
    if (visibleRows.length === 0) return [];
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
  }, [visibleRows, selectedCelltypes.size, total, grouped, secondaryColumn]);

  // Stable y-axis category order. Default Plotly behaviour can re-sort
  // alphabetically for some inputs (e.g. when uniqueValues comes pre-sorted
  // from a worker); pin it explicitly via categoryarray.
  const yCategoryArray = useMemo(() => {
    if (grouped) return [...grouped.primaryOrder].reverse();
    return [...visibleRows].reverse().map((r) => r.name);
  }, [grouped, visibleRows]);

  const layout = useMemo<Partial<Layout>>(
    () => ({
      autosize: true,
      margin: { l: 120, r: 24, t: grouped ? 24 : 8, b: 36 },
      paper_bgcolor: "rgba(0,0,0,0)",
      plot_bgcolor: "rgba(0,0,0,0)",
      font: { color: "#ddd", size: 11 },
      xaxis: {
        title: { text: "# cells", standoff: 4 },
        gridcolor: "rgba(255,255,255,0.08)",
        zerolinecolor: "rgba(255,255,255,0.2)",
        automargin: true,
        fixedrange: true,
        ...(xMax != null
          ? { range: [0, xMax], autorange: false }
          : { autorange: true }),
      },
      yaxis: {
        gridcolor: "rgba(255,255,255,0)",
        automargin: true,
        tickfont: { size: 10 },
        fixedrange: true,
        categoryorder: "array",
        categoryarray: yCategoryArray,
      },
      showlegend: !!grouped,
      legend: grouped
        ? {
            orientation: "h",
            y: 1.05,
            x: 0,
            xanchor: "left",
            yanchor: "bottom",
            font: { size: 10 },
          }
        : undefined,
      barmode: grouped ? "group" : undefined,
      bargap: 0.15,
      hoverlabel: {
        bgcolor: "rgba(15,15,17,0.92)",
        bordercolor: "rgba(255,255,255,0.2)",
        font: { color: "#f5f5f5", size: 11, family: "inherit" },
        align: "left",
      },
    }),
    [grouped, xMax, yCategoryArray],
  );

  const config = useMemo<Partial<Config>>(
    () => ({
      responsive: true,
      displayModeBar: false,
      doubleClick: false,
    }),
    [],
  );

  if (!grouped && !counts) {
    return (
      <div className="w-full h-full flex items-center justify-center text-xs text-default-400">
        Column not loaded
      </div>
    );
  }
  if (!grouped && counts && counts.length === 0) {
    return (
      <div className="w-full h-full flex items-center justify-center text-xs text-default-400">
        No data
      </div>
    );
  }
  if (grouped && grouped.primaryOrder.length === 0) {
    return (
      <div className="w-full h-full flex items-center justify-center text-xs text-default-400">
        Star a celltype to compare across {secondaryColumn}
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
