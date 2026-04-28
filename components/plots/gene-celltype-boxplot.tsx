"use client";

import { useEffect, useMemo, useRef, useState } from "react";
import type { Data, Layout, Config } from "plotly.js";

import type { StandardizedDataset } from "@/lib/StandardizedDataset";
import { getColorFromPalette } from "@/lib/utils/color-palette";

import { Plot } from "./plot-loader";

interface GeneCelltypeBoxplotProps {
  dataset: StandardizedDataset;
  gene: string;
  column: string;
  selectedCelltypes: Set<string>;
  topN: number;
  clusterVersion: number;
  onBarDoubleClick: (celltype: string) => void;
  // Optional secondary grouping (compare gene expression across e.g.
  // treatments within each starred celltype). Active when all three are set
  // and selectedCelltypes is non-empty.
  secondaryColumn?: string | null;
  selectedSecondaryValues?: Set<string>;
  secondaryPaletteOverrides?: Record<string, string>;
  // User-controllable y-axis cap. null = auto-range.
  yMax?: number | null;
  // When true, render `boxpoints: "outliers"` (Plotly's tukey outlier
  // markers) instead of suppressing points entirely.
  showOutliers?: boolean;
}

const DOUBLE_CLICK_MS = 350;

function quantile(sorted: Float32Array, q: number): number {
  const n = sorted.length;
  if (n === 0) return 0;
  const pos = (n - 1) * q;
  const lo = Math.floor(pos);
  const hi = Math.ceil(pos);
  if (lo === hi) return sorted[lo];
  return sorted[lo] + (pos - lo) * (sorted[hi] - sorted[lo]);
}

interface BoxStats {
  values: Float32Array;
  count: number;
  q1: number;
  median: number;
  q3: number;
  max: number;
}

function computeStats(values: Float32Array): BoxStats {
  const sorted = values.slice().sort();
  const n = sorted.length;
  return {
    values,
    count: n,
    q1: quantile(sorted, 0.25),
    median: quantile(sorted, 0.5),
    q3: quantile(sorted, 0.75),
    max: n > 0 ? sorted[n - 1] : 0,
  };
}

const fmt = (v: number) =>
  v.toLocaleString(undefined, { maximumFractionDigits: 2 });

export function GeneCelltypeBoxplot({
  dataset,
  gene,
  column,
  selectedCelltypes,
  topN,
  clusterVersion,
  onBarDoubleClick,
  secondaryColumn,
  selectedSecondaryValues,
  secondaryPaletteOverrides,
  yMax,
  showOutliers,
}: GeneCelltypeBoxplotProps) {
  const lastClickRef = useRef<{ name: string; time: number } | null>(null);
  const [expression, setExpression] = useState<number[] | null>(null);
  const [loading, setLoading] = useState(false);

  useEffect(() => {
    let cancelled = false;
    setLoading(true);
    dataset
      .getGeneExpression(gene)
      .then((vals) => {
        if (!cancelled) {
          setExpression(vals);
          setLoading(false);
        }
      })
      .catch(() => {
        if (!cancelled) {
          setExpression(null);
          setLoading(false);
        }
      });
    return () => {
      cancelled = true;
    };
  }, [dataset, gene]);

  const secondaryActive =
    !!secondaryColumn &&
    !!selectedSecondaryValues &&
    selectedSecondaryValues.size > 0 &&
    selectedCelltypes.size > 0 &&
    !!dataset.clusters?.find((c) => c.column === secondaryColumn);

  // ---------- Ungrouped path: single-axis boxplot per celltype ----------
  const groups = useMemo(() => {
    if (!expression || secondaryActive) return null;
    const cluster = dataset.clusters?.find((c) => c.column === column);
    if (!cluster) return null;

    const buckets = new Map<string, number[]>();
    if (cluster.valueIndices && cluster.uniqueValues) {
      const indices = cluster.valueIndices;
      const unique = cluster.uniqueValues;
      for (let i = 0; i < unique.length; i++) buckets.set(unique[i], []);
      for (let i = 0; i < indices.length; i++) {
        buckets.get(unique[indices[i]])!.push(i);
      }
    } else {
      for (let i = 0; i < cluster.values.length; i++) {
        const k = String(cluster.values[i]);
        const arr = buckets.get(k);
        if (arr) arr.push(i);
        else buckets.set(k, [i]);
      }
    }

    const palette = cluster.palette ?? {};
    return Array.from(buckets.entries()).map(([name, idxs], i) => {
      const vals = new Float32Array(idxs.length);
      for (let k = 0; k < idxs.length; k++) vals[k] = expression[idxs[k]] ?? 0;
      const stats = computeStats(vals);
      return {
        name,
        ...stats,
        color: palette[name] ?? getColorFromPalette(i),
        isSelected: selectedCelltypes.has(name),
      };
    });
  }, [expression, dataset, column, selectedCelltypes, clusterVersion, secondaryActive]);

  const visible = useMemo(() => {
    if (!groups) return [];
    const selected = groups
      .filter((r) => r.isSelected)
      .sort((a, b) => b.median - a.median);
    const nonSelected = groups
      .filter((r) => !r.isSelected)
      .sort((a, b) => b.median - a.median)
      .slice(0, Math.max(0, topN));
    return [...selected, ...nonSelected];
  }, [groups, topN]);

  // ---------- Grouped path: boxes per (primary × secondary) pair ----------
  const grouped = useMemo(() => {
    if (!expression || !secondaryActive || !secondaryColumn) return null;
    const primaryCluster = dataset.clusters?.find((c) => c.column === column);
    const secondaryCluster = dataset.clusters?.find(
      (c) => c.column === secondaryColumn,
    );
    if (!primaryCluster || !secondaryCluster) return null;

    const cellCount = primaryCluster.valueIndices
      ? primaryCluster.valueIndices.length
      : primaryCluster.values.length;

    // Per-cell label lookup for both columns.
    const primaryAt = (i: number): string =>
      primaryCluster.valueIndices && primaryCluster.uniqueValues
        ? primaryCluster.uniqueValues[primaryCluster.valueIndices[i]]
        : String(primaryCluster.values[i]);
    const secondaryAt = (i: number): string =>
      secondaryCluster.valueIndices && secondaryCluster.uniqueValues
        ? secondaryCluster.uniqueValues[secondaryCluster.valueIndices[i]]
        : String(secondaryCluster.values[i]);

    // Bucket cell indices per (primary, secondary) for starred + selected.
    const bucket = new Map<string, Map<string, number[]>>();
    for (let i = 0; i < cellCount; i++) {
      const p = primaryAt(i);
      if (!selectedCelltypes.has(p)) continue;
      const s = secondaryAt(i);
      if (!selectedSecondaryValues!.has(s)) continue;
      let inner = bucket.get(p);
      if (!inner) {
        inner = new Map();
        bucket.set(p, inner);
      }
      let arr = inner.get(s);
      if (!arr) {
        arr = [];
        inner.set(s, arr);
      }
      arr.push(i);
    }

    // Stats per pair.
    const pairs: Array<{ primary: string; secondary: string } & BoxStats> = [];
    for (const [primary, inner] of bucket) {
      for (const [secondary, idxs] of inner) {
        const vals = new Float32Array(idxs.length);
        for (let k = 0; k < idxs.length; k++)
          vals[k] = expression[idxs[k]] ?? 0;
        pairs.push({ primary, secondary, ...computeStats(vals) });
      }
    }

    // Sort primaries by mean of medians across selected secondaries.
    const aggMedian = new Map<string, number[]>();
    for (const p of pairs) {
      if (!aggMedian.has(p.primary)) aggMedian.set(p.primary, []);
      aggMedian.get(p.primary)!.push(p.median);
    }
    const primaryOrder = Array.from(selectedCelltypes)
      .filter((p) => aggMedian.has(p))
      .sort((a, b) => {
        const ma = aggMedian.get(a)!;
        const mb = aggMedian.get(b)!;
        const avgA = ma.reduce((s, v) => s + v, 0) / ma.length;
        const avgB = mb.reduce((s, v) => s + v, 0) / mb.length;
        return avgB - avgA;
      });

    // Secondary value order from the column's uniqueValues, filtered by
    // selection (per Q30).
    const allUnique =
      secondaryCluster.uniqueValues ??
      Array.from(new Set(secondaryCluster.values.map(String)));
    const secondaryOrder = allUnique.filter((v) =>
      selectedSecondaryValues!.has(v),
    );

    // Fill missing pairs with empty placeholders so the grid stays aligned
    // (per Q14).
    const lookup = new Map<string, BoxStats>();
    for (const p of pairs) lookup.set(`${p.primary} ${p.secondary}`, p);
    const allPairs: Array<{ primary: string; secondary: string } & BoxStats> =
      [];
    for (const primary of primaryOrder) {
      for (const secondary of secondaryOrder) {
        const found = lookup.get(`${primary} ${secondary}`);
        if (found) {
          allPairs.push({ primary, secondary, ...found });
        } else {
          allPairs.push({
            primary,
            secondary,
            values: new Float32Array(0),
            count: 0,
            q1: 0,
            median: 0,
            q3: 0,
            max: 0,
          });
        }
      }
    }

    // Secondary palette: column palette + per-value override.
    const secPalette = secondaryCluster.palette ?? {};
    const colorFor = (v: string, idx: number) =>
      secondaryPaletteOverrides?.[v] ??
      secPalette[v] ??
      getColorFromPalette(idx);

    return {
      primaryOrder,
      secondaryOrder,
      pairs: allPairs,
      colorFor,
    };
  }, [
    expression,
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
    // Grouped data
    if (grouped) {
      const traces: Data[] = [];
      grouped.secondaryOrder.forEach((secondary, i) => {
        const color = grouped.colorFor(secondary, i);
        const secPairs = grouped.pairs.filter((p) => p.secondary === secondary);

        // Visual box trace: all cells for this secondary, x = primary per cell.
        const xs: string[] = [];
        const ys: number[] = [];
        for (const p of secPairs) {
          for (const v of p.values) {
            xs.push(p.primary);
            ys.push(v);
          }
        }
        const boxTrace: Record<string, unknown> = {
          type: "box",
          name: secondary,
          x: xs,
          y: ys,
          boxpoints: showOutliers ? "outliers" : false,
          hoverinfo: "skip",
          line: { color: "rgba(255,255,255,0.85)", width: 1.5 },
          fillcolor: color,
          offsetgroup: secondary,
          legendgroup: secondary,
          showlegend: true,
        };
        if (showOutliers) {
          boxTrace.marker = {
            color,
            size: 3,
            opacity: 0.7,
            line: { width: 0 },
          };
        }
        traces.push(boxTrace as unknown as Data);

        // Hover overlay: one scatter point per primary at the median (or 0
        // for empty pairs, so the placeholder stays visible at the baseline).
        // Empty pairs show a faint dash so the slot stays visible.
        traces.push({
          type: "scatter",
          mode: "markers",
          x: secPairs.map((p) => p.primary),
          y: secPairs.map((p) => p.median),
          marker: {
            size: secPairs.map((p) => (p.count === 0 ? 12 : 24)),
            color: secPairs.map((p) =>
              p.count === 0 ? "rgba(255,255,255,0.35)" : "rgba(0,0,0,0)",
            ),
            symbol: secPairs.map((p) => (p.count === 0 ? "line-ew" : "circle")),
            line: {
              width: secPairs.map((p) => (p.count === 0 ? 1 : 0)),
              color: "rgba(255,255,255,0.35)",
            },
          },
          offsetgroup: secondary,
          legendgroup: secondary,
          showlegend: false,
          hovertemplate: secPairs.map((p) =>
            p.count === 0
              ? `${secondaryColumn}: ${secondary}<br>no cells<extra></extra>`
              : `${secondaryColumn}: ${secondary}<br>` +
                `n: ${p.count.toLocaleString()}<br>` +
                `max: ${fmt(p.max)}<br>` +
                `q3: ${fmt(p.q3)}<br>` +
                `median: ${fmt(p.median)}<br>` +
                `q1: ${fmt(p.q1)}<extra></extra>`,
          ),
        } as unknown as Data);
      });
      return traces;
    }

    // Ungrouped data (existing behaviour)
    if (visible.length === 0) return [];
    const hasSelection = selectedCelltypes.size > 0;
    const boxes: Data[] = visible.map((r) => {
      const trace: Record<string, unknown> = {
        type: "box",
        name: r.name,
        y: Array.from(r.values),
        boxpoints: showOutliers ? "outliers" : false,
        hoverinfo: "skip",
        line: { color: "rgba(255,255,255,0.85)", width: 1.5 },
        fillcolor: r.color,
        opacity: !hasSelection || r.isSelected ? 1.0 : 0.25,
        showlegend: false,
      };
      if (showOutliers) {
        trace.marker = {
          color: r.color,
          size: 3,
          opacity: 0.7,
          line: { width: 0 },
        };
      }
      return trace as unknown as Data;
    });

    const hoverTrace: Data = {
      type: "scatter",
      mode: "markers",
      x: visible.map((r) => r.name),
      y: visible.map((r) => r.median),
      marker: {
        size: 24,
        color: "rgba(0,0,0,0)",
        line: { width: 0 },
      },
      hovertemplate: visible.map(
        (r) =>
          `n: ${r.count.toLocaleString()}<br>` +
          `max: ${fmt(r.max)}<br>` +
          `q3: ${fmt(r.q3)}<br>` +
          `median: ${fmt(r.median)}<br>` +
          `q1: ${fmt(r.q1)}` +
          `<extra></extra>`,
      ),
      showlegend: false,
    };

    return [...boxes, hoverTrace];
  }, [visible, selectedCelltypes.size, grouped, secondaryColumn, showOutliers]);

  const layout = useMemo<Partial<Layout>>(() => {
    const xCategoryArray = grouped
      ? grouped.primaryOrder
      : visible.map((r) => r.name);
    return {
      autosize: true,
      margin: { l: 56, r: 16, t: grouped ? 24 : 8, b: 80 },
      paper_bgcolor: "rgba(0,0,0,0)",
      plot_bgcolor: "rgba(0,0,0,0)",
      font: { color: "#ddd", size: 11 },
      xaxis: {
        gridcolor: "rgba(255,255,255,0)",
        tickangle: -35,
        automargin: true,
        tickfont: { size: 10 },
        categoryorder: "array",
        categoryarray: xCategoryArray,
        fixedrange: true,
      },
      yaxis: {
        title: { text: gene, standoff: 4 },
        gridcolor: "rgba(255,255,255,0.08)",
        zerolinecolor: "rgba(255,255,255,0.2)",
        automargin: true,
        fixedrange: true,
        ...(yMax != null
          ? { range: [0, yMax], autorange: false }
          : { autorange: true }),
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
      boxmode: grouped ? "group" : undefined,
      boxgap: 0.25,
      hovermode: "closest",
      hoverlabel: {
        bgcolor: "rgba(15,15,17,0.92)",
        bordercolor: "rgba(255,255,255,0.2)",
        font: { color: "#f5f5f5", size: 11, family: "inherit" },
        align: "left",
      },
    };
  }, [gene, visible, grouped, yMax]);

  const config = useMemo<Partial<Config>>(
    () => ({
      responsive: true,
      displayModeBar: false,
      doubleClick: false,
    }),
    [],
  );

  if (loading) {
    return (
      <div className="w-full h-full flex items-center justify-center text-xs text-default-400">
        Loading {gene}…
      </div>
    );
  }
  if (!expression) {
    return (
      <div className="w-full h-full flex items-center justify-center text-xs text-default-400">
        Gene not found
      </div>
    );
  }
  if (!grouped && visible.length === 0) {
    return (
      <div className="w-full h-full flex items-center justify-center text-xs text-default-400">
        No data
      </div>
    );
  }
  if (grouped && grouped.pairs.length === 0) {
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
        if (!pt) return;
        const fromScatter = typeof pt.x === "string" ? pt.x : null;
        const fromBox = (pt as { data?: { name?: string } }).data?.name ?? null;
        // In grouped mode, both are set, but we want the primary celltype name
        // (the x category), not the secondary trace name. fromScatter wins.
        const name = fromScatter ?? fromBox;
        if (!name) return;
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
