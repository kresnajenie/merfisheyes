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
}

const DOUBLE_CLICK_MS = 350;

function median(sorted: Float32Array): number {
  const n = sorted.length;
  if (n === 0) return 0;
  const mid = n >> 1;
  return n % 2 ? sorted[mid] : (sorted[mid - 1] + sorted[mid]) / 2;
}

export function GeneCelltypeBoxplot({
  dataset,
  gene,
  column,
  selectedCelltypes,
  topN,
  clusterVersion,
  onBarDoubleClick,
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

  // Group expression values per celltype, compute medians, attach palette colors.
  const groups = useMemo(() => {
    if (!expression) return null;
    const cluster = dataset.clusters?.find((c) => c.column === column);
    if (!cluster) return null;

    // Bucket cell indices per category.
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
    const rows = Array.from(buckets.entries()).map(([name, idxs], i) => {
      const vals = new Float32Array(idxs.length);
      for (let k = 0; k < idxs.length; k++) vals[k] = expression[idxs[k]] ?? 0;
      const sorted = vals.slice().sort();
      return {
        name,
        values: vals,
        median: median(sorted),
        color: palette[name] ?? getColorFromPalette(i),
        isSelected: selectedCelltypes.has(name),
      };
    });
    return rows;
  }, [expression, dataset, column, selectedCelltypes, clusterVersion]);

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

  const data = useMemo<Data[]>(() => {
    if (visible.length === 0) return [];
    const hasSelection = selectedCelltypes.size > 0;
    return visible.map((r) => ({
      type: "box",
      name: r.name,
      y: Array.from(r.values),
      boxpoints: "all",
      jitter: 0.5,
      pointpos: 0,
      marker: {
        color: r.color,
        size: 3,
        opacity: !hasSelection || r.isSelected ? 0.6 : 0.15,
        line: { width: 0 },
      },
      line: { color: r.color, width: 1 },
      fillcolor: r.color,
      opacity: !hasSelection || r.isSelected ? 1.0 : 0.25,
      hovertemplate:
        `<b>${r.name}</b><br>` +
        gene +
        ": %{y}" +
        "<extra></extra>",
    }));
  }, [visible, selectedCelltypes.size, gene]);

  const layout = useMemo<Partial<Layout>>(
    () => ({
      autosize: true,
      margin: { l: 56, r: 16, t: 8, b: 80 },
      paper_bgcolor: "rgba(0,0,0,0)",
      plot_bgcolor: "rgba(0,0,0,0)",
      font: { color: "#ddd", size: 11 },
      xaxis: {
        gridcolor: "rgba(255,255,255,0)",
        tickangle: -35,
        automargin: true,
        tickfont: { size: 10 },
      },
      yaxis: {
        title: { text: gene, standoff: 4 },
        gridcolor: "rgba(255,255,255,0.08)",
        zerolinecolor: "rgba(255,255,255,0.2)",
        automargin: true,
      },
      showlegend: false,
      boxgap: 0.25,
      hoverlabel: {
        bgcolor: "rgba(15,15,17,0.92)",
        bordercolor: "rgba(255,255,255,0.2)",
        font: { color: "#f5f5f5", size: 11, family: "inherit" },
        align: "left",
      },
    }),
    [gene],
  );

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
  if (visible.length === 0) {
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
        if (!pt) return;
        const trace = (pt as { data?: { name?: string } }).data;
        const name = trace?.name;
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
