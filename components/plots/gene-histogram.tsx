"use client";

import { useEffect, useMemo, useState } from "react";
import type { Data, Layout, Config } from "plotly.js";

import type { StandardizedDataset } from "@/lib/StandardizedDataset";
import { colormapRgb } from "@/lib/utils/colormaps";
import { getColorFromPalette } from "@/lib/utils/color-palette";

import { Plot } from "./plot-loader";

interface GeneHistogramProps {
  dataset: StandardizedDataset;
  gene: string;
  colormap: string;
  // Optional secondary grouping. Active when all are set and
  // selectedCelltypes is non-empty — renders one overlay histogram per
  // selected secondary value, restricted to starred celltypes.
  primaryColumn?: string | null;
  selectedCelltypes?: Set<string>;
  secondaryColumn?: string | null;
  selectedSecondaryValues?: Set<string>;
  secondaryPaletteOverrides?: Record<string, string>;
  clusterVersion?: number;
  // When true and grouping is active, normalize each group to its own
  // probability-density (matplotlib-style density=True) so groups of
  // different sizes are visually comparable.
  density?: boolean;
}

export function GeneHistogram({
  dataset,
  gene,
  colormap,
  primaryColumn,
  selectedCelltypes,
  secondaryColumn,
  selectedSecondaryValues,
  secondaryPaletteOverrides,
  density,
}: GeneHistogramProps) {
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
    !!selectedCelltypes &&
    selectedCelltypes.size > 0 &&
    !!primaryColumn &&
    !!dataset.clusters?.find((c) => c.column === secondaryColumn);

  const baseColor = useMemo(() => {
    const [r, g, b] = colormapRgb(colormap, 0.7);
    return `rgb(${Math.round(r * 255)},${Math.round(g * 255)},${Math.round(b * 255)})`;
  }, [colormap]);

  // Per-secondary expression buckets (starred celltypes only).
  const buckets = useMemo(() => {
    if (!expression || !secondaryActive) return null;
    const primaryCluster = dataset.clusters?.find(
      (c) => c.column === primaryColumn,
    );
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

    const allUnique =
      secondaryCluster.uniqueValues ??
      Array.from(new Set(secondaryCluster.values.map(String)));
    const secondaryOrder = allUnique.filter((v) =>
      selectedSecondaryValues!.has(v),
    );

    // Two passes: count then fill, so we can size Float32Arrays exactly.
    const sizes = new Map<string, number>();
    for (const v of secondaryOrder) sizes.set(v, 0);
    for (let i = 0; i < cellCount; i++) {
      const p = primaryAt(i);
      if (!selectedCelltypes!.has(p)) continue;
      const s = secondaryAt(i);
      if (!sizes.has(s)) continue;
      sizes.set(s, sizes.get(s)! + 1);
    }
    const arrays = new Map<string, Float32Array>();
    const cursors = new Map<string, number>();
    for (const v of secondaryOrder) {
      arrays.set(v, new Float32Array(sizes.get(v) ?? 0));
      cursors.set(v, 0);
    }
    for (let i = 0; i < cellCount; i++) {
      const p = primaryAt(i);
      if (!selectedCelltypes!.has(p)) continue;
      const s = secondaryAt(i);
      const arr = arrays.get(s);
      if (!arr) continue;
      const c = cursors.get(s)!;
      arr[c] = expression[i] ?? 0;
      cursors.set(s, c + 1);
    }

    const secPalette = secondaryCluster.palette ?? {};
    const colorFor = (v: string, idx: number) =>
      secondaryPaletteOverrides?.[v] ??
      secPalette[v] ??
      getColorFromPalette(idx);

    return {
      secondaryOrder,
      arrays,
      colorFor,
    };
  }, [
    expression,
    dataset,
    primaryColumn,
    secondaryColumn,
    selectedCelltypes,
    selectedSecondaryValues,
    secondaryPaletteOverrides,
    secondaryActive,
  ]);

  const data = useMemo<Data[]>(() => {
    if (!expression) return [];
    if (buckets) {
      return buckets.secondaryOrder.map((secondary, i) => {
        const color = buckets.colorFor(secondary, i);
        const arr = buckets.arrays.get(secondary)!;
        return {
          type: "histogram",
          name: secondary,
          x: Array.from(arr),
          marker: { color, line: { width: 0 } },
          opacity: 0.55,
          nbinsx: 60,
          // probability density = bin counts / (total * bin width).
          histnorm: density ? "probability density" : undefined,
          legendgroup: secondary,
          showlegend: true,
          hovertemplate:
            `${secondaryColumn}: ${secondary}<br>` +
            `${gene}: %{x}<br>` +
            (density ? "density: %{y:.4f}" : "n: %{y:,}") +
            `<extra></extra>`,
        };
      });
    }
    return [
      {
        type: "histogram",
        x: expression,
        marker: { color: baseColor, line: { width: 0 } },
        nbinsx: 60,
        hovertemplate: "%{x}<br>Count: %{y:,}<extra></extra>",
      },
    ];
  }, [expression, buckets, baseColor, gene, secondaryColumn, density]);

  const layout = useMemo<Partial<Layout>>(
    () => ({
      autosize: true,
      margin: { l: 56, r: 24, t: buckets ? 24 : 8, b: 36 },
      paper_bgcolor: "rgba(0,0,0,0)",
      plot_bgcolor: "rgba(0,0,0,0)",
      font: { color: "#ddd", size: 11 },
      xaxis: {
        title: { text: gene, standoff: 4 },
        gridcolor: "rgba(255,255,255,0.08)",
        zerolinecolor: "rgba(255,255,255,0.2)",
        automargin: true,
        fixedrange: true,
      },
      yaxis: {
        title: {
          text: buckets && density ? "density" : "# cells",
          standoff: 4,
        },
        gridcolor: "rgba(255,255,255,0.08)",
        zerolinecolor: "rgba(255,255,255,0.2)",
        automargin: true,
        fixedrange: true,
      },
      showlegend: !!buckets,
      legend: buckets
        ? {
            orientation: "h",
            y: 1.05,
            x: 0,
            xanchor: "left",
            yanchor: "bottom",
            font: { size: 10 },
          }
        : undefined,
      barmode: buckets ? "overlay" : undefined,
      bargap: 0.05,
      hoverlabel: {
        bgcolor: "rgba(15,15,17,0.92)",
        bordercolor: "rgba(255,255,255,0.2)",
        font: { color: "#f5f5f5", size: 11, family: "inherit" },
        align: "left",
      },
    }),
    [gene, buckets, density],
  );

  const config = useMemo<Partial<Config>>(
    () => ({
      responsive: true,
      displayModeBar: false,
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
  if (secondaryActive && buckets && buckets.arrays.size === 0) {
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
    />
  );
}
