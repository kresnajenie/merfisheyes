"use client";

import { useMemo } from "react";
import type { Data, Layout, Config } from "plotly.js";

import type { StandardizedDataset } from "@/lib/StandardizedDataset";
import { getClusterValue } from "@/lib/StandardizedDataset";
import { colormapRgb } from "@/lib/utils/colormaps";

import { Plot } from "./plot-loader";

interface NumericalHistogramProps {
  dataset: StandardizedDataset;
  column: string;
  colormap: string;
  clusterVersion: number;
}

export function NumericalHistogram({
  dataset,
  column,
  colormap,
  clusterVersion,
}: NumericalHistogramProps) {
  const values = useMemo(() => {
    const cluster = dataset.clusters?.find((c) => c.column === column);
    if (!cluster) return null;

    const n = cluster.valueIndices
      ? cluster.valueIndices.length
      : cluster.values.length;
    const arr = new Float32Array(n);
    for (let i = 0; i < n; i++) arr[i] = Number(getClusterValue(cluster, i));
    return arr;
  }, [dataset, column, clusterVersion]);

  const barColor = useMemo(() => {
    const [r, g, b] = colormapRgb(colormap, 0.7);
    return `rgb(${Math.round(r * 255)},${Math.round(g * 255)},${Math.round(b * 255)})`;
  }, [colormap]);

  const data = useMemo<Data[]>(() => {
    if (!values) return [];
    return [
      {
        type: "histogram",
        x: Array.from(values),
        marker: { color: barColor, line: { width: 0 } },
        nbinsx: 60,
        hovertemplate: "%{x}<br>Count: %{y:,}<extra></extra>",
      },
    ];
  }, [values, barColor]);

  const layout = useMemo<Partial<Layout>>(
    () => ({
      autosize: true,
      margin: { l: 56, r: 24, t: 8, b: 36 },
      paper_bgcolor: "rgba(0,0,0,0)",
      plot_bgcolor: "rgba(0,0,0,0)",
      font: { color: "#ddd", size: 11 },
      xaxis: {
        title: { text: column, standoff: 4 },
        gridcolor: "rgba(255,255,255,0.08)",
        zerolinecolor: "rgba(255,255,255,0.2)",
        automargin: true,
        fixedrange: true,
      },
      yaxis: {
        title: { text: "# cells", standoff: 4 },
        gridcolor: "rgba(255,255,255,0.08)",
        zerolinecolor: "rgba(255,255,255,0.2)",
        automargin: true,
        fixedrange: true,
      },
      showlegend: false,
      bargap: 0.05,
      hoverlabel: {
        bgcolor: "rgba(15,15,17,0.92)",
        bordercolor: "rgba(255,255,255,0.2)",
        font: { color: "#f5f5f5", size: 11, family: "inherit" },
        align: "left",
      },
    }),
    [column],
  );

  const config = useMemo<Partial<Config>>(
    () => ({
      responsive: true,
      displayModeBar: false,
    }),
    [],
  );

  if (!values) {
    return (
      <div className="w-full h-full flex items-center justify-center text-xs text-default-400">
        Column not loaded
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
