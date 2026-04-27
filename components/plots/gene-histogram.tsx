"use client";

import { useEffect, useMemo, useState } from "react";
import type { Data, Layout, Config } from "plotly.js";

import type { StandardizedDataset } from "@/lib/StandardizedDataset";
import { colormapRgb } from "@/lib/utils/colormaps";

import { Plot } from "./plot-loader";

interface GeneHistogramProps {
  dataset: StandardizedDataset;
  gene: string;
  colormap: string;
}

export function GeneHistogram({ dataset, gene, colormap }: GeneHistogramProps) {
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

  const barColor = useMemo(() => {
    const [r, g, b] = colormapRgb(colormap, 0.7);
    return `rgb(${Math.round(r * 255)},${Math.round(g * 255)},${Math.round(b * 255)})`;
  }, [colormap]);

  const data = useMemo<Data[]>(() => {
    if (!expression) return [];
    return [
      {
        type: "histogram",
        x: expression,
        marker: { color: barColor, line: { width: 0 } },
        nbinsx: 60,
        hovertemplate: "%{x}<br>Count: %{y:,}<extra></extra>",
      },
    ];
  }, [expression, barColor]);

  const layout = useMemo<Partial<Layout>>(
    () => ({
      autosize: true,
      margin: { l: 56, r: 24, t: 8, b: 36 },
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
    [gene],
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
