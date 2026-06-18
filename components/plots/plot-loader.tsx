"use client";

import dynamic from "next/dynamic";
import type { PlotParams } from "react-plotly.js";

/**
 * Dynamically-loaded Plotly component.
 *
 * Uses plotly.js-cartesian-dist-min (~2.5MB), which includes bar, histogram,
 * box, scatter, contour, heatmap. When violin support is needed, switch to
 * plotly.js-dist-min (~3.5MB) — same import path otherwise.
 *
 * Wrapped in next/dynamic with ssr:false because Plotly touches `window`
 * at module init.
 */
export const Plot = dynamic<PlotParams>(
  async () => {
    const Plotly = (await import("plotly.js-cartesian-dist-min")).default;
    const create = (await import("react-plotly.js/factory")).default;
    return create(Plotly);
  },
  {
    ssr: false,
    loading: () => (
      <div className="w-full h-full flex items-center justify-center text-xs text-default-400">
        Loading plot…
      </div>
    ),
  },
);
