/**
 * Matplotlib-style colormaps for gene/numerical visualization.
 * Each cmap is a list of [stop, hex] entries; we linearly interpolate
 * between adjacent stops in sRGB. Stops sampled at ~5–7 points per cmap
 * for visual fidelity without pulling in d3-scale-chromatic.
 *
 * Convention: value=0 → first stop, value=1 → last stop.
 */

export type ColormapName =
  | "bwr"
  | "coolwarm"
  | "viridis"
  | "plasma"
  | "inferno"
  | "magma"
  | "cividis"
  | "turbo"
  | "Reds"
  | "Blues"
  | "Greens"
  | "Greys";

type ColormapDef = {
  name: ColormapName;
  label: string;
  stops: ReadonlyArray<readonly [number, string]>;
};

export const COLORMAPS: Record<ColormapName, ColormapDef> = {
  bwr: {
    name: "bwr",
    label: "bwr",
    stops: [
      [0.0, "#0000ff"],
      [0.5, "#ffffff"],
      [1.0, "#ff0000"],
    ],
  },
  coolwarm: {
    name: "coolwarm",
    label: "coolwarm",
    stops: [
      [0.0, "#3b4cc0"],
      [0.25, "#7ea7f1"],
      [0.5, "#dddddd"],
      [0.75, "#f59585"],
      [1.0, "#b40426"],
    ],
  },
  viridis: {
    name: "viridis",
    label: "viridis",
    stops: [
      [0.0, "#440154"],
      [0.25, "#3b528b"],
      [0.5, "#21918c"],
      [0.75, "#5ec962"],
      [1.0, "#fde725"],
    ],
  },
  plasma: {
    name: "plasma",
    label: "plasma",
    stops: [
      [0.0, "#0d0887"],
      [0.25, "#7e03a8"],
      [0.5, "#cc4778"],
      [0.75, "#f89540"],
      [1.0, "#f0f921"],
    ],
  },
  inferno: {
    name: "inferno",
    label: "inferno",
    stops: [
      [0.0, "#000004"],
      [0.25, "#420a68"],
      [0.5, "#932667"],
      [0.75, "#dd513a"],
      [1.0, "#fcffa4"],
    ],
  },
  magma: {
    name: "magma",
    label: "magma",
    stops: [
      [0.0, "#000004"],
      [0.25, "#3b0f70"],
      [0.5, "#8c2981"],
      [0.75, "#e85b6e"],
      [1.0, "#fcfdbf"],
    ],
  },
  cividis: {
    name: "cividis",
    label: "cividis",
    stops: [
      [0.0, "#00224e"],
      [0.25, "#2e486d"],
      [0.5, "#696e74"],
      [0.75, "#a4a279"],
      [1.0, "#f6e36b"],
    ],
  },
  turbo: {
    name: "turbo",
    label: "turbo",
    stops: [
      [0.0, "#30123b"],
      [0.166, "#4145ab"],
      [0.333, "#4675ed"],
      [0.5, "#1ae4b6"],
      [0.666, "#82e234"],
      [0.833, "#f3a83b"],
      [1.0, "#7a0403"],
    ],
  },
  Reds: {
    name: "Reds",
    label: "Reds",
    stops: [
      [0.0, "#fff5f0"],
      [0.5, "#fb6a4a"],
      [1.0, "#67000d"],
    ],
  },
  Blues: {
    name: "Blues",
    label: "Blues",
    stops: [
      [0.0, "#f7fbff"],
      [0.5, "#6baed6"],
      [1.0, "#08306b"],
    ],
  },
  Greens: {
    name: "Greens",
    label: "Greens",
    stops: [
      [0.0, "#f7fcf5"],
      [0.5, "#74c476"],
      [1.0, "#00441b"],
    ],
  },
  Greys: {
    name: "Greys",
    label: "Greys",
    stops: [
      [0.0, "#ffffff"],
      [0.5, "#969696"],
      [1.0, "#000000"],
    ],
  },
};

export const DIVERGING_COLORMAPS: ColormapName[] = ["bwr", "coolwarm"];

export const SEQUENTIAL_COLORMAPS: ColormapName[] = [
  "viridis",
  "plasma",
  "inferno",
  "magma",
  "cividis",
  "turbo",
  "Reds",
  "Blues",
  "Greens",
  "Greys",
];

export const DEFAULT_COLORMAP: ColormapName = "bwr";

function hexToRgb01(hex: string): [number, number, number] {
  const m = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex);
  if (!m) return [0, 0, 0];
  return [parseInt(m[1], 16) / 255, parseInt(m[2], 16) / 255, parseInt(m[3], 16) / 255];
}

function isValidColormap(name: string): name is ColormapName {
  return name in COLORMAPS;
}

export function getColormap(name: string): ColormapDef {
  return isValidColormap(name) ? COLORMAPS[name] : COLORMAPS[DEFAULT_COLORMAP];
}

/**
 * Returns RGB triplet (0–1) for a value in [0, 1] under the named colormap.
 * NaN → white. Out-of-range values are clamped.
 */
export function colormapRgb(
  name: string,
  value: number,
): [number, number, number] {
  if (Number.isNaN(value)) return [1, 1, 1];
  const v = Math.max(0, Math.min(1, value));
  const cmap = getColormap(name);
  const stops = cmap.stops;

  for (let i = 0; i < stops.length - 1; i++) {
    const [t0, hex0] = stops[i];
    const [t1, hex1] = stops[i + 1];
    if (v <= t1) {
      const t = (v - t0) / (t1 - t0);
      const c0 = hexToRgb01(hex0);
      const c1 = hexToRgb01(hex1);
      return [
        c0[0] + (c1[0] - c0[0]) * t,
        c0[1] + (c1[1] - c0[1]) * t,
        c0[2] + (c1[2] - c0[2]) * t,
      ];
    }
  }
  return hexToRgb01(stops[stops.length - 1][1]);
}

/**
 * CSS linear-gradient string for the named colormap.
 * Direction: "to top" puts value=1 at the top (used by the vertical scalebar);
 * "to right" for horizontal preview swatches.
 */
export function colormapCss(
  name: string,
  direction: "to top" | "to bottom" | "to right" | "to left" = "to top",
): string {
  const cmap = getColormap(name);
  // For "to top": low values at bottom (0%), high at top (100%) — matches stop order.
  // For "to bottom": reverse the stops so high is at top.
  const reversed = direction === "to bottom" || direction === "to left";
  const stops = reversed
    ? cmap.stops.map(([t, hex]) => [1 - t, hex] as const).reverse()
    : cmap.stops;
  const css = stops.map(([t, hex]) => `${hex} ${(t * 100).toFixed(1)}%`).join(", ");
  return `linear-gradient(${direction}, ${css})`;
}
