/**
 * Viewer defaults per technology: what a dataset of that technology opens
 * with, before its owner's saved ViewerConfig and before the link (v=), which
 * each override only the fields they carry. Keys are the catalog's platform
 * spellings (docs/specs/conventions.json, `catalog.platforms`).
 *
 * Presets are edited at /presets and read from /api/presets
 * (docs/specs/ADR-technology-presets.md, contract in
 * docs/specs/schemas/preset.schema.json). The table below is what a
 * technology uses until one is published, and what the viewer falls back to
 * when the API cannot be reached.
 */
import type { CellFootprint } from "@/lib/spatialdata/footprint";

import { COLORMAPS } from "@/lib/utils/colormaps";

export interface TechnologyPreset {
  /** Colormap of gene expression and numerical columns. */
  colormap?: string;
  /** Cells with no count of the shown gene are not drawn. */
  hideZeroExpression?: boolean;
  /** The cells' opacity as a whole, 0-1. */
  cellOpacity?: number;
}

export const TECHNOLOGY_PRESETS: Record<string, TechnologyPreset> = {
  "Visium HD": {
    colormap: "viridis",
    hideZeroExpression: true,
    cellOpacity: 0.8,
  },
};

/** The technologies the viewer can tell from the data (see technologyOf). */
export const RECOGNISED_TECHNOLOGIES = ["Visium HD"];

const NONE: TechnologyPreset = {};

/** Published presets by technology; null = turned off. Filled by load. */
let published: Record<string, TechnologyPreset | null> = {};
/** The version named by ?preset=<hash> in the link: used for every dataset
 * opened from that link, whatever is published. */
let preview: TechnologyPreset | null = null;

/** The fields this viewer knows, when their values are usable. */
export function readPreset(raw: unknown): TechnologyPreset {
  const r = (raw ?? {}) as Record<string, unknown>;
  const out: TechnologyPreset = {};

  if (typeof r.colormap === "string" && r.colormap in COLORMAPS)
    out.colormap = r.colormap;
  if (typeof r.hideZeroExpression === "boolean")
    out.hideZeroExpression = r.hideZeroExpression;
  if (
    typeof r.cellOpacity === "number" &&
    r.cellOpacity > 0 &&
    r.cellOpacity <= 1
  )
    out.cellOpacity = r.cellOpacity;

  return out;
}

/**
 * Read the published presets, and the version the link previews if it names
 * one. Called before a dataset is opened. Never throws and never waits long:
 * on any failure the viewer keeps what it had (at first, the table above).
 */
export async function loadTechnologyPresets(): Promise<void> {
  if (typeof window === "undefined") return;
  const get = async (path: string) => {
    const r = await fetch(path, { signal: AbortSignal.timeout(4000) });

    if (!r.ok) throw new Error(String(r.status));

    return r.json();
  };
  const hash = new URLSearchParams(window.location.search).get("preset");

  const [all, one] = await Promise.allSettled([
    get("/api/presets"),
    hash ? get(`/api/presets/versions/${encodeURIComponent(hash)}`) : null,
  ]);

  if (all.status === "fulfilled" && all.value?.presets) {
    published = Object.fromEntries(
      Object.entries(all.value.presets as Record<string, any>).map(
        ([technology, p]) => [technology, p ? readPreset(p.config) : null],
      ),
    );
  }
  preview =
    one.status === "fulfilled" && one.value
      ? readPreset(one.value.config)
      : null;
}

/**
 * The technology the viewer can tell from the data itself: cells that are a
 * grid of square bins are Visium HD. Null for everything else, for now.
 */
export function technologyOf(dataset: {
  cellFootprint?: CellFootprint | null;
}): string | null {
  return dataset.cellFootprint?.shape === "square" ? "Visium HD" : null;
}

/** The preset a dataset opens with; empty when its technology has none. */
export function technologyPresetOf(dataset: {
  cellFootprint?: CellFootprint | null;
}): TechnologyPreset {
  if (preview) return preview;
  const technology = technologyOf(dataset);

  if (!technology) return NONE;
  if (technology in published) return published[technology] ?? NONE;

  return TECHNOLOGY_PRESETS[technology] ?? NONE;
}

/** Set a preset's fields on a visualization store (after its reset). */
export function applyTechnologyPreset(
  preset: TechnologyPreset,
  store: {
    setColormap: (name: string) => void;
    setHideZeroExpression: (hide: boolean) => void;
    setAlphaScale: (alpha: number) => void;
  },
): void {
  if (preset.colormap) store.setColormap(preset.colormap);
  if (preset.hideZeroExpression !== undefined)
    store.setHideZeroExpression(preset.hideZeroExpression);
  if (preset.cellOpacity !== undefined) store.setAlphaScale(preset.cellOpacity);
}
