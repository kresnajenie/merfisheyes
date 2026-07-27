import type { StandardizedDataset } from "@/lib/StandardizedDataset";
import type { VisualizationState } from "@/lib/stores/createVisualizationStore";
import type { SampleTransform } from "@/lib/utils/sample-transforms";

import { selectBestClusterColumn } from "@/lib/utils/dataset-utils";

/**
 * Owner-saved viewer defaults for a dataset. Stored as the `viewer_config`
 * JSON column on `Dataset` and applied on load (see applyViewerConfig).
 *
 * This is the params-only, reversible representation — it overrides the
 * dataset's stored palette/coords at render time rather than requiring the
 * S3 files to be rewritten (though owned-bucket datasets may also bake it in).
 */
export interface ViewerConfig {
  version: 1;
  sceneRotation?: number; // degrees
  flipX?: boolean;
  flipY?: boolean;
  priorityColumn?: string | null; // default cluster column on load
  transformColumn?: string | null; // sample column for per-sample align
  sampleTransforms?: Record<string, SampleTransform>; // keyed by sample VALUE
  colorOverrides?: Record<string, Record<string, string>>; // column -> (value -> hex)
}

const HEX_RE = /^#([0-9a-fA-F]{3}|[0-9a-fA-F]{6})$/;

/**
 * Server-safe. Sanitize arbitrary input into a ViewerConfig, dropping unknown
 * keys and coercing/clamping values. Returns null if nothing valid remains.
 */
export function validateViewerConfig(input: unknown): ViewerConfig | null {
  if (!input || typeof input !== "object") return null;
  const src = input as Record<string, unknown>;
  const out: ViewerConfig = { version: 1 };

  if (typeof src.sceneRotation === "number" && Number.isFinite(src.sceneRotation)) {
    // Normalize to [0, 360)
    out.sceneRotation = ((src.sceneRotation % 360) + 360) % 360;
  }
  if (typeof src.flipX === "boolean") out.flipX = src.flipX;
  if (typeof src.flipY === "boolean") out.flipY = src.flipY;
  if (src.priorityColumn === null || typeof src.priorityColumn === "string") {
    out.priorityColumn = src.priorityColumn as string | null;
  }
  if (src.transformColumn === null || typeof src.transformColumn === "string") {
    out.transformColumn = src.transformColumn as string | null;
  }

  if (src.sampleTransforms && typeof src.sampleTransforms === "object") {
    const st: Record<string, SampleTransform> = {};

    for (const [k, v] of Object.entries(src.sampleTransforms as object)) {
      const t = v as Record<string, unknown>;

      if (
        t &&
        typeof t.dx === "number" &&
        typeof t.dy === "number" &&
        typeof t.theta === "number" &&
        Number.isFinite(t.dx) &&
        Number.isFinite(t.dy) &&
        Number.isFinite(t.theta)
      ) {
        st[k] = { dx: t.dx, dy: t.dy, theta: t.theta };
      }
    }
    if (Object.keys(st).length > 0) out.sampleTransforms = st;
  }

  if (src.colorOverrides && typeof src.colorOverrides === "object") {
    const co: Record<string, Record<string, string>> = {};

    for (const [col, palette] of Object.entries(src.colorOverrides as object)) {
      if (!palette || typeof palette !== "object") continue;
      const clean: Record<string, string> = {};

      for (const [value, hex] of Object.entries(palette as object)) {
        if (typeof hex === "string" && HEX_RE.test(hex)) clean[value] = hex;
      }
      if (Object.keys(clean).length > 0) co[col] = clean;
    }
    if (Object.keys(co).length > 0) out.colorOverrides = co;
  }

  // "version" is the only guaranteed key; treat a config with nothing else as empty.
  return Object.keys(out).length > 1 ? out : null;
}

/**
 * Server-safe. Canonicalize an S3 base URL so the same dataset dedups to one
 * row. Used identically in register-s3, by-url, and dedup lookups.
 */
export function normalizeS3Url(raw: string): string {
  let s = (raw ?? "").trim();

  try {
    const u = new URL(s);

    u.hash = "";
    u.search = "";
    u.hostname = u.hostname.toLowerCase();
    // Strip trailing slashes from the path.
    u.pathname = u.pathname.replace(/\/+$/, "");
    s = u.toString();
  } catch {
    // Not a parseable URL — fall back to a plain trailing-slash strip.
    s = s.replace(/[?#].*$/, "").replace(/\/+$/, "");
  }

  return s;
}

function paletteStorageKey(datasetId: string, column: string): string {
  return `cluster_palette_${datasetId}_${column}`;
}

/**
 * Client-only. Apply an owner-saved config to a visualization store (works for
 * both the live store and the split-panel factory store — same setter surface).
 * Mirrors applyCellVizState in lib/hooks/useUrlVizSync.ts.
 */
export function applyViewerConfig(
  config: ViewerConfig | null | undefined,
  store: VisualizationState,
  dataset: StandardizedDataset,
): void {
  if (!config) return;

  // Colors first, so the palette is in place before the column selection
  // triggers a re-render.
  if (config.colorOverrides) {
    for (const [column, palette] of Object.entries(config.colorOverrides)) {
      const cluster = dataset.clusters?.find((c) => c.column === column);

      if (cluster) {
        cluster.palette = { ...(cluster.palette ?? {}), ...palette };
      }
      // Persist to the same localStorage key the legends read, so the
      // palette-sync effect keeps these colors instead of the S3 defaults.
      if (typeof window !== "undefined") {
        try {
          const key = paletteStorageKey(dataset.id, column);
          const existing = window.localStorage.getItem(key);
          const merged = { ...(existing ? JSON.parse(existing) : {}), ...palette };

          window.localStorage.setItem(key, JSON.stringify(merged));
        } catch {
          // localStorage may be unavailable / full — non-fatal.
        }
      }
    }
  }

  if (config.priorityColumn) {
    const allColumnNames =
      dataset.allClusterColumnNames && dataset.allClusterColumnNames.length > 0
        ? dataset.allClusterColumnNames
        : (dataset.clusters?.map((c) => c.column) ?? []);
    const column = allColumnNames.includes(config.priorityColumn)
      ? config.priorityColumn
      : selectBestClusterColumn(dataset);
    const isNumerical =
      dataset.allClusterColumnTypes?.[column ?? ""] === "numerical";

    store.setSelectedColumn(column, isNumerical);

    if (column && config.colorOverrides?.[column]) {
      const cluster = dataset.clusters?.find((c) => c.column === column);

      store.setColorPalette({ ...(cluster?.palette ?? {}) });
    }
  }

  if (config.sceneRotation !== undefined) store.setSceneRotation(config.sceneRotation);
  if (config.flipX !== undefined) store.setFlipX(config.flipX);
  if (config.flipY !== undefined) store.setFlipY(config.flipY);

  if (config.transformColumn !== undefined) {
    store.setTransformColumn(config.transformColumn);
  }
  if (config.sampleTransforms) {
    for (const [sampleId, t] of Object.entries(config.sampleTransforms)) {
      store.setSampleTransform(sampleId, t);
    }
  }
}

/**
 * Client-only. Read the current store + dataset into a ViewerConfig for saving.
 * Captures rotation/flips, the active cluster column as the priority column,
 * per-sample transforms, and per-column custom colors (from localStorage).
 */
export function extractViewerConfig(
  store: VisualizationState,
  dataset: StandardizedDataset,
): ViewerConfig {
  const config: ViewerConfig = {
    version: 1,
    sceneRotation: store.sceneRotation,
    flipX: store.flipX,
    flipY: store.flipY,
    priorityColumn: store.selectedColumn ?? null,
  };

  if (store.transformColumn) config.transformColumn = store.transformColumn;

  if (store.sampleTransforms && store.sampleTransforms.size > 0) {
    config.sampleTransforms = Object.fromEntries(store.sampleTransforms);
  }

  // Collect any per-column color overrides the user has set (legends persist
  // them to localStorage under cluster_palette_<datasetId>_<column>).
  if (typeof window !== "undefined") {
    const colorOverrides: Record<string, Record<string, string>> = {};
    const prefix = `cluster_palette_${dataset.id}_`;

    try {
      for (let i = 0; i < window.localStorage.length; i++) {
        const key = window.localStorage.key(i);

        if (!key || !key.startsWith(prefix)) continue;
        const column = key.slice(prefix.length);
        const raw = window.localStorage.getItem(key);

        if (!raw) continue;
        const palette = JSON.parse(raw) as Record<string, string>;

        if (palette && typeof palette === "object" && Object.keys(palette).length) {
          colorOverrides[column] = palette;
        }
      }
    } catch {
      // ignore malformed entries
    }
    if (Object.keys(colorOverrides).length > 0) {
      config.colorOverrides = colorOverrides;
    }
  }

  return config;
}
