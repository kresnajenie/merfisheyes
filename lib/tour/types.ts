export interface TourStopOverrides {
  flyInMs?: number;
  orbitMs?: number;
  orbitRevolutions?: number;
  zoomOutMs?: number;
  zoomRadius?: number;
  // Camera radius (microns) during the zoom-out rest position between stops.
  // Independent of zoomRadius (fly-in radius). Default 1500µm.
  zoomOutRadius?: number;
  // Distance-from-target filter radius (in microns) applied during the
  // orbit phase only. Default 30µm.
  filterRadius?: number;
  // Brief hold between fly-in and orbit during which every gene
  // (persistent + per-stop) becomes visible. Default 150ms.
  waitMs?: number;
  // Brief hold between orbit end and zoom-out start. Per-stop genes are
  // already deselected when this wait begins. Default 150ms.
  waitAfterMs?: number;
}

export interface TourStop extends TourStopOverrides {
  name: string;
  x: number;
  y: number;
  z: number;
  genes: string[];
}

// Genes loaded at every stop in addition to the per-stop `genes` array.
// Useful for background channels (e.g. somatic, total_cell_RNA) that should
// remain visible throughout the tour. Optional `color` is applied on Play
// only — user edits during the tour are preserved between stops.
export interface AlwaysOnGene {
  name: string;
  color?: string;
}

export interface TourConfig {
  defaults?: TourStopOverrides;
  alwaysOnGenes?: AlwaysOnGene[];
  stops: TourStop[];
}

export const TOUR_DEFAULTS: Required<TourStopOverrides> = {
  flyInMs: 1000,
  orbitMs: 1000,
  orbitRevolutions: 1,
  zoomOutMs: 1000,
  zoomRadius: 300,
  zoomOutRadius: 1500,
  filterRadius: 30,
  waitMs: 150,
  waitAfterMs: 150,
};

export function resolveStopTiming(
  stop: TourStop,
  config: TourConfig,
): Required<TourStopOverrides> {
  return {
    flyInMs: stop.flyInMs ?? config.defaults?.flyInMs ?? TOUR_DEFAULTS.flyInMs,
    orbitMs: stop.orbitMs ?? config.defaults?.orbitMs ?? TOUR_DEFAULTS.orbitMs,
    orbitRevolutions:
      stop.orbitRevolutions ??
      config.defaults?.orbitRevolutions ??
      TOUR_DEFAULTS.orbitRevolutions,
    zoomOutMs:
      stop.zoomOutMs ?? config.defaults?.zoomOutMs ?? TOUR_DEFAULTS.zoomOutMs,
    zoomRadius:
      stop.zoomRadius ??
      config.defaults?.zoomRadius ??
      TOUR_DEFAULTS.zoomRadius,
    zoomOutRadius:
      stop.zoomOutRadius ??
      config.defaults?.zoomOutRadius ??
      TOUR_DEFAULTS.zoomOutRadius,
    filterRadius:
      stop.filterRadius ??
      config.defaults?.filterRadius ??
      TOUR_DEFAULTS.filterRadius,
    waitMs: stop.waitMs ?? config.defaults?.waitMs ?? TOUR_DEFAULTS.waitMs,
    waitAfterMs:
      stop.waitAfterMs ??
      config.defaults?.waitAfterMs ??
      TOUR_DEFAULTS.waitAfterMs,
  };
}

export function isValidTourConfig(value: unknown): value is TourConfig {
  if (!value || typeof value !== "object") return false;
  const cfg = value as Partial<TourConfig>;

  if (!Array.isArray(cfg.stops) || cfg.stops.length === 0) return false;

  if (cfg.alwaysOnGenes !== undefined) {
    if (!Array.isArray(cfg.alwaysOnGenes)) return false;
    const allValid = cfg.alwaysOnGenes.every(
      (g) =>
        g &&
        typeof g === "object" &&
        typeof (g as AlwaysOnGene).name === "string" &&
        ((g as AlwaysOnGene).color === undefined ||
          typeof (g as AlwaysOnGene).color === "string"),
    );

    if (!allValid) return false;
  }

  return cfg.stops.every(
    (s) =>
      s &&
      typeof s.name === "string" &&
      typeof s.x === "number" &&
      typeof s.y === "number" &&
      typeof s.z === "number" &&
      Array.isArray(s.genes) &&
      s.genes.every((g) => typeof g === "string"),
  );
}
