/**
 * Production-safe E2E instrumentation hooks.
 *
 * Exposes a small `window.__merfish` object that Playwright tests poll to get
 * deterministic timing (time-to-interactive, gene-selection latency) and to
 * read authoritative dataset stats for correctness assertions — instead of
 * relying on `<canvas>` presence plus fixed sleeps.
 *
 * Everything here is a no-op unless `NEXT_PUBLIC_E2E === "1"`, so production
 * builds carry no behavioural change and effectively no cost. The env var is
 * inlined by Next at build time for client bundles, and read at runtime on the
 * server, so the guard is evaluated once and dead-code-eliminated when off.
 *
 * All timestamps use the page's monotonic `performance.now()` clock so deltas
 * computed across separate `markX()` calls are directly comparable.
 */

export interface MerfishRenderStats {
  /** Number of rendered points/cells (single cell). */
  pointCount?: number;
  /** Number of genes available in the dataset. */
  geneCount?: number;
  /** Genes with expression actually loaded (null when unknown). */
  exprGeneCount?: number | null;
  /** Spatial dimensionality (2 or 3). */
  dimensions?: number;
  /** Which pipeline produced this dataset. */
  dataType?: "single_cell" | "single_molecule";
  /** Number of active point clouds (single molecule: one per selected gene). */
  pointClouds?: number;
}

export interface MerfishGeneRender {
  gene: string;
  /** Molecules (single molecule) or cells with expression (single cell). */
  count: number;
  /** performance.now() at the time the gene finished rendering. */
  at: number;
}

/**
 * What one layer (an element kind drawn in the scene) has done so far. The
 * budgets in docs/specs/BUDGETS.md are read off these: time to first content,
 * time to a settled view, what is drawn and what is held.
 */
export interface MerfishLayerStats {
  /** performance.now() when the layer first had content on screen. */
  firstContentAt: number | null;
  /** performance.now() when everything the current view wants was drawn. */
  settledAt: number | null;
  /** performance.now() of the last change to what the view wants. */
  changedAt: number | null;
  /** Molecules / tiles / polygons on screen now. */
  drawn: number;
  /** What the current view asks for, in the same unit as `drawn`. */
  wanted: number;
  /** Kept in memory for this layer, in the same unit as `drawn`. */
  resident: number;
  /** Network reads the layer itself counts, when it does (points element). */
  requests: number;
  bytes: number;
  /** Monotonic; bumps on every settled view. */
  settledCount: number;
}

/**
 * Frame timing since the last `resetFrames()`. Intervals are between ticks
 * of the scene loop (every animation frame, drawn or not), so they measure
 * main-thread jank; an on-demand scene that skips a draw while still is
 * not a slow frame.
 */
export interface MerfishFrameStats {
  /** Loop ticks (animation frames). */
  count: number;
  /** Of those, frames actually drawn. */
  drawn: number;
  /** Ticks whose interval to the previous tick exceeded 20 ms (a dropped 60 Hz frame). */
  over16ms: number;
  /** ... exceeded 34 ms (two dropped frames). */
  over33ms: number;
  /** Longest interval, ms. */
  worstMs: number;
  /** Longest CPU time inside renderer.render, ms. */
  worstRenderMs: number;
  /** Sum of CPU time inside renderer.render, ms. */
  renderMs: number;
  since: number;
}

export interface MerfishMemoryStats {
  /** Chrome only: performance.memory.usedJSHeapSize, bytes. */
  jsHeap: number | null;
  /** three.js renderer.info.memory. */
  geometries: number;
  textures: number;
  /** three.js renderer.info.render of the last frame. */
  drawCalls: number;
  points: number;
  triangles: number;
  at: number | null;
}

/**
 * One molecule gene file's journey through `MoleculeLayer`, keyed by gene
 * (assigned) or `${gene}:unassigned`. Lets a browser test prove a gene was
 * streamed (several batches) rather than loaded whole, or was cancelled.
 */
export interface MerfishGeneLoad {
  /** Distinguishes this load from one that replaced it (scene rebuilt). */
  id: number;
  /** Streamed into the scene as it downloaded, or loaded whole first. */
  mode: "stream" | "whole";
  /** Batches appended so far (1 for a whole load once it lands). */
  batches: number;
  /** Molecules drawn so far. */
  loaded: number;
  /** Molecules the file holds, when known before the download. */
  expected: number | null;
  done: boolean;
  aborted: boolean;
  failed: boolean;
}

export interface MerfishTestHooks {
  enabled: true;
  /** Per element kind: "table", "points", "shapes", "image", "molecules-legacy". */
  layers: Record<string, MerfishLayerStats>;
  frames: MerfishFrameStats;
  memory: MerfishMemoryStats;
  /** Molecule gene loads of the current scene(s); see MerfishGeneLoad. */
  geneLoads: Record<string, MerfishGeneLoad>;
  /** performance.now() when the Three.js scene finished initializing. */
  sceneReady: number | null;
  /** performance.now() of the most recent point-cloud attribute update. */
  renderComplete: number | null;
  /** Monotonic counter incremented on every render-complete. */
  renderCount: number;
  /** performance.now() when the dataset entered the store (pre-render). */
  datasetLoadedAt: number | null;
  /** Authoritative dataset stats, seeded at load and refined at render. */
  stats: MerfishRenderStats | null;
  /** Most recently rendered gene + count. */
  lastGene: MerfishGeneRender | null;
  /** performance.now() when DE stats for the most recent column became ready. */
  deStatsReadyAt: number | null;
  /** Cluster column the most recent DE stats are for. */
  deStatsColumn: string | null;
  /** Monotonic counter incremented every time DE stats resolve. */
  deStatsCount: number;
}

const E2E_ENABLED =
  typeof process !== "undefined" && process.env.NEXT_PUBLIC_E2E === "1";

function nowMs(): number {
  return typeof performance !== "undefined" &&
    typeof performance.now === "function"
    ? performance.now()
    : Date.now();
}

/** Returns the live hooks object, lazily created, or null when disabled. */
function getHooks(): MerfishTestHooks | null {
  if (!E2E_ENABLED || typeof window === "undefined") return null;
  const w = window as unknown as { __merfish?: MerfishTestHooks };

  if (!w.__merfish) {
    w.__merfish = {
      enabled: true,
      layers: {},
      frames: freshFrames(),
      memory: {
        jsHeap: null,
        geometries: 0,
        textures: 0,
        drawCalls: 0,
        points: 0,
        triangles: 0,
        at: null,
      },
      geneLoads: {},
      sceneReady: null,
      renderComplete: null,
      renderCount: 0,
      datasetLoadedAt: null,
      stats: null,
      lastGene: null,
      deStatsReadyAt: null,
      deStatsColumn: null,
      deStatsCount: 0,
    };
  }

  return w.__merfish;
}

function freshFrames(): MerfishFrameStats {
  return {
    count: 0,
    drawn: 0,
    over16ms: 0,
    over33ms: 0,
    worstMs: 0,
    worstRenderMs: 0,
    renderMs: 0,
    since: nowMs(),
  };
}

function freshLayer(): MerfishLayerStats {
  return {
    firstContentAt: null,
    settledAt: null,
    changedAt: null,
    drawn: 0,
    wanted: 0,
    resident: 0,
    requests: 0,
    bytes: 0,
    settledCount: 0,
  };
}

/** Clear all markers — call when a brand-new dataset starts loading. */
export function resetHooks(): void {
  const h = getHooks();

  if (!h) return;
  h.layers = {};
  h.frames = freshFrames();
  h.geneLoads = {};
  h.sceneReady = null;
  h.renderComplete = null;
  h.renderCount = 0;
  h.datasetLoadedAt = null;
  h.stats = null;
  h.lastGene = null;
  h.deStatsReadyAt = null;
  h.deStatsColumn = null;
  h.deStatsCount = 0;
}

/** Record that a dataset has been parsed and added to the store. */
export function markDatasetLoaded(stats?: MerfishRenderStats): void {
  const h = getHooks();

  if (!h) return;
  h.datasetLoadedAt = nowMs();
  if (stats) h.stats = { ...(h.stats ?? {}), ...stats };
}

/** Record that the Three.js scene has been initialized and is animating. */
export function markSceneReady(): void {
  const h = getHooks();

  if (!h) return;
  h.sceneReady = nowMs();
}

/** Record that point-cloud attributes have been (re)applied to the GPU. */
export function markRenderComplete(stats?: MerfishRenderStats): void {
  const h = getHooks();

  if (!h) return;
  h.renderComplete = nowMs();
  h.renderCount += 1;
  if (stats) h.stats = { ...(h.stats ?? {}), ...stats };
  // The cells are a layer like the others: content when first drawn.
  if (stats?.dataType === "single_cell" && stats.pointCount) {
    markLayerContent("table", {
      drawn: stats.pointCount,
      wanted: stats.pointCount,
    });
    markLayerSettled("table");
  }
}

/** Record that a specific gene finished rendering. */
export function markGeneRendered(gene: string, count: number): void {
  const h = getHooks();

  if (!h) return;
  h.lastGene = { gene, count, at: nowMs() };
  markLayerContent("molecules-legacy");
}

// ── layers, frames, memory ───────────────────────────────────────────────

/** The stats slot of a layer kind, created on first use. */
export function layerStats(kind: string): MerfishLayerStats | null {
  const h = getHooks();

  if (!h) return null;

  return (h.layers[kind] ??= freshLayer());
}

/** Update a layer's counters (drawn, wanted, resident, bytes, requests). */
export function markLayer(
  kind: string,
  patch: Partial<
    Omit<MerfishLayerStats, "firstContentAt" | "settledAt" | "settledCount">
  >,
): void {
  const l = layerStats(kind);

  if (!l) return;
  Object.assign(l, patch);
}

/** The view wants something new (a pan, a zoom, a gene): the clock restarts. */
export function markLayerChanged(kind: string, wanted?: number): void {
  const l = layerStats(kind);

  if (!l) return;
  l.changedAt = nowMs();
  l.settledAt = null;
  if (wanted !== undefined) l.wanted = wanted;
}

/** The layer has content on screen (first time only is recorded). */
export function markLayerContent(
  kind: string,
  patch?: Partial<Pick<MerfishLayerStats, "drawn" | "wanted" | "resident">>,
): void {
  const l = layerStats(kind);

  if (!l) return;
  if (patch) Object.assign(l, patch);
  if (l.firstContentAt === null) l.firstContentAt = nowMs();
}

/** Everything the current view wants is drawn. */
export function markLayerSettled(kind: string): void {
  const l = layerStats(kind);

  if (!l) return;
  if (l.settledAt !== null) return;
  l.settledAt = nowMs();
  l.settledCount += 1;
}

/** Start a fresh frame window (a test calls this before a scripted pan). */
export function resetFrames(): void {
  const h = getHooks();

  if (!h) return;
  h.frames = freshFrames();
}

/**
 * One tick of the scene loop: `intervalMs` since the previous tick,
 * `renderMs` spent inside renderer.render (null when nothing was drawn).
 * Called by the scene loop; a no-op in production builds like everything
 * here.
 */
export function markFrame(intervalMs: number, renderMs: number | null): void {
  const h = getHooks();

  if (!h) return;
  const f = h.frames;

  f.count += 1;
  if (intervalMs > 20) f.over16ms += 1;
  if (intervalMs > 34) f.over33ms += 1;
  if (intervalMs > f.worstMs) f.worstMs = intervalMs;
  if (renderMs !== null) {
    f.drawn += 1;
    if (renderMs > f.worstRenderMs) f.worstRenderMs = renderMs;
    f.renderMs += renderMs;
  }
}

/** Sample memory and draw counts from a renderer (about once a second). */
export function markMemory(info: {
  geometries: number;
  textures: number;
  drawCalls: number;
  points: number;
  triangles: number;
}): void {
  const h = getHooks();

  if (!h) return;
  const perf = performance as Performance & {
    memory?: { usedJSHeapSize: number };
  };

  h.memory = {
    jsHeap: perf.memory?.usedJSHeapSize ?? null,
    ...info,
    at: nowMs(),
  };
}

/** True when the hooks are live (NEXT_PUBLIC_E2E=1), for callers that pay to compute a stat. */
export function hooksEnabled(): boolean {
  return getHooks() !== null;
}

/** Key of a gene load in `geneLoads`. */
export function geneLoadKey(
  gene: string,
  kind: "assigned" | "unassigned",
): string {
  return kind === "assigned" ? gene : `${gene}:unassigned`;
}

let nextGeneLoadId = 1;

/**
 * Record that a molecule gene file started loading, and how. Returns the id
 * the load passes back with its progress and end, so a load that was
 * replaced (its scene rebuilt, a new layer started the same gene) can no
 * longer write over the entry of the load that replaced it.
 */
export function markGeneLoadStart(
  key: string,
  mode: MerfishGeneLoad["mode"],
  expected: number | null,
): number {
  const h = getHooks();
  if (!h) return 0;
  const id = nextGeneLoadId++;
  h.geneLoads[key] = {
    id,
    mode,
    batches: 0,
    loaded: 0,
    expected,
    done: false,
    aborted: false,
    failed: false,
  };
  return id;
}

/** Record molecules drawn so far for a gene load (one call per batch). */
export function markGeneLoadProgress(
  key: string,
  id: number,
  loaded: number,
): void {
  const h = getHooks();
  const load = h?.geneLoads[key];
  if (!load || load.id !== id) return;
  load.batches += 1;
  load.loaded = loaded;
}

/** Record how a gene load ended. */
export function markGeneLoadEnd(
  key: string,
  id: number,
  outcome: "done" | "aborted" | "failed",
): void {
  const h = getHooks();
  const load = h?.geneLoads[key];
  if (!load || load.id !== id) return;
  load.done = true;
  load.aborted = outcome === "aborted";
  load.failed = outcome === "failed";
}

/** Record that DE stats for a cluster column have finished computing/loading. */
export function markDeStatsReady(column: string): void {
  const h = getHooks();

  if (!h) return;
  h.deStatsReadyAt = nowMs();
  h.deStatsColumn = column;
  h.deStatsCount += 1;
}
