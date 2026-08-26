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

export interface MerfishTestHooks {
  enabled: true;
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
  return typeof performance !== "undefined" && typeof performance.now === "function"
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

/** Clear all markers — call when a brand-new dataset starts loading. */
export function resetHooks(): void {
  const h = getHooks();
  if (!h) return;
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
}

/** Record that a specific gene finished rendering. */
export function markGeneRendered(gene: string, count: number): void {
  const h = getHooks();
  if (!h) return;
  h.lastGene = { gene, count, at: nowMs() };
}

/** Record that DE stats for a cluster column have finished computing/loading. */
export function markDeStatsReady(column: string): void {
  const h = getHooks();
  if (!h) return;
  h.deStatsReadyAt = nowMs();
  h.deStatsColumn = column;
  h.deStatsCount += 1;
}
