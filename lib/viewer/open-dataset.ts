/**
 * The ONE way a dataset gets opened into a visualization store — used by the
 * full-page viewers AND the split-screen panels, so a dataset always loads
 * identically no matter which surface opened it.
 *
 * Fresh-open semantics, in order:
 *   1. `store.reset()` — clear every residue of a previously shown dataset
 *      (a panel/page reuses its store across dataset switches; without the
 *      reset, the old rotation/column/selection bleeds into the new dataset
 *      and the URL sync then writes that stale state into the share link).
 *   2. Owner-saved ViewerConfig — colors, rotation/flips, priority column,
 *      per-sample transforms (with the transform column lazy-loaded).
 *      Heuristic column selection when there is no config.
 *   3. URL viz state (v= / rv=), if present, overlays ON TOP for exactly the
 *      fields it carries. Callers guarantee this ordering by calling
 *      apply*OpenState BEFORE flipping the state that triggers the URL-sync
 *      hook's read (setDataset / setDatasetReady). Per-sample transforms and
 *      color palettes are not URL-encodable, so they always come from config —
 *      this layering is what keeps a shared link aligned instead of stacked.
 */
import type { StandardizedDataset } from "@/lib/StandardizedDataset";
import type { SingleMoleculeDataset } from "@/lib/SingleMoleculeDataset";
import type { VisualizationState } from "@/lib/stores/createVisualizationStore";
import type { SingleMoleculeVisualizationState } from "@/lib/stores/createSingleMoleculeVisualizationStore";

import {
  applyViewerConfig,
  normalizeS3Url,
  type ViewerConfig,
} from "@/lib/utils/viewer-config";
import { loadClusterColumn } from "@/lib/utils/load-cluster-column";
import { selectBestClusterColumn } from "@/lib/utils/dataset-utils";
import { pickDefaultGenes } from "@/lib/utils/auto-select-genes";
import { tryReadSMVizFromUrl } from "@/lib/hooks/useUrlVizSync";

/**
 * Owner-saved viewer defaults for a dataset, wherever it lives: registered
 * S3 URLs resolve through /api/datasets/by-url; app datasets through their
 * meta route. Best-effort — a failure just opens the dataset unconfigured.
 */
export async function fetchViewerConfig(opts: {
  s3Url?: string | null;
  datasetId?: string | null;
  kind: "cell" | "sm";
}): Promise<ViewerConfig | null> {
  try {
    if (opts.s3Url) {
      const r = await fetch(
        `/api/datasets/by-url?url=${encodeURIComponent(normalizeS3Url(opts.s3Url))}`,
      );

      if (!r.ok) return null;
      const j = await r.json();

      return (j?.viewerConfig as ViewerConfig | null) ?? null;
    }
    if (opts.datasetId) {
      const r = await fetch(
        opts.kind === "sm"
          ? `/api/single-molecule/${opts.datasetId}`
          : `/api/datasets/${opts.datasetId}`,
      );

      if (!r.ok) return null;
      const j = await r.json();

      return (j?.viewerConfig as ViewerConfig | null) ?? null;
    }
  } catch {
    // best-effort; the viewer still loads without owner defaults.
  }

  return null;
}

/**
 * Open a single-cell dataset into a visualization store (steps 1–2 above).
 * Call (and AWAIT) BEFORE setDataset/setDatasetReady so any URL viz state
 * overlays after — and so the loading screen only clears once the per-sample
 * transform data is in place. Without the await, the scene paints the raw
 * (unaligned) coordinates first and the slices visibly jump once the
 * transform column arrives.
 */
export async function applyCellOpenState(opts: {
  dataset: StandardizedDataset;
  config: ViewerConfig | null;
  store: VisualizationState;
}): Promise<void> {
  const { dataset, config, store } = opts;

  store.reset();

  if (config) {
    applyViewerConfig(config, store, dataset);

    // Per-sample align keys off the transform column's per-cell values, which
    // are lazy-loaded on S3 datasets. applyViewerConfig only sets the column
    // name, so fetch its data (if absent) and bump clusterVersion — the
    // scene's transform effect then applies the saved sample transforms.
    // Awaited so the first painted frame is already aligned; a fetch failure
    // just opens unaligned rather than blocking the viewer.
    const tc = config.transformColumn;

    if (tc) {
      try {
        const fetched = await loadClusterColumn(dataset, tc);

        if (fetched) store.incrementClusterVersion();
      } catch {
        // best-effort — transforms simply don't apply.
      }
    }
  } else {
    store.setSelectedColumn(selectBestClusterColumn(dataset));
  }
}

/**
 * Open a single-molecule dataset into its visualization store. Camera config
 * (rotation/flips/view mode) applies unconditionally — a URL state overlays
 * afterwards for what it carries. Default genes come from the owner's config
 * when valid, else the heuristic — skipped entirely when a URL state exists
 * (it brings its own gene selection).
 */
export function applySMOpenState(opts: {
  dataset: SingleMoleculeDataset;
  config: ViewerConfig | null;
  store: SingleMoleculeVisualizationState;
  panel: "left" | "right";
}): void {
  const { dataset, config, store, panel } = opts;

  store.reset();

  if (config) {
    if (config.sceneRotation !== undefined) {
      store.setSceneRotation(config.sceneRotation);
    }
    if (config.flipX !== undefined) store.setFlipX(config.flipX);
    if (config.flipY !== undefined) store.setFlipY(config.flipY);
    if (config.viewMode) store.setViewMode(config.viewMode);
  }

  const urlState = tryReadSMVizFromUrl(panel);

  if (!urlState) {
    const ownerGenes = (config?.defaultGenes ?? []).filter((g) =>
      dataset.uniqueGenes.includes(g),
    );
    const genesToSelect =
      ownerGenes.length > 0
        ? ownerGenes
        : pickDefaultGenes(dataset.uniqueGenes);

    genesToSelect.forEach((gene) => {
      store.addGene(gene);
    });
  }
}
