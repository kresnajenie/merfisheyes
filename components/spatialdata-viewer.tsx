"use client";

import type { ReactNode } from "react";

import { ClaimDatasetBanner } from "@/components/claim-dataset-banner";
import {
  DatasetErrorScreen,
  DatasetLoadingScreen,
} from "@/components/dataset-loading-screen";
import { ImageScene } from "@/components/image-viewer";
import { ImageViewerControls } from "@/components/image-viewer-controls";
import { SingleMoleculeControls } from "@/components/single-molecule-controls";
import { SingleMoleculeLegends } from "@/components/single-molecule-legends";
import { SingleMoleculeThreeScene } from "@/components/single-molecule-three-scene";
import { ThreeScene } from "@/components/three-scene";
import UMAPPanel from "@/components/umap-panel";
import { VisualizationControls } from "@/components/visualization-controls";
import { SplitScreenPage } from "@/components/split-screen-container";
import { useDatasetLoader } from "@/lib/hooks/useDatasetLoader";
import {
  usePanelSingleMoleculeStore,
  usePanelSingleMoleculeVisualizationStore,
  usePanelVisualizationStore,
} from "@/lib/hooks/usePanelStores";
import {
  useCellVizUrlSync,
  useSMOverlayUrlSync,
} from "@/lib/hooks/useUrlVizSync";

/**
 * Viewer for a SpatialData zarr store: by dataset id (uploaded to our S3) or
 * by public URL. A store may hold cells (a table), molecules (a points
 * element, or a linked single molecule dataset named by the mapping.json
 * beside it), cell polygons (a shapes element) and an image, in any
 * combination:
 *   - with cells: the single cell scene, with the rest drawn in it
 *   - molecules but no cells: the single molecule scene, image underneath
 *   - an image alone: just the image
 * The store is opened by the one dataset loader; which elements it holds
 * is decided by the element registry (lib/scene/registry.tsx).
 *
 * `panel`: the right side of a split view (components/split-panel-viewer.tsx),
 * in that panel's stores, without the page's split container, claim banner or
 * a counted view.
 */
export function SpatialDataViewer({
  datasetId,
  url,
  panel = false,
}: {
  datasetId?: string;
  url?: string;
  panel?: boolean;
}) {
  const { loaded, loading, progress, message, error, retry } = useDatasetLoader(
    { kind: "spatialdata", datasetId, s3Url: url },
    { countViews: !panel },
  );

  // The scene loads the molecule dataset into the panel's SM store; this
  // picks its default genes (or restores them from the URL).
  const smOverlayDataset = usePanelSingleMoleculeStore((s) =>
    s.currentDatasetId ? (s.datasets.get(s.currentDatasetId) ?? null) : null,
  );
  const smVizStore = usePanelSingleMoleculeVisualizationStore();

  useSMOverlayUrlSync(!!smOverlayDataset, smOverlayDataset, smVizStore);

  // The cells' view in the link (v=, rv= in the right split panel), as the
  // cell viewer keeps it: 2D/3D, gene, column, cell types, so a shared link
  // opens as it was left
  const vizStore = usePanelVisualizationStore();
  const cellsForUrl = loaded?.cells ?? null;

  useCellVizUrlSync(!!cellsForUrl, cellsForUrl, vizStore);

  if (error) {
    return panel ? (
      <DatasetErrorScreen error={error} variant="panel" />
    ) : (
      <DatasetErrorScreen error={error} homeHref="/" onRetry={retry} />
    );
  }
  if (loading || !loaded) {
    return (
      <DatasetLoadingScreen
        label="Loading SpatialData store..."
        message={message}
        progress={progress}
        variant={panel ? "panel" : undefined}
      />
    );
  }
  const { cells, molecules, elements, storeUrl } = loaded;
  // The page wraps its scene for the split; a panel is already inside one
  const frame = (scene: ReactNode) =>
    panel ? (
      scene
    ) : (
      <SplitScreenPage>
        {storeUrl && <ClaimDatasetBanner />}
        {scene}
      </SplitScreenPage>
    );

  if (cells) {
    return frame(
      <>
        <VisualizationControls />
        <ThreeScene dataset={cells} elements={elements} />
        {/* A table's embeddings (obsm/X_umap, ...), as the cell viewer
            shows them: a saved store opens here, its preview in /viewer */}
        <UMAPPanel />
      </>,
    );
  }
  if (molecules) {
    return frame(
      <>
        <SingleMoleculeControls />
        <SingleMoleculeLegends />
        <SingleMoleculeThreeScene elements={elements} />
      </>,
    );
  }

  // An image alone: the image scene, with the buttons the other scenes have
  return (
    elements.image &&
    frame(
      <>
        <ImageViewerControls />
        <ImageScene image={elements.image} />
      </>,
    )
  );
}
