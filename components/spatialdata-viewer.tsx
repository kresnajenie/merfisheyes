"use client";

import { ClaimDatasetBanner } from "@/components/claim-dataset-banner";
import {
  DatasetErrorScreen,
  DatasetLoadingScreen,
} from "@/components/dataset-loading-screen";
import { ImageScene } from "@/components/image-viewer";
import { SingleMoleculeControls } from "@/components/single-molecule-controls";
import { SingleMoleculeLegends } from "@/components/single-molecule-legends";
import { SingleMoleculeThreeScene } from "@/components/single-molecule-three-scene";
import { ThreeScene } from "@/components/three-scene";
import { VisualizationControls } from "@/components/visualization-controls";
import { SplitScreenContainer } from "@/components/split-screen-container";
import { useDatasetLoader } from "@/lib/hooks/useDatasetLoader";
import { useSMOverlayUrlSync } from "@/lib/hooks/useUrlVizSync";
import { useSingleMoleculeStore } from "@/lib/stores/singleMoleculeStore";
import { useSingleMoleculeVisualizationStore } from "@/lib/stores/singleMoleculeVisualizationStore";

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
 */
export function SpatialDataViewer({
  datasetId,
  url,
}: {
  datasetId?: string;
  url?: string;
}) {
  const { loaded, loading, progress, message, error, retry } = useDatasetLoader(
    { kind: "spatialdata", datasetId, s3Url: url },
  );

  // The scene loads the molecule dataset into the global SM store; this
  // picks its default genes (or restores them from the URL).
  const smOverlayDataset = useSingleMoleculeStore((s) =>
    s.currentDatasetId ? (s.datasets.get(s.currentDatasetId) ?? null) : null,
  );
  const smVizStore = useSingleMoleculeVisualizationStore();

  useSMOverlayUrlSync(!!smOverlayDataset, smOverlayDataset, smVizStore);

  if (error) {
    return <DatasetErrorScreen error={error} homeHref="/" onRetry={retry} />;
  }
  if (loading || !loaded) {
    return (
      <DatasetLoadingScreen
        label="Loading SpatialData store..."
        message={message}
        progress={progress}
      />
    );
  }
  const { cells, molecules, elements, storeUrl } = loaded;

  if (cells) {
    return (
      <SplitScreenContainer>
        {storeUrl && <ClaimDatasetBanner />}
        <VisualizationControls />
        <ThreeScene dataset={cells} elements={elements} />
      </SplitScreenContainer>
    );
  }
  if (molecules) {
    return (
      <SplitScreenContainer>
        {storeUrl && <ClaimDatasetBanner />}
        <SingleMoleculeControls />
        <SingleMoleculeLegends />
        <SingleMoleculeThreeScene elements={elements} />
      </SplitScreenContainer>
    );
  }

  return (
    elements.image && (
      <ImageScene image={elements.image}>
        {storeUrl && <ClaimDatasetBanner />}
      </ImageScene>
    )
  );
}
