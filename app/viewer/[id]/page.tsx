"use client";

import { Suspense } from "react";
import { useParams } from "next/navigation";
import { Spinner } from "@heroui/react";

import { ThreeScene } from "@/components/three-scene";
import { VisualizationControls } from "@/components/visualization-controls";
import UMAPPanel from "@/components/umap-panel";
import { SplitScreenContainer } from "@/components/split-screen-container";
import {
  DatasetErrorScreen,
  DatasetLoadingScreen,
} from "@/components/dataset-loading-screen";
import { useVisualizationStore } from "@/lib/stores/visualizationStore";
import { useViewerRegistrationStore } from "@/lib/stores/viewerRegistrationStore";
import { useSingleMoleculeStore } from "@/lib/stores/singleMoleculeStore";
import { useSingleMoleculeVisualizationStore } from "@/lib/stores/singleMoleculeVisualizationStore";
import { useDatasetLoader } from "@/lib/hooks/useDatasetLoader";
import { useSplitUrlParams } from "@/lib/hooks/useSplitUrlParams";
import {
  useCellVizUrlSync,
  useSMOverlayUrlSync,
} from "@/lib/hooks/useUrlVizSync";
import { overlayDefaultGenes } from "@/lib/viewer/overlay-default-genes";

function ViewerByIdContent() {
  const datasetId = useParams().id as string;
  const vizStore = useVisualizationStore();
  const { loaded, loading, progress, message, error, retry } = useDatasetLoader(
    { kind: "cell", datasetId },
  );
  const dataset = loaded?.cells ?? null;
  const viewerConfig = useViewerRegistrationStore((s) => s.viewerConfig);

  useSplitUrlParams();
  // URL visualization state sync
  useCellVizUrlSync(!!dataset, dataset, vizStore);

  // SM overlay URL state sync — uses a dedicated ov= slot so the overlay's
  // selected molecule genes persist in the link. The SM dataset is loaded
  // asynchronously by three-scene.tsx when mapping.json has linkColumn
  // "__all__" and pushed into the global SM store, so we read it from there.
  const smOverlayDataset = useSingleMoleculeStore((s) => {
    const id = s.currentDatasetId;

    return id ? (s.datasets.get(id) ?? null) : null;
  });
  const smVizStore = useSingleMoleculeVisualizationStore();

  useSMOverlayUrlSync(
    !!smOverlayDataset,
    smOverlayDataset,
    smVizStore,
    overlayDefaultGenes(viewerConfig),
  );

  if (loading) {
    return <DatasetLoadingScreen message={message} progress={progress} />;
  }
  if (error) {
    return (
      <DatasetErrorScreen error={error} homeHref="/viewer" onRetry={retry} />
    );
  }
  if (!dataset) return null;

  return (
    <SplitScreenContainer>
      <VisualizationControls />
      <ThreeScene dataset={dataset} />
      <UMAPPanel />
    </SplitScreenContainer>
  );
}

export default function ViewerByIdPage() {
  return (
    <Suspense
      fallback={
        <div className="flex items-center justify-center h-full">
          <Spinner size="lg" />
        </div>
      }
    >
      <ViewerByIdContent />
    </Suspense>
  );
}
