"use client";

import { Suspense } from "react";
import { useSearchParams } from "next/navigation";
import { Spinner } from "@heroui/react";

import { ThreeScene } from "@/components/three-scene";
import { VisualizationControls } from "@/components/visualization-controls";
import UMAPPanel from "@/components/umap-panel";
import { SplitScreenContainer } from "@/components/split-screen-container";
import { ClaimDatasetBanner } from "@/components/claim-dataset-banner";
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

function ViewerFromS3Content() {
  const s3Url = useSearchParams().get("url");
  const vizStore = useVisualizationStore();
  const { loaded, loading, progress, message, error, retry } = useDatasetLoader(
    {
      kind: "cell",
      s3Url: s3Url ? decodeURIComponent(s3Url) : null,
    },
  );
  const dataset = loaded?.cells ?? null;
  const viewerConfig = useViewerRegistrationStore((s) => s.viewerConfig);

  useSplitUrlParams();
  // URL visualization state sync (SC viz writes to v=)
  useCellVizUrlSync(!!dataset, dataset, vizStore);

  // SM overlay URL state sync — uses a dedicated ov= slot so it doesn't
  // collide with the SC viz. The SM dataset is loaded asynchronously by
  // three-scene.tsx when mapping.json has linkColumn "__all__" and pushed
  // into the global SM store, so we read it from there.
  const smDataset = useSingleMoleculeStore((s) => {
    const id = s.currentDatasetId;

    return id ? (s.datasets.get(id) ?? null) : null;
  });
  const smVizStore = useSingleMoleculeVisualizationStore();

  useSMOverlayUrlSync(
    !!smDataset,
    smDataset,
    smVizStore,
    overlayDefaultGenes(viewerConfig),
  );

  if (loading) {
    return (
      <DatasetLoadingScreen
        label="Loading dataset from S3..."
        message={message}
        progress={progress}
      />
    );
  }
  if (error) {
    return (
      <DatasetErrorScreen
        error={error}
        heading="Failed to load dataset from S3"
        homeHref="/viewer"
        onRetry={retry}
      />
    );
  }
  if (!dataset) return null;

  return (
    <SplitScreenContainer>
      <ClaimDatasetBanner />
      <VisualizationControls />
      <ThreeScene dataset={dataset} />
      <UMAPPanel />
    </SplitScreenContainer>
  );
}

export default function ViewerFromS3Page() {
  return (
    <Suspense
      fallback={
        <div className="flex items-center justify-center h-full">
          <Spinner size="lg" />
        </div>
      }
    >
      <ViewerFromS3Content />
    </Suspense>
  );
}
