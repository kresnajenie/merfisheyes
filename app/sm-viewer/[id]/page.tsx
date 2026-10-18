"use client";

import { Suspense } from "react";
import { useParams } from "next/navigation";
import { Spinner } from "@heroui/react";

import { SingleMoleculeThreeScene } from "@/components/single-molecule-three-scene";
import { SingleMoleculeControls } from "@/components/single-molecule-controls";
import { SingleMoleculeLegends } from "@/components/single-molecule-legends";
import { SplitScreenContainer } from "@/components/split-screen-container";
import {
  DatasetErrorScreen,
  DatasetLoadingScreen,
} from "@/components/dataset-loading-screen";
import { useSingleMoleculeVisualizationStore } from "@/lib/stores/singleMoleculeVisualizationStore";
import { useDatasetLoader } from "@/lib/hooks/useDatasetLoader";
import { useSplitUrlParams } from "@/lib/hooks/useSplitUrlParams";
import { useSMVizUrlSync } from "@/lib/hooks/useUrlVizSync";

function SingleMoleculeViewerByIdContent() {
  const datasetId = useParams().id as string;
  const smVizStore = useSingleMoleculeVisualizationStore();
  const { loaded, loading, progress, message, error, retry } = useDatasetLoader(
    { kind: "sm", datasetId },
  );
  const dataset = loaded?.molecules ?? null;

  useSplitUrlParams();
  // URL visualization state sync
  useSMVizUrlSync(!!dataset, dataset, smVizStore);

  if (loading) {
    return (
      <DatasetLoadingScreen
        label="Loading single molecule dataset..."
        message={message}
        progress={progress}
        tone="sm"
      />
    );
  }
  if (error) {
    return (
      <DatasetErrorScreen
        error={error}
        homeHref="/sm-viewer"
        tone="sm"
        onRetry={retry}
      />
    );
  }
  if (!dataset) return null;

  return (
    <SplitScreenContainer>
      <SingleMoleculeControls />
      <SingleMoleculeLegends />
      <SingleMoleculeThreeScene />
    </SplitScreenContainer>
  );
}

export default function SingleMoleculeViewerByIdPage() {
  return (
    <Suspense
      fallback={
        <div className="flex items-center justify-center h-full">
          <Spinner size="lg" />
        </div>
      }
    >
      <SingleMoleculeViewerByIdContent />
    </Suspense>
  );
}
