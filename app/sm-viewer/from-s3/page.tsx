"use client";

import { Suspense } from "react";
import { useSearchParams } from "next/navigation";
import { Spinner } from "@heroui/react";

import { SingleMoleculeThreeScene } from "@/components/single-molecule-three-scene";
import { SingleMoleculeControls } from "@/components/single-molecule-controls";
import { SingleMoleculeLegends } from "@/components/single-molecule-legends";
import { SplitScreenContainer } from "@/components/split-screen-container";
import { ClaimDatasetBanner } from "@/components/claim-dataset-banner";
import {
  DatasetErrorScreen,
  DatasetLoadingScreen,
} from "@/components/dataset-loading-screen";
import { useSingleMoleculeVisualizationStore } from "@/lib/stores/singleMoleculeVisualizationStore";
import { useDatasetLoader } from "@/lib/hooks/useDatasetLoader";
import { useSplitUrlParams } from "@/lib/hooks/useSplitUrlParams";
import { useSMVizUrlSync } from "@/lib/hooks/useUrlVizSync";

function SingleMoleculeViewerFromS3Content() {
  const s3Url = useSearchParams().get("url");
  const smVizStore = useSingleMoleculeVisualizationStore();
  const { loaded, loading, progress, message, error, retry } = useDatasetLoader(
    {
      kind: "sm",
      s3Url: s3Url ? decodeURIComponent(s3Url) : null,
    },
  );
  const dataset = loaded?.molecules ?? null;

  useSplitUrlParams();
  // URL visualization state sync
  useSMVizUrlSync(!!dataset, dataset, smVizStore);

  if (loading) {
    return (
      <DatasetLoadingScreen
        label="Loading dataset from S3..."
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
        heading="Failed to load dataset from S3"
        homeHref="/sm-viewer"
        tone="sm"
        onRetry={retry}
      />
    );
  }
  if (!dataset) return null;

  return (
    <SplitScreenContainer>
      <ClaimDatasetBanner />
      <SingleMoleculeControls />
      <SingleMoleculeLegends />
      <SingleMoleculeThreeScene />
    </SplitScreenContainer>
  );
}

export default function SingleMoleculeViewerFromS3Page() {
  return (
    <Suspense
      fallback={
        <div className="flex items-center justify-center h-full">
          <Spinner size="lg" />
        </div>
      }
    >
      <SingleMoleculeViewerFromS3Content />
    </Suspense>
  );
}
