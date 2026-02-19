"use client";

import type { StandardizedDataset } from "@/lib/StandardizedDataset";
import type { PanelType } from "@/lib/stores/splitScreenStore";

import { Suspense, useEffect, useState } from "react";
import { useRouter, useSearchParams } from "next/navigation";
import { Button, Progress, Spinner } from "@heroui/react";

import { ThreeScene } from "@/components/three-scene";
import { VisualizationControls } from "@/components/visualization-controls";
import UMAPPanel from "@/components/umap-panel";
import { SplitScreenContainer } from "@/components/split-screen-container";
import { useVisualizationStore } from "@/lib/stores/visualizationStore";
import { useDatasetStore } from "@/lib/stores/datasetStore";
import { useSplitScreenStore } from "@/lib/stores/splitScreenStore";
import { selectBestClusterColumn } from "@/lib/utils/dataset-utils";
import { useCellVizUrlSync } from "@/lib/hooks/useUrlVizSync";
import { useBackgroundClusterLoader } from "@/lib/hooks/useBackgroundClusterLoader";
import LightRays from "@/components/react-bits/LightRays";
import { subtitle, title } from "@/components/primitives";

function ViewerFromS3Content() {
  const searchParams = useSearchParams();
  const router = useRouter();
  const vizStore = useVisualizationStore();
  const { addDataset } = useDatasetStore();
  const {
    isSplitMode,
    rightPanelDatasetId,
    rightPanelS3Url,
    rightPanelType,
    syncEnabled,
    enableSplit,
    setRightPanel,
    setRightPanelS3,
    setSyncEnabled,
    setSyncFromUrl,
  } = useSplitScreenStore();
  const [dataset, setDataset] = useState<StandardizedDataset | null>(null);
  const [isLoading, setIsLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [loadingProgress, setLoadingProgress] = useState(0);
  const [loadingMessage, setLoadingMessage] = useState("Initializing...");

  const s3Url = searchParams.get("url");

  // URL visualization state sync
  const { hasUrlStateRef } = useCellVizUrlSync(!!dataset, dataset, vizStore);

  // Background-load remaining cluster columns after priority column
  useBackgroundClusterLoader(dataset, vizStore.incrementClusterVersion);

  // Read split params from URL on mount
  useEffect(() => {
    const splitId = searchParams.get("split");
    const splitS3Url = searchParams.get("splitS3Url");
    const splitType = searchParams.get("splitType") as PanelType | null;

    if (splitS3Url && splitType) {
      enableSplit();
      setRightPanelS3(decodeURIComponent(splitS3Url), splitType);
    } else if (splitId && splitType) {
      enableSplit();
      setRightPanel(splitId, splitType);
    }

    if (searchParams.get("sync") === "1") {
      setSyncEnabled(true);
      setSyncFromUrl(true);
    }
  }, []);

  // Write split params to URL when split state changes
  useEffect(() => {
    // Use window.location.search as base to avoid stale Next.js searchParams
    const newParams = new URLSearchParams(window.location.search);

    if (isSplitMode && rightPanelType) {
      if (rightPanelS3Url) {
        newParams.set("splitS3Url", encodeURIComponent(rightPanelS3Url));
        newParams.delete("split");
      } else if (rightPanelDatasetId) {
        newParams.set("split", rightPanelDatasetId);
        newParams.delete("splitS3Url");
      }
      newParams.set("splitType", rightPanelType);
      if (syncEnabled) {
        newParams.set("sync", "1");
      } else {
        newParams.delete("sync");
      }
      router.replace(`?${newParams.toString()}`, { scroll: false });
    } else if (!isSplitMode) {
      newParams.delete("split");
      newParams.delete("splitS3Url");
      newParams.delete("splitType");
      newParams.delete("sync");
      const paramStr = newParams.toString();

      router.replace(paramStr ? `?${paramStr}` : window.location.pathname, {
        scroll: false,
      });
    }
  }, [
    isSplitMode,
    rightPanelDatasetId,
    rightPanelS3Url,
    rightPanelType,
    syncEnabled,
  ]);

  useEffect(() => {
    if (!s3Url) {
      setError("No S3 URL provided");
      setIsLoading(false);

      return;
    }

    loadDatasetFromCustomS3(decodeURIComponent(s3Url));
  }, [s3Url]);

  const loadDatasetFromCustomS3 = async (baseUrl: string) => {
    try {
      setIsLoading(true);
      setError(null);
      setLoadingProgress(0);
      setLoadingMessage("Initializing...");

      console.log("Loading dataset from custom S3:", baseUrl);

      // Import StandardizedDataset
      const { StandardizedDataset } = await import("@/lib/StandardizedDataset");

      // Load dataset using fromCustomS3 method
      const standardizedDataset = await StandardizedDataset.fromCustomS3(
        baseUrl,
        (progress, message) => {
          console.log(`${progress}%: ${message}`);
          setLoadingProgress(progress);
          setLoadingMessage(message);
        },
      );

      console.log(
        "StandardizedDataset created from custom S3:",
        standardizedDataset,
      );

      // Store dataset in both local state and global store
      setDataset(standardizedDataset);
      addDataset(standardizedDataset);
      console.log("Dataset added to datasetStore");

      setIsLoading(false);
    } catch (err) {
      console.error("Error loading dataset from custom S3:", err);
      setError(
        err instanceof Error
          ? err.message
          : "Failed to load dataset from custom S3",
      );
      setIsLoading(false);
    }
  };

  // Auto-select best cluster column when dataset changes (skip if URL state was applied)
  useEffect(() => {
    if (dataset && !hasUrlStateRef.current) {
      const bestColumn = selectBestClusterColumn(dataset);

      vizStore.setSelectedColumn(bestColumn);
      console.log("Auto-selected column:", bestColumn);
    }
  }, [dataset]);

  // Loading state
  if (isLoading) {
    return (
      <>
        <div className="fixed inset-0 w-full h-full z-0">
          <LightRays
            lightSpread={1.0}
            mouseInfluence={0.1}
            pulsating={false}
            rayLength={10}
            raysColor="#667eea"
            raysOrigin="top-left"
            raysSpeed={0.8}
          />
        </div>
        <div className="fixed inset-0 w-full h-full z-0">
          <LightRays
            lightSpread={1.0}
            mouseInfluence={0.1}
            pulsating={false}
            rayLength={10}
            raysColor="#764ba2"
            raysOrigin="top-right"
            raysSpeed={0.8}
          />
        </div>
        <div className="relative z-10 flex items-center justify-center h-full">
          <div className="flex flex-col items-center gap-4 w-full max-w-md px-4">
            <Spinner color="secondary" size="lg" />
            <p className={subtitle()}>Loading dataset from S3...</p>
            <Progress
              aria-label="Loading progress"
              className="w-full"
              color="secondary"
              size="md"
              value={loadingProgress}
            />
            <p className="text-sm text-default-500">{loadingMessage}</p>
          </div>
        </div>
      </>
    );
  }

  // Error state
  if (error) {
    return (
      <>
        <div className="fixed inset-0 w-full h-full z-0">
          <LightRays
            lightSpread={1.0}
            mouseInfluence={0.1}
            pulsating={false}
            rayLength={10}
            raysColor="#FF72E1"
            raysOrigin="top-center"
            raysSpeed={0.8}
          />
        </div>
        <div className="relative z-10 flex items-center justify-center h-full p-8">
          <div className="flex flex-col items-center gap-6 max-w-2xl w-full">
            <div className="text-center">
              <h2 className={title({ size: "md", color: "pink" })}>
                Failed to load dataset from S3
              </h2>
              <p className={subtitle({ class: "mt-4" })}>{error}</p>
            </div>
            <div className="flex gap-4">
              <Button color="secondary" onPress={() => router.push("/viewer")}>
                Go to Home
              </Button>
              <Button
                color="default"
                variant="bordered"
                onPress={() =>
                  s3Url && loadDatasetFromCustomS3(decodeURIComponent(s3Url))
                }
              >
                Retry
              </Button>
            </div>
          </div>
        </div>
      </>
    );
  }

  // Dataset loaded
  if (!dataset) {
    return null;
  }

  return (
    <SplitScreenContainer>
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
