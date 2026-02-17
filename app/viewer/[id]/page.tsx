"use client";

import type { StandardizedDataset } from "@/lib/StandardizedDataset";
import type { PanelType } from "@/lib/stores/splitScreenStore";
import type { LocalDatasetMetadata } from "@/lib/services/localDatasetDB";

import { Suspense, useEffect, useState } from "react";
import { useParams, useRouter, useSearchParams } from "next/navigation";
import { Button, Progress, Spinner } from "@heroui/react";

import { ThreeScene } from "@/components/three-scene";
import { VisualizationControls } from "@/components/visualization-controls";
import UMAPPanel from "@/components/umap-panel";
import { SplitScreenContainer } from "@/components/split-screen-container";
import { LocalDatasetReuploadModal } from "@/components/local-dataset-reupload-modal";
import { useVisualizationStore } from "@/lib/stores/visualizationStore";
import { useDatasetStore } from "@/lib/stores/datasetStore";
import { useSplitScreenStore } from "@/lib/stores/splitScreenStore";
import { selectBestClusterColumn } from "@/lib/utils/dataset-utils";
import { useCellVizUrlSync } from "@/lib/hooks/useUrlVizSync";
import {
  isLocalDatasetId,
  getLocalDatasetMeta,
  saveLocalDatasetMeta,
} from "@/lib/services/localDatasetDB";
import LightRays from "@/components/react-bits/LightRays";
import { subtitle, title } from "@/components/primitives";

function ViewerByIdContent() {
  const params = useParams();
  const router = useRouter();
  const searchParams = useSearchParams();
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
  } = useSplitScreenStore();
  const [dataset, setDataset] = useState<StandardizedDataset | null>(null);
  const [isLoading, setIsLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [loadingProgress, setLoadingProgress] = useState(0);
  const [loadingMessage, setLoadingMessage] = useState("Initializing...");
  const [localMetadata, setLocalMetadata] =
    useState<LocalDatasetMetadata | null>(null);
  const [showReuploadModal, setShowReuploadModal] = useState(false);

  const datasetId = params.id as string;

  // URL visualization state sync
  const { hasUrlStateRef } = useCellVizUrlSync(!!dataset, dataset, vizStore);

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
    }
  }, []);

  // Write split params to URL when split state changes
  useEffect(() => {
    // Preserve v/rv params written by replaceState (not in Next.js searchParams)
    const currentUrl = new URLSearchParams(window.location.search);
    const currentV = currentUrl.get("v");
    const currentRv = currentUrl.get("rv");

    if (isSplitMode && rightPanelType) {
      const newParams = new URLSearchParams(searchParams.toString());

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
      if (currentV) newParams.set("v", currentV);
      if (currentRv) newParams.set("rv", currentRv);
      router.replace(`?${newParams.toString()}`, { scroll: false });
    } else if (!isSplitMode) {
      const newParams = new URLSearchParams(searchParams.toString());

      newParams.delete("split");
      newParams.delete("splitS3Url");
      newParams.delete("splitType");
      newParams.delete("sync");
      if (currentV) newParams.set("v", currentV);
      if (currentRv) newParams.set("rv", currentRv);
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
    if (!datasetId) {
      setError("No dataset ID provided");
      setIsLoading(false);

      return;
    }

    resolveDataset(datasetId);
  }, [datasetId]);

  const resolveDataset = async (id: string) => {
    // Step 1: Check Zustand store first
    const storeDataset = useDatasetStore.getState().datasets.get(id);

    if (storeDataset && "spatial" in storeDataset) {
      console.log("Dataset found in store:", id);
      setDataset(storeDataset as StandardizedDataset);
      setIsLoading(false);

      return;
    }

    // Step 2: Check if this is a local dataset
    if (isLocalDatasetId(id)) {
      const meta = await getLocalDatasetMeta(id);

      if (meta) {
        console.log("Local dataset metadata found in IndexedDB:", meta);
        setLocalMetadata(meta);
        setShowReuploadModal(true);
        setIsLoading(false);

        return;
      }

      // Metadata not found (evicted)
      setError(
        "This local dataset is no longer available. The metadata was evicted after loading newer datasets.",
      );
      setIsLoading(false);

      return;
    }

    // Step 3: Load from S3
    loadDatasetFromS3(id);
  };

  const loadDatasetFromS3 = async (id: string) => {
    try {
      setIsLoading(true);
      setError(null);
      setLoadingProgress(0);
      setLoadingMessage("Initializing...");

      console.log("Loading dataset from server:", id);

      const { StandardizedDataset } = await import("@/lib/StandardizedDataset");

      const standardizedDataset = await StandardizedDataset.fromS3(
        id,
        (progress, message) => {
          console.log(`${progress}%: ${message}`);
          setLoadingProgress(progress);
          setLoadingMessage(message);
        },
      );

      console.log("StandardizedDataset created:", standardizedDataset);

      setDataset(standardizedDataset);
      addDataset(standardizedDataset);
      console.log("Dataset added to datasetStore");

      setIsLoading(false);
    } catch (err) {
      console.error("Error loading dataset:", err);
      setError(err instanceof Error ? err.message : "Failed to load dataset");
      setIsLoading(false);
    }
  };

  const handleLocalDatasetLoaded = (ds: any) => {
    const standardizedDataset = ds as StandardizedDataset;

    setDataset(standardizedDataset);
    addDataset(standardizedDataset);
    setShowReuploadModal(false);

    // Refresh timestamp in IndexedDB
    if (localMetadata) {
      saveLocalDatasetMeta({ ...localMetadata, createdAt: Date.now() });
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

  // Re-upload modal for local datasets
  if (showReuploadModal && localMetadata) {
    return (
      <LocalDatasetReuploadModal
        expectedDatasetId={datasetId}
        isOpen={showReuploadModal}
        metadata={localMetadata}
        onClose={() => router.push("/")}
        onDatasetLoaded={handleLocalDatasetLoaded}
      />
    );
  }

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
            raysColor="#FFD700"
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
            raysColor="#FFD700"
            raysOrigin="top-right"
            raysSpeed={0.8}
          />
        </div>
        <div className="relative z-10 flex items-center justify-center h-full">
          <div className="flex flex-col items-center gap-4 w-full max-w-md px-4">
            <Spinner color="primary" size="lg" />
            <p className={subtitle()}>Loading dataset...</p>
            <Progress
              aria-label="Loading progress"
              className="w-full"
              color="primary"
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
                Failed to load dataset
              </h2>
              <p className={subtitle({ class: "mt-4" })}>{error}</p>
            </div>
            <div className="flex gap-4">
              <Button color="primary" onPress={() => router.push("/viewer")}>
                Go to Home
              </Button>
              <Button
                color="default"
                variant="bordered"
                onPress={() => resolveDataset(datasetId)}
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
