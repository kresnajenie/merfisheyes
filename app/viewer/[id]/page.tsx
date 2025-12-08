"use client";

import type { StandardizedDataset } from "@/lib/StandardizedDataset";

import { Suspense, useEffect, useState } from "react";
import { useParams, useRouter } from "next/navigation";
import { Button, Progress, Spinner } from "@heroui/react";

import { ThreeScene } from "@/components/three-scene";
import { VisualizationControls } from "@/components/visualization-controls";
import UMAPPanel from "@/components/umap-panel";
import { useVisualizationStore } from "@/lib/stores/visualizationStore";
import { useDatasetStore } from "@/lib/stores/datasetStore";
import { selectBestClusterColumn } from "@/lib/utils/dataset-utils";
import LightRays from "@/components/react-bits/LightRays";
import { subtitle, title } from "@/components/primitives";

function ViewerByIdContent() {
  const params = useParams();
  const router = useRouter();
  const { setSelectedColumn } = useVisualizationStore();
  const { addDataset } = useDatasetStore();
  const [dataset, setDataset] = useState<StandardizedDataset | null>(null);
  const [isLoading, setIsLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [loadingProgress, setLoadingProgress] = useState(0);
  const [loadingMessage, setLoadingMessage] = useState("Initializing...");

  const datasetId = params.id as string;

  useEffect(() => {
    if (!datasetId) {
      setError("No dataset ID provided");
      setIsLoading(false);

      return;
    }

    loadDatasetFromServer(datasetId);
  }, [datasetId]);

  const loadDatasetFromServer = async (id: string) => {
    try {
      setIsLoading(true);
      setError(null);
      setLoadingProgress(0);
      setLoadingMessage("Initializing...");

      console.log("Loading dataset from server:", id);

      // Import StandardizedDataset
      const { StandardizedDataset } = await import("@/lib/StandardizedDataset");

      // Load dataset using fromS3 method
      const standardizedDataset = await StandardizedDataset.fromS3(
        id,
        (progress, message) => {
          console.log(`${progress}%: ${message}`);
          setLoadingProgress(progress);
          setLoadingMessage(message);
        },
      );

      console.log("StandardizedDataset created:", standardizedDataset);

      // Store dataset in both local state and global store
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

  // Auto-select best cluster column when dataset changes
  useEffect(() => {
    if (dataset) {
      const bestColumn = selectBestClusterColumn(dataset);

      setSelectedColumn(bestColumn);
      console.log("Auto-selected column:", bestColumn);
    }
  }, [dataset, setSelectedColumn]);

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
                onPress={() => loadDatasetFromServer(datasetId)}
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
    <>
      <VisualizationControls />
      <ThreeScene dataset={dataset} />
      <UMAPPanel />
    </>
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
