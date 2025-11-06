"use client";

import type { StandardizedDataset } from "@/lib/StandardizedDataset";

import { Suspense, useEffect } from "react";
import { Spinner } from "@heroui/react";
import { useRouter } from "next/navigation";

import { ThreeScene } from "@/components/three-scene";
import { VisualizationControls } from "@/components/visualization-controls";
import { subtitle } from "@/components/primitives";
import { selectBestClusterColumn } from "@/lib/utils/dataset-utils";
import { useDatasetStore } from "@/lib/stores/datasetStore";
import { useVisualizationStore } from "@/lib/stores/visualizationStore";
import  LightRays  from "@/components/react-bits/LightRays";

function ViewerContent() {
  const router = useRouter();
  const dataset = useDatasetStore((state) => state.getCurrentDataset());
  const isLoading = useDatasetStore((state) => state.isLoading);
  const { setSelectedColumn } = useVisualizationStore();

  useEffect(() => {
    if (isLoading) return;

    if (!dataset || !("spatial" in dataset) || !("embeddings" in dataset)) {
      router.replace("/");
    }
  }, [dataset, isLoading, router]);

  useEffect(() => {
    if (!dataset || !("clusters" in dataset)) return;

    const standardized = dataset as StandardizedDataset;
    const bestColumn = selectBestClusterColumn(standardized);

    if (bestColumn) {
      const cluster = standardized.clusters?.find(
        (c) => c.column === bestColumn,
      );
      const isNumerical = cluster?.type === "numerical";

      setSelectedColumn(bestColumn, isNumerical);
    }
  }, [dataset, setSelectedColumn]);

  if (isLoading) {
    return (
      <>
        <div className="fixed inset-0 w-full h-full z-0">
          <LightRays
            lightSpread={1.0}
            mouseInfluence={0.1}
            pulsating={false}
            rayLength={10}
            raysColor="#5EA2EF"
            raysOrigin="top-center"
            raysSpeed={1.0}
          />
        </div>
        <div className="relative z-10 flex flex-col items-center justify-center h-full gap-4">
          <Spinner color="primary" size="lg" />
          <p className={subtitle()}>Preparing your dataset viewer...</p>
        </div>
      </>
    );
  }

  if (!dataset || !("spatial" in dataset) || !("embeddings" in dataset)) {
    return (
      <div className="relative z-10 flex flex-col items-center justify-center h-full gap-4">
        <Spinner color="primary" size="lg" />
        <p className={subtitle()}>Redirecting to uploads…</p>
      </div>
    );
  }

  const standardized = dataset as StandardizedDataset;

  return (
    <>
      <VisualizationControls />
      <ThreeScene dataset={standardized} />
    </>
  );
}

export default function ViewerPage() {
  return (
    <Suspense
      fallback={
        <div className="flex items-center justify-center h-full">
          <Spinner size="lg" />
        </div>
      }
    >
      <ViewerContent />
    </Suspense>
  );
}
