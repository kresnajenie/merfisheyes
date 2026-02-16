"use client";

import type { StandardizedDataset } from "@/lib/StandardizedDataset";

import { Suspense, useEffect, useState } from "react";
import { useRouter, useSearchParams } from "next/navigation";

import { ThreeScene } from "@/components/three-scene";
import { VisualizationControls } from "@/components/visualization-controls";
import UMAPPanel from "@/components/umap-panel";
import { SplitScreenContainer } from "@/components/split-screen-container";
import { useDatasetStore } from "@/lib/stores/datasetStore";
import { useVisualizationStore } from "@/lib/stores/visualizationStore";
import { selectBestClusterColumn } from "@/lib/utils/dataset-utils";
import { useCellVizUrlSync } from "@/lib/hooks/useUrlVizSync";

function ViewerContent() {
  const router = useRouter();
  const searchParams = useSearchParams();
  const { datasets, currentDatasetId, getCurrentDataset } = useDatasetStore();
  const vizStore = useVisualizationStore();
  const [dataset, setDataset] = useState<StandardizedDataset | null>(null);

  // URL visualization state sync
  const { hasUrlStateRef } = useCellVizUrlSync(!!dataset, dataset, vizStore);

  useEffect(() => {
    // Helper function to check if dataset is StandardizedDataset
    const isStandardizedDataset = (ds: any): ds is StandardizedDataset => {
      return ds && "spatial" in ds && "embeddings" in ds;
    };

    // Check for dataset ID in URL params
    const urlDatasetId = searchParams.get("dataset");

    if (urlDatasetId) {
      // Try to get dataset from URL param
      const datasetFromUrl = datasets.get(urlDatasetId);

      if (datasetFromUrl && isStandardizedDataset(datasetFromUrl)) {
        console.log("Loading dataset from URL param:", urlDatasetId);
        setDataset(datasetFromUrl);
      } else {
        console.warn(
          `Dataset with ID ${urlDatasetId} not found in store or is not a StandardizedDataset`,
        );
        // Fall back to current dataset
        const current = getCurrentDataset();

        if (current && isStandardizedDataset(current)) {
          setDataset(current);
        }
      }
    } else {
      // Use current dataset from store
      const current = getCurrentDataset();

      console.log("Loading current dataset:", current?.id);
      if (current && isStandardizedDataset(current)) {
        setDataset(current);
      }
    }
  }, [searchParams, datasets, currentDatasetId, getCurrentDataset]);

  // Auto-select best cluster column when dataset changes (skip if URL state was applied)
  useEffect(() => {
    if (!hasUrlStateRef.current) {
      const bestColumn = selectBestClusterColumn(dataset);

      vizStore.setSelectedColumn(bestColumn);
      console.log("Auto-selected column:", bestColumn);
      console.log("Dataset:", dataset);
    }
  }, [dataset]);

  const isLoading = useDatasetStore((state) => state.isLoading);
  const datasetCount = useDatasetStore((state) => state.datasets.size);

  useEffect(() => {
    // If we have a current dataset, redirect to /viewer/{id} so the URL is shareable
    if (currentDatasetId && datasets.has(currentDatasetId)) {
      router.replace(`/viewer/${currentDatasetId}`);

      return;
    }

    if (!dataset && !isLoading && datasetCount === 0) {
      router.replace("/");
    }
  }, [dataset, datasetCount, isLoading, router, currentDatasetId, datasets]);

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

export default function ViewerPage() {
  return (
    <Suspense fallback={<div>Loading...</div>}>
      <ViewerContent />
    </Suspense>
  );
}
