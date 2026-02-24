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


function ViewerFromLocalContent() {
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
      const datasetFromUrl = datasets.get(urlDatasetId);

      if (datasetFromUrl && isStandardizedDataset(datasetFromUrl)) {
        setDataset(datasetFromUrl);
      } else {
        const current = getCurrentDataset();

        if (current && isStandardizedDataset(current)) {
          setDataset(current);
        }
      }
    } else {
      const current = getCurrentDataset();

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
    }
  }, [dataset]);

  const isLoading = useDatasetStore((state) => state.isLoading);
  const datasetCount = useDatasetStore((state) => state.datasets.size);

  useEffect(() => {
    if (!dataset && !isLoading && datasetCount === 0) {
      router.replace("/");
    }
  }, [dataset, datasetCount, isLoading, router]);

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

export default function ViewerFromLocalPage() {
  return (
    <Suspense fallback={<div>Loading...</div>}>
      <ViewerFromLocalContent />
    </Suspense>
  );
}
