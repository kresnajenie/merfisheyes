"use client";

import { Suspense, useEffect, useState } from "react";
import { useSearchParams } from "next/navigation";
import { ThreeScene } from "@/components/three-scene";
import { FileUpload } from "@/components/file-upload";
import { VisualizationControls } from "@/components/visualization-controls";
import { useDatasetStore } from "@/lib/stores/datasetStore";
import { useVisualizationStore } from "@/lib/stores/visualizationStore";
import { selectBestClusterColumn } from "@/lib/utils/dataset-utils";
import type { StandardizedDataset } from "@/lib/StandardizedDataset";
import LightRays from "@/components/react-bits/LightRays";
import { subtitle, title } from "@/components/primitives";

function ViewerContent() {
  const searchParams = useSearchParams();
  const { datasets, currentDatasetId, getCurrentDataset } = useDatasetStore();
  const { setSelectedColumn } = useVisualizationStore();
  const [dataset, setDataset] = useState<StandardizedDataset | null>(null);

  useEffect(() => {
    // Check for dataset ID in URL params
    const urlDatasetId = searchParams.get("dataset");

    if (urlDatasetId) {
      // Try to get dataset from URL param
      const datasetFromUrl = datasets.get(urlDatasetId);
      if (datasetFromUrl) {
        console.log("Loading dataset from URL param:", urlDatasetId);
        setDataset(datasetFromUrl);
      } else {
        console.warn(`Dataset with ID ${urlDatasetId} not found in store`);
        // Fall back to current dataset
        const current = getCurrentDataset();
        setDataset(current);
      }
    } else {
      // Use current dataset from store
      const current = getCurrentDataset();
      console.log("Loading current dataset:", current?.id);
      setDataset(current);
    }
  }, [searchParams, datasets, currentDatasetId, getCurrentDataset]);

  // Auto-select best cluster column when dataset changes
  useEffect(() => {
    const bestColumn = selectBestClusterColumn(dataset);
    setSelectedColumn(bestColumn);
    console.log("Auto-selected column:", bestColumn);
    console.log("Dataset:", dataset);
  }, [dataset, setSelectedColumn]);

  if (!dataset) {
    return (
      <>
        <div className="fixed inset-0 w-full h-full z-0">
          <LightRays
            raysOrigin="top-left"
            raysColor="#FFD700"
            rayLength={10}
            raysSpeed={0.8}
            lightSpread={1.0}
            pulsating={false}
            mouseInfluence={0.1}
          />
        </div>
        <div className="fixed inset-0 w-full h-full z-0">
          <LightRays
            raysOrigin="top-right"
            raysColor="#FFD700"
            rayLength={10}
            raysSpeed={0.8}
            lightSpread={1.0}
            pulsating={false}
            mouseInfluence={0.1}
          />
        </div>
        <div className="relative z-10 flex items-center justify-center h-full p-8">
          <div className="flex flex-col items-center gap-6 max-w-5xl w-full">
            <div className="text-center">
              <h2 className={title({ size: "md", color: "yellow" })}>
                No dataset loaded
              </h2>
              <p className={subtitle({ class: "mt-2" })}>
                Upload a dataset to start visualizing
              </p>
            </div>
            <div className="grid grid-cols-3 gap-4 w-full">
              <FileUpload
                type="h5ad"
                title="H5AD File"
                description="Single .h5ad file"
              />
              <FileUpload
                type="xenium"
                title="Xenium Folder"
                description="Select Xenium output folder"
              />
              <FileUpload
                type="merscope"
                title="MERSCOPE Folder"
                description="Select MERSCOPE output folder"
              />
            </div>
          </div>
        </div>
      </>
    );
  }

  return (
    <>
      <VisualizationControls />
      <ThreeScene dataset={dataset} />
    </>
  );
}

export default function ViewerPage() {
  return (
    <Suspense fallback={<div>Loading...</div>}>
      <ViewerContent />
    </Suspense>
  );
}
