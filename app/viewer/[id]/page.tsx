"use client";

import { Suspense, useEffect, useState } from "react";
import { useParams, useRouter } from "next/navigation";
import { ThreeScene } from "@/components/three-scene";
import { VisualizationControls } from "@/components/visualization-controls";
import { useVisualizationStore } from "@/lib/stores/visualizationStore";
import { useDatasetStore } from "@/lib/stores/datasetStore";
import { selectBestClusterColumn } from "@/lib/utils/dataset-utils";
import type { StandardizedDataset } from "@/lib/StandardizedDataset";
import LightRays from "@/components/react-bits/LightRays";
import { subtitle, title } from "@/components/primitives";
import { Button, Spinner } from "@heroui/react";

function ViewerByIdContent() {
  const params = useParams();
  const router = useRouter();
  const { setSelectedColumn } = useVisualizationStore();
  const { addDataset } = useDatasetStore();
  const [dataset, setDataset] = useState<StandardizedDataset | null>(null);
  const [isLoading, setIsLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

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

      console.log("Loading dataset from server:", id);

      // Import adapter dynamically
      const { ChunkedDataAdapter } = await import(
        "@/lib/adapters/ChunkedDataAdapter"
      );
      const { StandardizedDataset } = await import(
        "@/lib/StandardizedDataset"
      );

      // Create adapter with dataset ID
      const adapter = new ChunkedDataAdapter(id);

      // Initialize adapter (fetches URLs and loads manifest/index)
      await adapter.initialize();
      console.log("Adapter initialized");

      // Load spatial coordinates
      const spatial = await adapter.loadSpatialCoordinates();
      console.log("Loaded spatial coordinates:", spatial.coordinates.length);

      // Load embeddings
      const embeddings = await adapter.loadEmbeddings();
      console.log("Loaded embeddings:", Object.keys(embeddings));

      // Load genes
      const genes = await adapter.loadGenes();
      console.log("Loaded genes:", genes.length);

      // Load clusters
      const clusters = await adapter.loadClusters();
      console.log("Loaded clusters:", clusters);

      // Get dataset info
      const dataInfo = adapter.getDatasetInfo();
      console.log("Dataset info:", dataInfo);

      // Create StandardizedDataset
      const standardizedDataset = new StandardizedDataset({
        id: dataInfo.id,
        name: dataInfo.name,
        type: dataInfo.type,
        spatial: {
          coordinates: spatial.coordinates,
          dimensions: spatial.dimensions,
        },
        embeddings: embeddings,
        genes: genes,
        clusters: clusters,
        metadata: {
          numCells: dataInfo.numCells,
          numGenes: dataInfo.numGenes,
          spatialDimensions: dataInfo.spatialDimensions,
          availableEmbeddings: dataInfo.availableEmbeddings,
          clusterCount: dataInfo.clusterCount,
        },
        adapter: adapter,
      });

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
        <div className="relative z-10 flex items-center justify-center h-full">
          <div className="flex flex-col items-center gap-4">
            <Spinner size="lg" color="primary" />
            <p className={subtitle()}>Loading dataset...</p>
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
            raysOrigin="top-center"
            raysColor="#FF72E1"
            rayLength={10}
            raysSpeed={0.8}
            lightSpread={1.0}
            pulsating={false}
            mouseInfluence={0.1}
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
                Go to Viewer
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

      {/* Back button in top right */}
      <div className="fixed top-4 right-4 z-50">
        <Button
          color="default"
          variant="bordered"
          onPress={() => router.push("/viewer")}
          size="sm"
        >
          Back to Viewer
        </Button>
      </div>
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
