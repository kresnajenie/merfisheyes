"use client";

import { useEffect, useState } from "react";
import { Spinner, Progress } from "@heroui/react";

import type { StandardizedDataset } from "@/lib/StandardizedDataset";
import type { PanelType } from "@/lib/stores/splitScreenStore";
import { usePanelDatasetStore, usePanelVisualizationStore, usePanelSingleMoleculeStore, usePanelSingleMoleculeVisualizationStore } from "@/lib/hooks/usePanelStores";
import { selectBestClusterColumn } from "@/lib/utils/dataset-utils";

import { ThreeScene } from "./three-scene";
import { VisualizationControls } from "./visualization-controls";
import UMAPPanel from "./umap-panel";
import { SingleMoleculeThreeScene } from "./single-molecule-three-scene";
import { SingleMoleculeControls } from "./single-molecule-controls";
import { SingleMoleculeLegends } from "./single-molecule-legends";

interface SplitPanelViewerProps {
  datasetId: string | null;
  s3Url: string | null;
  type: PanelType;
}

export function SplitPanelViewer({ datasetId, s3Url, type }: SplitPanelViewerProps) {
  if (type === "cell") {
    return <CellViewer datasetId={datasetId} s3Url={s3Url} />;
  }

  return <SingleMoleculeViewer datasetId={datasetId} s3Url={s3Url} />;
}

function CellViewer({ datasetId, s3Url }: { datasetId: string | null; s3Url: string | null }) {
  const { addDataset } = usePanelDatasetStore();
  const { setSelectedColumn } = usePanelVisualizationStore();
  const [dataset, setDataset] = useState<StandardizedDataset | null>(null);
  const [isLoading, setIsLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [progress, setProgress] = useState(0);
  const [message, setMessage] = useState("Initializing...");

  // Use a stable key to track which source to load
  const sourceKey = s3Url || datasetId;

  useEffect(() => {
    if (!sourceKey) return;
    let cancelled = false;

    async function loadDataset() {
      try {
        setIsLoading(true);
        setError(null);
        setProgress(0);

        const { StandardizedDataset } = await import("@/lib/StandardizedDataset");

        let ds: StandardizedDataset;

        if (s3Url) {
          ds = await StandardizedDataset.fromCustomS3(s3Url, (p, msg) => {
            if (!cancelled) {
              setProgress(p);
              setMessage(msg);
            }
          });
        } else {
          ds = await StandardizedDataset.fromS3(datasetId!, (p, msg) => {
            if (!cancelled) {
              setProgress(p);
              setMessage(msg);
            }
          });
        }

        if (!cancelled) {
          setDataset(ds);
          addDataset(ds);

          const bestColumn = selectBestClusterColumn(ds);

          setSelectedColumn(bestColumn);
          setIsLoading(false);
        }
      } catch (err) {
        if (!cancelled) {
          setError(err instanceof Error ? err.message : "Failed to load dataset");
          setIsLoading(false);
        }
      }
    }

    loadDataset();

    return () => {
      cancelled = true;
    };
  }, [sourceKey]);

  if (isLoading) {
    return (
      <div className="absolute inset-0 flex items-center justify-center bg-black">
        <div className="flex flex-col items-center gap-4 w-full max-w-xs px-4">
          <Spinner color="primary" size="lg" />
          <p className="text-sm text-white/70">Loading dataset...</p>
          <Progress
            aria-label="Loading progress"
            className="w-full"
            color="primary"
            size="sm"
            value={progress}
          />
          <p className="text-xs text-white/50">{message}</p>
        </div>
      </div>
    );
  }

  if (error) {
    return (
      <div className="absolute inset-0 flex items-center justify-center bg-black">
        <div className="text-center px-4">
          <p className="text-sm text-danger">{error}</p>
        </div>
      </div>
    );
  }

  if (!dataset) return null;

  return (
    <>
      <VisualizationControls />
      <ThreeScene dataset={dataset} />
      <UMAPPanel />
    </>
  );
}

function SingleMoleculeViewer({ datasetId, s3Url }: { datasetId: string | null; s3Url: string | null }) {
  const { addDataset } = usePanelSingleMoleculeStore();
  const { addGene } = usePanelSingleMoleculeVisualizationStore();
  const [isLoading, setIsLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [progress, setProgress] = useState(0);
  const [message, setMessage] = useState("Initializing...");
  const [loaded, setLoaded] = useState(false);

  const sourceKey = s3Url || datasetId;

  useEffect(() => {
    if (!sourceKey) return;
    let cancelled = false;

    async function loadDataset() {
      try {
        setIsLoading(true);
        setError(null);
        setProgress(0);

        const { SingleMoleculeDataset } = await import("@/lib/SingleMoleculeDataset");

        let ds;

        if (s3Url) {
          ds = await SingleMoleculeDataset.fromCustomS3(s3Url, (p, msg) => {
            if (!cancelled) {
              setProgress(p);
              setMessage(msg);
            }
          });
        } else {
          ds = await SingleMoleculeDataset.fromS3(datasetId!, (p, msg) => {
            if (!cancelled) {
              setProgress(p);
              setMessage(msg);
            }
          });
        }

        if (!cancelled) {
          addDataset(ds);

          // Auto-select first 3 genes
          const genesToSelect = ds.uniqueGenes.slice(0, 3);

          genesToSelect.forEach((gene) => {
            // Use geneColors if available (from-s3 datasets have them)
            const geneProps = ds.geneColors?.[gene];

            if (geneProps) {
              addGene(gene, geneProps.color, geneProps.size);
            } else {
              addGene(gene);
            }
          });

          setLoaded(true);
          setIsLoading(false);
        }
      } catch (err) {
        if (!cancelled) {
          setError(err instanceof Error ? err.message : "Failed to load dataset");
          setIsLoading(false);
        }
      }
    }

    loadDataset();

    return () => {
      cancelled = true;
    };
  }, [sourceKey]);

  if (isLoading) {
    return (
      <div className="absolute inset-0 flex items-center justify-center bg-black">
        <div className="flex flex-col items-center gap-4 w-full max-w-xs px-4">
          <Spinner color="primary" size="lg" />
          <p className="text-sm text-white/70">Loading dataset...</p>
          <Progress
            aria-label="Loading progress"
            className="w-full"
            color="primary"
            size="sm"
            value={progress}
          />
          <p className="text-xs text-white/50">{message}</p>
        </div>
      </div>
    );
  }

  if (error) {
    return (
      <div className="absolute inset-0 flex items-center justify-center bg-black">
        <div className="text-center px-4">
          <p className="text-sm text-danger">{error}</p>
        </div>
      </div>
    );
  }

  if (!loaded) return null;

  return (
    <>
      <SingleMoleculeControls />
      <SingleMoleculeLegends />
      <SingleMoleculeThreeScene />
    </>
  );
}
