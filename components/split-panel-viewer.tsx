"use client";

import type { StandardizedDataset } from "@/lib/StandardizedDataset";
import type { SingleMoleculeDataset } from "@/lib/SingleMoleculeDataset";
import type { PanelType } from "@/lib/stores/splitScreenStore";
import type { LocalDatasetMetadata } from "@/lib/services/localDatasetDB";

import { Spinner, Progress } from "@heroui/react";
import { useEffect, useState } from "react";

import { LocalDatasetReuploadModal } from "./local-dataset-reupload-modal";
import { ThreeScene } from "./three-scene";
import { VisualizationControls } from "./visualization-controls";
import UMAPPanel from "./umap-panel";
import { SingleMoleculeThreeScene } from "./single-molecule-three-scene";
import { SingleMoleculeControls } from "./single-molecule-controls";
import { SingleMoleculeLegends } from "./single-molecule-legends";

import { pickDefaultGenes } from "@/lib/utils/auto-select-genes";
import {
  isLocalDatasetId,
  getLocalDatasetMeta,
  saveLocalDatasetMeta,
} from "@/lib/services/localDatasetDB";
import {
  tryReadCellVizFromUrl,
  tryReadSMVizFromUrl,
  useCellVizUrlSync,
  useSMVizUrlSync,
} from "@/lib/hooks/useUrlVizSync";
import { selectBestClusterColumn } from "@/lib/utils/dataset-utils";
import { useSingleMoleculeStore } from "@/lib/stores/singleMoleculeStore";
import { useDatasetStore } from "@/lib/stores/datasetStore";
import {
  usePanelDatasetStore,
  usePanelVisualizationStore,
  usePanelSingleMoleculeStore,
  usePanelSingleMoleculeVisualizationStore,
} from "@/lib/hooks/usePanelStores";
import { useBackgroundClusterLoader } from "@/lib/hooks/useBackgroundClusterLoader";

interface SplitPanelViewerProps {
  datasetId: string | null;
  s3Url: string | null;
  type: PanelType;
}

export function SplitPanelViewer({
  datasetId,
  s3Url,
  type,
}: SplitPanelViewerProps) {
  if (type === "cell") {
    return <CellViewer datasetId={datasetId} s3Url={s3Url} />;
  }

  return <SingleMoleculeViewer datasetId={datasetId} s3Url={s3Url} />;
}

function CellViewer({
  datasetId,
  s3Url,
}: {
  datasetId: string | null;
  s3Url: string | null;
}) {
  const { addDataset } = usePanelDatasetStore();
  const vizStore = usePanelVisualizationStore();
  const [dataset, setDataset] = useState<StandardizedDataset | null>(null);
  const [isLoading, setIsLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [progress, setProgress] = useState(0);
  const [message, setMessage] = useState("Initializing...");
  const [datasetReady, setDatasetReady] = useState(false);
  const [localMetadata, setLocalMetadata] =
    useState<LocalDatasetMetadata | null>(null);
  const [showReupload, setShowReupload] = useState(false);

  // URL sync hook (handles reading after datasetReady + writing)
  useCellVizUrlSync(datasetReady, dataset, vizStore);

  // Background-load remaining cluster columns after priority column
  useBackgroundClusterLoader(dataset, vizStore.incrementClusterVersion);

  // Use a stable key to track which source to load
  const sourceKey = s3Url || datasetId;

  useEffect(() => {
    if (!sourceKey) return;
    let cancelled = false;

    async function resolveDataset() {
      // If no S3 URL, check store first
      if (!s3Url && datasetId) {
        const storeDataset = useDatasetStore.getState().datasets.get(datasetId);

        if (storeDataset && "spatial" in storeDataset) {
          const ds = storeDataset as StandardizedDataset;

          setDataset(ds);
          addDataset(ds);

          const urlState = tryReadCellVizFromUrl("right");

          if (!urlState) {
            vizStore.setSelectedColumn(selectBestClusterColumn(ds));
          }

          setDatasetReady(true);
          setIsLoading(false);

          return;
        }

        // Check if local dataset
        if (isLocalDatasetId(datasetId)) {
          const meta = await getLocalDatasetMeta(datasetId);

          if (meta) {
            setLocalMetadata(meta);
            setShowReupload(true);
            setIsLoading(false);

            return;
          }

          setError("Local dataset metadata not found (evicted).");
          setIsLoading(false);

          return;
        }
      }

      // Load from S3
      try {
        setIsLoading(true);
        setError(null);
        setProgress(0);

        const { StandardizedDataset } = await import(
          "@/lib/StandardizedDataset"
        );

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

          const urlState = tryReadCellVizFromUrl("right");

          if (!urlState) {
            vizStore.setSelectedColumn(selectBestClusterColumn(ds));
          }

          setDatasetReady(true);
          setIsLoading(false);
        }
      } catch (err) {
        if (!cancelled) {
          setError(
            err instanceof Error ? err.message : "Failed to load dataset",
          );
          setIsLoading(false);
        }
      }
    }

    resolveDataset();

    return () => {
      cancelled = true;
    };
  }, [sourceKey]);

  const handleLocalDatasetLoaded = (ds: any) => {
    const standardizedDataset = ds as StandardizedDataset;

    setDataset(standardizedDataset);
    addDataset(standardizedDataset);
    setShowReupload(false);

    const urlState = tryReadCellVizFromUrl("right");

    if (!urlState) {
      vizStore.setSelectedColumn(selectBestClusterColumn(standardizedDataset));
    }

    setDatasetReady(true);

    if (localMetadata) {
      saveLocalDatasetMeta({ ...localMetadata, createdAt: Date.now() });
    }
  };

  if (showReupload && localMetadata && datasetId) {
    return (
      <div className="absolute inset-0 flex items-center justify-center bg-black">
        <LocalDatasetReuploadModal
          expectedDatasetId={datasetId}
          isOpen={showReupload}
          metadata={localMetadata}
          onClose={() => setShowReupload(false)}
          onDatasetLoaded={handleLocalDatasetLoaded}
        />
      </div>
    );
  }

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

function SingleMoleculeViewer({
  datasetId,
  s3Url,
}: {
  datasetId: string | null;
  s3Url: string | null;
}) {
  const { addDataset } = usePanelSingleMoleculeStore();
  const smVizStore = usePanelSingleMoleculeVisualizationStore();
  const [smDataset, setSmDataset] = useState<SingleMoleculeDataset | null>(
    null,
  );
  const [isLoading, setIsLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [progress, setProgress] = useState(0);
  const [message, setMessage] = useState("Initializing...");
  const [datasetReady, setDatasetReady] = useState(false);
  const [localMetadata, setLocalMetadata] =
    useState<LocalDatasetMetadata | null>(null);
  const [showReupload, setShowReupload] = useState(false);

  // URL sync hook (handles reading after datasetReady + writing)
  useSMVizUrlSync(datasetReady, smDataset, smVizStore);

  const sourceKey = s3Url || datasetId;

  const autoSelectGenes = (ds: SingleMoleculeDataset) => {
    const urlState = tryReadSMVizFromUrl("right");

    if (!urlState) {
      const genesToSelect = pickDefaultGenes(ds.uniqueGenes);

      genesToSelect.forEach((gene) => {
        smVizStore.addGene(gene);
      });
    }
  };

  useEffect(() => {
    if (!sourceKey) return;
    let cancelled = false;

    async function resolveDataset() {
      // If no S3 URL, check store first
      if (!s3Url && datasetId) {
        const storeDataset = useSingleMoleculeStore
          .getState()
          .datasets.get(datasetId);

        if (storeDataset && "uniqueGenes" in storeDataset) {
          const ds = storeDataset as SingleMoleculeDataset;

          addDataset(ds);
          setSmDataset(ds);
          autoSelectGenes(ds);
          setDatasetReady(true);
          setIsLoading(false);

          return;
        }

        if (isLocalDatasetId(datasetId)) {
          const meta = await getLocalDatasetMeta(datasetId);

          if (meta) {
            setLocalMetadata(meta);
            setShowReupload(true);
            setIsLoading(false);

            return;
          }

          setError("Local dataset metadata not found (evicted).");
          setIsLoading(false);

          return;
        }
      }

      // Load from S3
      try {
        setIsLoading(true);
        setError(null);
        setProgress(0);

        const { SingleMoleculeDataset } = await import(
          "@/lib/SingleMoleculeDataset"
        );

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
          setSmDataset(ds);
          autoSelectGenes(ds);
          setDatasetReady(true);
          setIsLoading(false);
        }
      } catch (err) {
        if (!cancelled) {
          setError(
            err instanceof Error ? err.message : "Failed to load dataset",
          );
          setIsLoading(false);
        }
      }
    }

    resolveDataset();

    return () => {
      cancelled = true;
    };
  }, [sourceKey]);

  const handleLocalDatasetLoaded = (ds: any) => {
    const smDs = ds as SingleMoleculeDataset;

    addDataset(smDs);
    setSmDataset(smDs);
    setShowReupload(false);
    autoSelectGenes(smDs);
    setDatasetReady(true);

    if (localMetadata) {
      saveLocalDatasetMeta({ ...localMetadata, createdAt: Date.now() });
    }
  };

  if (showReupload && localMetadata && datasetId) {
    return (
      <div className="absolute inset-0 flex items-center justify-center bg-black">
        <LocalDatasetReuploadModal
          expectedDatasetId={datasetId}
          isOpen={showReupload}
          metadata={localMetadata}
          onClose={() => setShowReupload(false)}
          onDatasetLoaded={handleLocalDatasetLoaded}
        />
      </div>
    );
  }

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

  if (!datasetReady) return null;

  return (
    <>
      <SingleMoleculeControls />
      <SingleMoleculeLegends />
      <SingleMoleculeThreeScene />
    </>
  );
}
