"use client";

import type { SingleMoleculeDataset } from "@/lib/SingleMoleculeDataset";
import type { PanelType } from "@/lib/stores/splitScreenStore";
import type { LocalDatasetMetadata } from "@/lib/services/localDatasetDB";

import { Suspense, useEffect, useState } from "react";
import { useParams, useRouter, useSearchParams } from "next/navigation";
import { Button, Spinner } from "@heroui/react";

import { SingleMoleculeThreeScene } from "@/components/single-molecule-three-scene";
import { SingleMoleculeControls } from "@/components/single-molecule-controls";
import { SingleMoleculeLegends } from "@/components/single-molecule-legends";
import { SplitScreenContainer } from "@/components/split-screen-container";
import { LocalDatasetReuploadModal } from "@/components/local-dataset-reupload-modal";
import { useSingleMoleculeStore } from "@/lib/stores/singleMoleculeStore";
import { pickDefaultGenes } from "@/lib/utils/auto-select-genes";
import { useSingleMoleculeVisualizationStore } from "@/lib/stores/singleMoleculeVisualizationStore";
import { useSplitScreenStore } from "@/lib/stores/splitScreenStore";
import {
  useSMVizUrlSync,
  tryReadSMVizFromUrl,
  applySMVizState,
} from "@/lib/hooks/useUrlVizSync";
import {
  isLocalDatasetId,
  getLocalDatasetMeta,
  saveLocalDatasetMeta,
} from "@/lib/services/localDatasetDB";
import LightRays from "@/components/react-bits/LightRays";
import { subtitle, title } from "@/components/primitives";

function SingleMoleculeViewerByIdContent() {
  const params = useParams();
  const router = useRouter();
  const searchParams = useSearchParams();
  const { addDataset } = useSingleMoleculeStore();
  const smVizStore = useSingleMoleculeVisualizationStore();
  const { addGene, loadFromLocalStorage, saveToLocalStorage } = smVizStore;
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
  const [dataset, setDataset] = useState<SingleMoleculeDataset | null>(null);
  const [isLoading, setIsLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [localMetadata, setLocalMetadata] =
    useState<LocalDatasetMetadata | null>(null);
  const [showReuploadModal, setShowReuploadModal] = useState(false);
  const datasetId = params.id as string;

  // URL visualization state sync
  useSMVizUrlSync(!!dataset, dataset, smVizStore);

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
  const selectedGenesLegend = useSingleMoleculeVisualizationStore(
    (state) => state.selectedGenesLegend,
  );

  useEffect(() => {
    if (!datasetId) {
      setError("No dataset ID provided");
      setIsLoading(false);

      return;
    }

    resolveDataset(datasetId);
  }, [datasetId]);

  // Save visibility state to localStorage whenever it changes
  useEffect(() => {
    if (datasetId && dataset) {
      saveToLocalStorage(datasetId);
    }
  }, [selectedGenesLegend, datasetId, dataset, saveToLocalStorage]);

  const resolveDataset = async (id: string) => {
    // Step 1: Check Zustand store first
    const storeDataset = useSingleMoleculeStore.getState().datasets.get(id);

    if (storeDataset && "uniqueGenes" in storeDataset) {
      console.log("SM dataset found in store:", id);
      setDataset(storeDataset as SingleMoleculeDataset);
      applyVizStateForDataset(storeDataset as SingleMoleculeDataset, id);
      setIsLoading(false);

      return;
    }

    // Step 2: Check if this is a local dataset
    if (isLocalDatasetId(id)) {
      const meta = await getLocalDatasetMeta(id);

      if (meta) {
        console.log("Local SM dataset metadata found in IndexedDB:", meta);
        setLocalMetadata(meta);
        setShowReuploadModal(true);
        setIsLoading(false);

        return;
      }

      setError(
        "This local dataset is no longer available. The metadata was evicted after loading newer datasets.",
      );
      setIsLoading(false);

      return;
    }

    // Step 3: Load from S3
    loadDatasetFromS3(id);
  };

  const applyVizStateForDataset = async (
    smDataset: SingleMoleculeDataset,
    id: string,
  ) => {
    const urlVizState = tryReadSMVizFromUrl("left");

    if (urlVizState) {
      console.log("Applying visualization state from URL");
      applySMVizState(urlVizState, smVizStore, smDataset);
    } else {
      await new Promise((resolve) => setTimeout(resolve, 0));
      loadFromLocalStorage(id);
      await new Promise((resolve) => setTimeout(resolve, 10));

      const { selectedGenesLegend, clearGenes } =
        useSingleMoleculeVisualizationStore.getState();
      const validGenes = Array.from(selectedGenesLegend).filter((gene) =>
        smDataset.uniqueGenes.includes(gene),
      );

      if (
        validGenes.length !== selectedGenesLegend.size &&
        selectedGenesLegend.size > 0
      ) {
        clearGenes();
      }

      const { selectedGenesLegend: currentSelection } =
        useSingleMoleculeVisualizationStore.getState();

      if (currentSelection.size === 0) {
        const genesToSelect = pickDefaultGenes(smDataset.uniqueGenes);

        genesToSelect.forEach((gene) => {
          addGene(gene);
        });
      }
    }
  };

  const loadDatasetFromS3 = async (id: string) => {
    try {
      setIsLoading(true);
      setError(null);

      console.log("Loading single molecule dataset from S3:", id);

      const { SingleMoleculeDataset } = await import(
        "@/lib/SingleMoleculeDataset"
      );

      const smDataset = await SingleMoleculeDataset.fromS3(
        id,
        (progress, message) => {
          console.log(`${progress}%: ${message}`);
        },
      );

      console.log(
        "SingleMoleculeDataset loaded from S3:",
        smDataset.getSummary(),
      );

      setDataset(smDataset);
      addDataset(smDataset);
      console.log("Dataset added to singleMoleculeStore");

      await applyVizStateForDataset(smDataset, id);

      setIsLoading(false);
    } catch (err) {
      console.error("Error loading single molecule dataset:", err);
      setError(err instanceof Error ? err.message : "Failed to load dataset");
      setIsLoading(false);
    }
  };

  const handleLocalDatasetLoaded = (ds: any) => {
    const smDataset = ds as SingleMoleculeDataset;

    setDataset(smDataset);
    addDataset(smDataset);
    setShowReuploadModal(false);

    // Refresh timestamp in IndexedDB
    if (localMetadata) {
      saveLocalDatasetMeta({ ...localMetadata, createdAt: Date.now() });
    }

    // Apply viz state
    applyVizStateForDataset(smDataset, datasetId);
  };

  // Re-upload modal for local datasets
  if (showReuploadModal && localMetadata) {
    return (
      <LocalDatasetReuploadModal
        expectedDatasetId={datasetId}
        isOpen={showReuploadModal}
        metadata={localMetadata}
        onClose={() => router.push("/?mode=sm")}
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
          <div className="flex flex-col items-center gap-4">
            <Spinner color="secondary" size="lg" />
            <p className={subtitle()}>Loading single molecule dataset...</p>
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
              <Button
                color="secondary"
                onPress={() => router.push("/sm-viewer")}
              >
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
      <SingleMoleculeControls />
      <SingleMoleculeLegends />
      <SingleMoleculeThreeScene />
    </SplitScreenContainer>
  );
}

export default function SingleMoleculeViewerByIdPage() {
  return (
    <Suspense
      fallback={
        <div className="flex items-center justify-center h-full">
          <Spinner size="lg" />
        </div>
      }
    >
      <SingleMoleculeViewerByIdContent />
    </Suspense>
  );
}
