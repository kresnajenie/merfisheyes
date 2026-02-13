"use client";

import type { SingleMoleculeDataset } from "@/lib/SingleMoleculeDataset";

import { Suspense, useEffect, useState } from "react";
import { useParams, useRouter, useSearchParams } from "next/navigation";
import { Button, Spinner } from "@heroui/react";

import { SingleMoleculeThreeScene } from "@/components/single-molecule-three-scene";
import { SingleMoleculeControls } from "@/components/single-molecule-controls";
import { SingleMoleculeLegends } from "@/components/single-molecule-legends";
import { SplitScreenContainer } from "@/components/split-screen-container";
import { useSingleMoleculeStore } from "@/lib/stores/singleMoleculeStore";
import { useSingleMoleculeVisualizationStore } from "@/lib/stores/singleMoleculeVisualizationStore";
import { useSplitScreenStore } from "@/lib/stores/splitScreenStore";
import type { PanelType } from "@/lib/stores/splitScreenStore";
import LightRays from "@/components/react-bits/LightRays";
import { subtitle, title } from "@/components/primitives";

function SingleMoleculeViewerByIdContent() {
  const params = useParams();
  const router = useRouter();
  const searchParams = useSearchParams();
  const { addDataset } = useSingleMoleculeStore();
  const { addGene, loadFromLocalStorage, saveToLocalStorage } =
    useSingleMoleculeVisualizationStore();
  const { isSplitMode, rightPanelDatasetId, rightPanelS3Url, rightPanelType, enableSplit, setRightPanel, setRightPanelS3 } =
    useSplitScreenStore();
  const [dataset, setDataset] = useState<SingleMoleculeDataset | null>(null);
  const [isLoading, setIsLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  const datasetId = params.id as string;

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
  }, []);

  // Write split params to URL when split state changes
  useEffect(() => {
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
      router.replace(`?${newParams.toString()}`, { scroll: false });
    } else if (!isSplitMode) {
      const newParams = new URLSearchParams(searchParams.toString());

      newParams.delete("split");
      newParams.delete("splitS3Url");
      newParams.delete("splitType");
      const paramStr = newParams.toString();

      router.replace(paramStr ? `?${paramStr}` : window.location.pathname, {
        scroll: false,
      });
    }
  }, [isSplitMode, rightPanelDatasetId, rightPanelS3Url, rightPanelType]);
  const selectedGenesLegend = useSingleMoleculeVisualizationStore(
    (state) => state.selectedGenesLegend,
  );

  useEffect(() => {
    if (!datasetId) {
      setError("No dataset ID provided");
      setIsLoading(false);

      return;
    }

    loadDatasetFromServer(datasetId);
  }, [datasetId]);

  // Save visibility state to localStorage whenever it changes
  useEffect(() => {
    if (datasetId && dataset) {
      saveToLocalStorage(datasetId);
    }
  }, [selectedGenesLegend, datasetId, dataset, saveToLocalStorage]);

  const loadDatasetFromServer = async (id: string) => {
    try {
      setIsLoading(true);
      setError(null);

      console.log("Loading single molecule dataset from S3:", id);

      // Import SingleMoleculeDataset
      const { SingleMoleculeDataset } = await import(
        "@/lib/SingleMoleculeDataset"
      );

      // Load dataset using fromS3 method with lazy loading
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

      // Store dataset in both local state and global store
      setDataset(smDataset);
      addDataset(smDataset);
      console.log("Dataset added to singleMoleculeStore");

      // Wait a tick for store to update
      await new Promise((resolve) => setTimeout(resolve, 0));

      // Try to load visibility state from localStorage first
      loadFromLocalStorage(id);

      // Wait a bit to see if anything was loaded
      await new Promise((resolve) => setTimeout(resolve, 10));

      // Validate that loaded genes exist in current dataset
      const { selectedGenesLegend, clearGenes } =
        useSingleMoleculeVisualizationStore.getState();

      const validGenes = Array.from(selectedGenesLegend).filter((gene) =>
        smDataset.uniqueGenes.includes(gene)
      );

      // If loaded genes are invalid, clear them
      if (validGenes.length !== selectedGenesLegend.size && selectedGenesLegend.size > 0) {
        console.warn(
          `Loaded ${selectedGenesLegend.size} genes from localStorage, but only ${validGenes.length} exist in current dataset. Clearing invalid genes.`
        );
        clearGenes();
      }

      // If nothing was loaded or all genes were invalid, auto-select first 3 genes
      const { selectedGenesLegend: currentSelection } =
        useSingleMoleculeVisualizationStore.getState();

      if (currentSelection.size === 0) {
        console.log("No valid saved state, auto-selecting first 3 genes");
        const genesToSelect = smDataset.uniqueGenes.slice(0, 3);

        console.log("Auto-selecting genes:", genesToSelect);

        genesToSelect.forEach((gene) => {
          const geneProps = smDataset.geneColors[gene];

          if (!geneProps) {
            console.error(`Missing geneProps for gene: ${gene}`);
            console.error("This should never happen - gene is in uniqueGenes but not in geneColors");
            return;
          }

          console.log(
            `Adding gene to visualization: ${gene} with color ${geneProps.color}`,
          );
          addGene(gene, geneProps.color, geneProps.size);
        });
      } else {
        console.log(
          "Loaded valid visibility state from localStorage:",
          currentSelection.size,
          "genes",
        );
      }

      setIsLoading(false);
    } catch (err) {
      console.error("Error loading single molecule dataset:", err);
      setError(err instanceof Error ? err.message : "Failed to load dataset");
      setIsLoading(false);
    }
  };

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
