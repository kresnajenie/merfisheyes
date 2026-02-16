"use client";

import { Suspense, useEffect } from "react";
import { useRouter } from "next/navigation";

import { SingleMoleculeThreeScene } from "@/components/single-molecule-three-scene";
import { SingleMoleculeControls } from "@/components/single-molecule-controls";
import { SingleMoleculeLegends } from "@/components/single-molecule-legends";
import { SplitScreenContainer } from "@/components/split-screen-container";
import { useSingleMoleculeStore } from "@/lib/stores/singleMoleculeStore";
import { useSingleMoleculeVisualizationStore } from "@/lib/stores/singleMoleculeVisualizationStore";
import { useSMVizUrlSync } from "@/lib/hooks/useUrlVizSync";
import { pickDefaultGenes } from "@/lib/utils/auto-select-genes";

function ViewerContent() {
  const router = useRouter();
  const currentDatasetId = useSingleMoleculeStore(
    (state) => state.currentDatasetId,
  );
  // Get dataset from single molecule store
  const dataset = useSingleMoleculeStore((state) => {
    const id = state.currentDatasetId;

    return id ? state.datasets.get(id) : null;
  });

  const smVizStore = useSingleMoleculeVisualizationStore();
  const { addGene, clearGenes } = smVizStore;

  // URL visualization state sync
  const { hasUrlStateRef } = useSMVizUrlSync(
    !!dataset,
    dataset ?? null,
    smVizStore,
  );

  // Auto-select default genes when dataset loads (skip if URL state was applied)
  useEffect(() => {
    if (!dataset || hasUrlStateRef.current) return;

    clearGenes();

    const genesToSelect = pickDefaultGenes(dataset.uniqueGenes);

    genesToSelect.forEach((gene) => {
      const geneProps = dataset.geneColors[gene];

      if (geneProps) {
        addGene(gene, geneProps.color, geneProps.size);
      }
    });
  }, [dataset, addGene, clearGenes]);

  useEffect(() => {
    // If we have a current dataset, redirect to /sm-viewer/{id} so the URL is shareable
    if (currentDatasetId && dataset) {
      router.replace(`/sm-viewer/${currentDatasetId}`);

      return;
    }

    if (!dataset) {
      router.replace("/?mode=sm");
    }
  }, [dataset, router, currentDatasetId]);

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

export default function ViewerPage() {
  return (
    <Suspense fallback={<div>Loading...</div>}>
      <ViewerContent />
    </Suspense>
  );
}
