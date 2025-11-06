"use client";

import { Suspense, useEffect, useRef } from "react";
import { Spinner } from "@heroui/react";
import { useRouter } from "next/navigation";

import LightRays from "@/components/react-bits/LightRays";
import { subtitle } from "@/components/primitives";
import { SingleMoleculeControls } from "@/components/single-molecule-controls";
import { SingleMoleculeLegends } from "@/components/single-molecule-legends";
import { SingleMoleculeThreeScene } from "@/components/single-molecule-three-scene";
import { useSingleMoleculeStore } from "@/lib/stores/singleMoleculeStore";
import { useSingleMoleculeVisualizationStore } from "@/lib/stores/singleMoleculeVisualizationStore";

function SingleMoleculeViewerContent() {
  const router = useRouter();
  const dataset = useSingleMoleculeStore((state) => state.getCurrentDataset());
  const isLoading = useSingleMoleculeStore((state) => state.isLoading);
  const addGene = useSingleMoleculeVisualizationStore(
    (state) => state.addGene,
  );
  const clearGenes = useSingleMoleculeVisualizationStore(
    (state) => state.clearGenes,
  );
  const lastDatasetIdRef = useRef<string | null>(null);

  useEffect(() => {
    if (!dataset) return;

    if (lastDatasetIdRef.current === dataset.id) return;

    clearGenes();

    const genesToSelect = dataset.uniqueGenes.slice(0, 3);

    genesToSelect.forEach((gene) => {
      const geneProps = dataset.geneColors[gene];
      addGene(gene, geneProps?.color, geneProps?.size);
    });

    lastDatasetIdRef.current = dataset.id;
  }, [dataset, addGene, clearGenes]);

  useEffect(() => {
    if (isLoading) return;

    if (!dataset) {
      router.replace("/?mode=sm");
    }
  }, [dataset, isLoading, router]);

  if (isLoading) {
    return (
      <>
        <div className="fixed inset-0 w-full h-full z-0">
          <LightRays
            lightSpread={1.0}
            mouseInfluence={0.1}
            pulsating={false}
            rayLength={10}
            raysColor="#FF1CF7"
            raysOrigin="top-center"
            raysSpeed={1.0}
          />
        </div>
        <div className="relative z-10 flex flex-col items-center justify-center h-full gap-4">
          <Spinner color="secondary" size="lg" />
          <p className={subtitle()}>Preparing your single molecule viewer...</p>
        </div>
      </>
    );
  }

  if (!dataset) {
    return (
      <div className="relative z-10 flex flex-col items-center justify-center h-full gap-4">
        <Spinner color="secondary" size="lg" />
        <p className={subtitle()}>Redirecting to single molecule uploads…</p>
      </div>
    );
  }

  return (
    <>
      <SingleMoleculeControls />
      <SingleMoleculeLegends />
      <SingleMoleculeThreeScene key={dataset.id} />
    </>
  );
}

export default function SingleMoleculeViewerPage() {
  return (
    <Suspense
      fallback={
        <div className="flex items-center justify-center h-full">
          <Spinner size="lg" />
        </div>
      }
    >
      <SingleMoleculeViewerContent />
    </Suspense>
  );
}
