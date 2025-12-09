"use client";

import { Suspense, useEffect } from "react";
import { useRouter } from "next/navigation";

import { SingleMoleculeThreeScene } from "@/components/single-molecule-three-scene";
import { SingleMoleculeControls } from "@/components/single-molecule-controls";
import { SingleMoleculeLegends } from "@/components/single-molecule-legends";
import { FileUpload } from "@/components/file-upload";
import { useSingleMoleculeStore } from "@/lib/stores/singleMoleculeStore";
import { useSingleMoleculeVisualizationStore } from "@/lib/stores/singleMoleculeVisualizationStore";
import LightRays from "@/components/react-bits/LightRays";
import { subtitle, title } from "@/components/primitives";

function ViewerContent() {
  const router = useRouter();
  // Get dataset from single molecule store
  const dataset = useSingleMoleculeStore((state) => {
    const id = state.currentDatasetId;

    return id ? state.datasets.get(id) : null;
  });

  const { addGene, clearGenes } = useSingleMoleculeVisualizationStore();

  // Auto-select LEF1 + 4 random genes when dataset loads
  useEffect(() => {
    if (!dataset) return;

    console.log("[sm-viewer] Dataset loaded, selecting LEF1 + 4 random genes");
    console.log("Available genes:", dataset.uniqueGenes.length);

    // Clear previous selections
    clearGenes();

    const selectedGenes: string[] = [];

    // Always add LEF1 if it exists
    if (dataset.uniqueGenes.includes("LEF1")) {
      selectedGenes.push("LEF1");
      console.log("Added LEF1 (always included)");
    } else {
      console.warn("LEF1 not found in dataset");
    }

    // Pick 4 additional random genes
    const availableGenes = dataset.uniqueGenes.filter(
      (gene) => gene !== "LEF1",
    );
    const numRandomGenes = Math.min(4, availableGenes.length);

    for (let i = 0; i < numRandomGenes; i++) {
      const randomIndex = Math.floor(Math.random() * availableGenes.length);

      selectedGenes.push(availableGenes[randomIndex]);
      availableGenes.splice(randomIndex, 1); // Remove to avoid duplicates
    }

    console.log("Selected genes:", selectedGenes);

    // Add genes to visualization store with their persistent colors
    selectedGenes.forEach((gene) => {
      const geneProps = dataset.geneColors[gene];
      addGene(gene, geneProps.color, geneProps.size);
    });
  }, [dataset, addGene, clearGenes]);

  useEffect(() => {
    if (!dataset) {
      router.replace("/?mode=sm");
    }
  }, [dataset, router]);

  if (!dataset) {
    return null;
  }

  return (
    <>
      <SingleMoleculeControls />
      <SingleMoleculeLegends />
      <SingleMoleculeThreeScene />
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
