"use client";

import React, { useEffect } from "react";
import { X } from "lucide-react";
import { useSingleMoleculeVisualizationStore } from "@/lib/stores/singleMoleculeVisualizationStore";

export const SingleMoleculeLegends: React.FC = () => {
  const { selectedGenes, removeGene } = useSingleMoleculeVisualizationStore();

  // Debug logging
  useEffect(() => {
    console.log("[SingleMoleculeLegends] Component mounted/updated");
    console.log("[SingleMoleculeLegends] selectedGenes size:", selectedGenes.size);
    console.log("[SingleMoleculeLegends] selectedGenes entries:", Array.from(selectedGenes.entries()));
  }, [selectedGenes]);

  // Don't render if no genes are selected
  if (selectedGenes.size === 0) {
    console.log("[SingleMoleculeLegends] No genes selected, returning null");
    return null;
  }

  const genesArray = Array.from(selectedGenes.entries());
  console.log("[SingleMoleculeLegends] Rendering", genesArray.length, "gene badges");

  return (
    <div className="fixed right-6 top-24 z-10 flex flex-col items-end gap-4 max-w-xs">
      {/* Selected Genes Badges */}
      <div className="flex flex-col items-end gap-2">
        <div
          className="text-xs text-white/70 font-medium cursor-pointer hover:text-white transition-colors"
          onClick={() => {
            console.log("[SingleMoleculeLegends] Clear All clicked");
            // Clear all selected genes
            Array.from(selectedGenes.keys()).forEach((gene) => {
              removeGene(gene);
            });
          }}
        >
          Selected Genes ({selectedGenes.size}) - Clear All
        </div>
        <div className="flex flex-col items-end gap-2 max-h-96 overflow-y-auto">
          {genesArray.map(([gene, geneViz]) => {
            console.log("[SingleMoleculeLegends] Rendering badge for gene:", gene, "color:", geneViz.color);
            return (
              <div
                key={gene}
                className="group flex items-center gap-2 px-4 py-2 rounded-full transition-colors cursor-pointer"
                style={{
                  backgroundColor: geneViz.color, // 100% opacity
                }}
                onClick={() => {
                  console.log("[SingleMoleculeLegends] Removing gene:", gene);
                  removeGene(gene);
                }}
              >
                <span className="text-xs font-medium text-black">
                  {gene}
                </span>
                <X className="w-3 h-3 text-black/70 group-hover:text-black" />
              </div>
            );
          })}
        </div>
      </div>
    </div>
  );
};
