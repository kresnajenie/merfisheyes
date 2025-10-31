"use client";

import React, { useEffect } from "react";
import { X } from "lucide-react";
import { Checkbox } from "@heroui/checkbox";
import { useSingleMoleculeVisualizationStore } from "@/lib/stores/singleMoleculeVisualizationStore";

export const SingleMoleculeLegends: React.FC = () => {
  const { selectedGenes, selectedGenesLegend, geneDataCache, removeGene, toggleGeneVisibility } = useSingleMoleculeVisualizationStore();

  // Debug logging
  useEffect(() => {
    console.log("[SingleMoleculeLegends] Component mounted/updated");
    console.log("[SingleMoleculeLegends] selectedGenesLegend size:", selectedGenesLegend.size);
    console.log("[SingleMoleculeLegends] selectedGenes size:", selectedGenes.size);
  }, [selectedGenesLegend, selectedGenes]);

  // Don't render if no genes are in legend
  if (selectedGenesLegend.size === 0) {
    console.log("[SingleMoleculeLegends] No genes in legend, returning null");
    return null;
  }

  // Create array of legend genes with their data and visibility state
  const legendGenesArray = Array.from(selectedGenesLegend).map((gene) => {
    const geneViz = geneDataCache.get(gene);
    const isVisible = selectedGenes.has(gene);
    return { gene, geneViz, isVisible };
  });

  console.log("[SingleMoleculeLegends] Rendering", legendGenesArray.length, "gene badges");

  return (
    <div className="fixed right-6 top-24 z-10 flex flex-col items-end gap-4 max-w-xs">
      {/* Selected Genes Badges */}
      <div className="flex flex-col items-end gap-2">
        <div
          className="text-xs text-white/70 font-medium cursor-pointer hover:text-white transition-colors"
          onClick={() => {
            console.log("[SingleMoleculeLegends] Clear All clicked");
            // Clear all legend genes
            Array.from(selectedGenesLegend).forEach((gene) => {
              removeGene(gene);
            });
          }}
        >
          Selected Genes ({selectedGenesLegend.size}) - Clear All
        </div>
        <div className="flex flex-col items-end gap-2 max-h-96 overflow-y-auto">
          {legendGenesArray.map(({ gene, geneViz, isVisible }) => {
            if (!geneViz) return null;

            console.log("[SingleMoleculeLegends] Rendering badge for gene:", gene, "color:", geneViz.color, "visible:", isVisible);

            return (
              <div
                key={gene}
                className="group flex items-center gap-2 px-4 py-2 rounded-full transition-all"
                style={{
                  backgroundColor: isVisible ? geneViz.color : `${geneViz.color}80`, // 50% opacity when hidden
                }}
              >
                {/* Checkbox for visibility toggle */}
                <Checkbox
                  size="sm"
                  isSelected={isVisible}
                  onValueChange={() => {
                    console.log("[SingleMoleculeLegends] Toggling visibility for gene:", gene);
                    toggleGeneVisibility(gene);
                  }}
                  classNames={{
                    wrapper: "bg-white/20",
                  }}
                />

                {/* Gene name with strikethrough when hidden */}
                <span
                  className={`text-xs font-medium text-black ${!isVisible ? 'line-through' : ''}`}
                >
                  {gene}
                </span>

                {/* X button to remove from legend */}
                <X
                  className="w-3 h-3 text-black/70 group-hover:text-black cursor-pointer transition-colors"
                  onClick={() => {
                    console.log("[SingleMoleculeLegends] Removing gene from legend:", gene);
                    removeGene(gene);
                  }}
                />
              </div>
            );
          })}
        </div>
      </div>
    </div>
  );
};
