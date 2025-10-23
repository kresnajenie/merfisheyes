"use client";

import { Button } from "@heroui/button";
import { Input } from "@heroui/input";
import { Checkbox } from "@heroui/checkbox";
import { Chip } from "@heroui/chip";
import { useState, useMemo } from "react";

import { useSingleMoleculeStore } from "@/lib/stores/singleMoleculeStore";
import { useSingleMoleculeVisualizationStore } from "@/lib/stores/singleMoleculeVisualizationStore";
import { glassButton } from "@/components/primitives";

export function SingleMoleculeControls() {
  const [isPanelOpen, setIsPanelOpen] = useState(false);
  const [searchTerm, setSearchTerm] = useState("");

  // Get dataset
  const dataset = useSingleMoleculeStore((state) => {
    const id = state.currentDatasetId;

    return id ? state.datasets.get(id) : null;
  });

  // Get visualization state
  const { selectedGenes, addGene, removeGene, clearGenes } =
    useSingleMoleculeVisualizationStore();

  // Filter genes based on search
  const filteredGenes = useMemo(() => {
    if (!dataset) return [];

    return dataset.uniqueGenes.filter((gene) =>
      gene.toLowerCase().includes(searchTerm.toLowerCase()),
    );
  }, [dataset, searchTerm]);

  // Convert Map to array for easier display
  const selectedGenesArray = Array.from(selectedGenes.keys());

  return (
    <div className="fixed top-20 left-4 z-50 flex flex-col gap-2">
      {/* Gene Selection Button */}
      <Button
        className={`w-14 h-14 min-w-0 rounded-full font-medium text-xs ${
          isPanelOpen ? "" : glassButton()
        }`}
        color={isPanelOpen ? "primary" : "default"}
        variant={isPanelOpen ? "shadow" : "light"}
        onPress={() => setIsPanelOpen(!isPanelOpen)}
      >
        Genes
      </Button>

      {/* Gene Selection Panel */}
      {isPanelOpen && (
        <div className="fixed top-20 left-20 z-50 w-[320px] bg-content1 border-2 border-default-200 rounded-lg shadow-lg">
          <div className="p-4 space-y-3">
            {/* Title */}
            <div className="flex items-center justify-between">
              <h3 className="text-lg font-semibold">Select Genes</h3>
              <span className="text-sm text-default-500">
                {selectedGenesArray.length} selected
              </span>
            </div>

            {/* Selected Genes Display */}
            {selectedGenesArray.length > 0 && (
              <div className="flex flex-wrap gap-1">
                {selectedGenesArray.map((gene) => {
                  const geneViz = selectedGenes.get(gene);

                  return (
                    <Chip
                      key={gene}
                      size="sm"
                      style={{
                        backgroundColor: geneViz?.color || "#808080",
                        color: "#000",
                      }}
                      variant="flat"
                      onClose={() => removeGene(gene)}
                    >
                      {gene}
                    </Chip>
                  );
                })}
              </div>
            )}

            {/* Search Input */}
            <Input
              classNames={{
                input: "text-sm",
              }}
              placeholder="Search genes..."
              size="sm"
              value={searchTerm}
              onValueChange={setSearchTerm}
            />

            {/* Clear Button */}
            <Button
              className="w-full"
              color="danger"
              size="sm"
              variant="ghost"
              onPress={() => {
                clearGenes();
                setSearchTerm("");
              }}
            >
              Clear All
            </Button>

            {/* Gene List */}
            <div className="max-h-[400px] overflow-y-auto flex flex-col gap-0">
              {filteredGenes.length > 0 ? (
                filteredGenes.map((gene) => (
                  <Checkbox
                    key={gene}
                    className="w-full"
                    isSelected={selectedGenes.has(gene)}
                    size="sm"
                    onValueChange={() => {
                      if (selectedGenes.has(gene)) {
                        removeGene(gene);
                      } else {
                        addGene(gene);
                      }
                    }}
                  >
                    <span className="text-sm">{gene}</span>
                  </Checkbox>
                ))
              ) : (
                <p className="text-sm text-default-400 text-center py-4">
                  {searchTerm ? "No genes found" : "No genes available"}
                </p>
              )}
            </div>

            {/* Gene count info */}
            <div className="text-xs text-default-400 text-center pt-2 border-t border-default-200">
              Showing {filteredGenes.length} of{" "}
              {dataset?.uniqueGenes.length || 0} genes
            </div>
          </div>
        </div>
      )}
    </div>
  );
}
