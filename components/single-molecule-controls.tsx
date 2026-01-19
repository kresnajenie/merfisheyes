"use client";

import { Button } from "@heroui/button";
import { Input } from "@heroui/input";
import { Checkbox } from "@heroui/checkbox";
import { Slider } from "@heroui/react";
import { Tooltip } from "@heroui/tooltip";
import { useState, useMemo } from "react";

import { useSingleMoleculeStore } from "@/lib/stores/singleMoleculeStore";
import { useSingleMoleculeVisualizationStore } from "@/lib/stores/singleMoleculeVisualizationStore";
import { glassButton } from "@/components/primitives";
import { VISUALIZATION_CONFIG } from "@/lib/config/visualization.config";

export function SingleMoleculeControls() {
  const [isPanelOpen, setIsPanelOpen] = useState(false);
  const [searchTerm, setSearchTerm] = useState("");

  // Get dataset
  const dataset = useSingleMoleculeStore((state) => {
    const id = state.currentDatasetId;

    return id ? state.datasets.get(id) : null;
  });

  // Get visualization state
  const {
    selectedGenes,
    selectedGenesLegend,
    addGene,
    removeGene,
    clearGenes,
    globalScale,
    setGlobalScale,
    viewMode,
    setViewMode,
  } = useSingleMoleculeVisualizationStore();

  // Filter genes based on search
  const filteredGenes = useMemo(() => {
    if (!dataset) return [];

    return dataset.uniqueGenes.filter((gene) =>
      gene.toLowerCase().includes(searchTerm.toLowerCase()),
    );
  }, [dataset, searchTerm]);

  return (
    <div className="fixed top-28 left-4 z-50 flex flex-col gap-2">
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

      {/* Dot Size Slider */}
      <Tooltip content="Change dotsize" placement="right">
        <div
          className={`w-14 h-32 rounded-full border-2 border-default-200 p-2 flex flex-col items-center justify-center ${glassButton()}`}
        >
          <Slider
            aria-label="Dot size"
            className="h-full"
            maxValue={VISUALIZATION_CONFIG.SINGLE_MOLECULE_GLOBAL_SCALE_MAX}
            minValue={VISUALIZATION_CONFIG.SINGLE_MOLECULE_GLOBAL_SCALE_MIN}
            orientation="vertical"
            size="sm"
            step={VISUALIZATION_CONFIG.SINGLE_MOLECULE_GLOBAL_SCALE_STEP}
            value={globalScale}
            onChange={(value) => setGlobalScale(value as number)}
          />
        </div>
      </Tooltip>

      {/* 2D/3D View Toggle */}
      <Button
        className={`w-14 h-14 min-w-0 rounded-full font-medium text-xs ${glassButton()}`}
        color="default"
        variant="light"
        onPress={() => setViewMode(viewMode === "2D" ? "3D" : "2D")}
      >
        {viewMode}
      </Button>

      {/* Gene Selection Panel */}
      {isPanelOpen && (
        <div
          className={`fixed top-28 left-20 z-50 w-[320px] border-2 border-white/20 rounded-3xl shadow-lg ${glassButton()}`}
        >
          <div className="p-4 space-y-3">
            {/* Title */}
            <div className="flex items-center justify-between">
              <h3 className="text-lg font-semibold">Select Genes</h3>
              <span className="text-sm text-default-500">
                {selectedGenesLegend.size} selected
              </span>
            </div>

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
                        const geneProps = dataset?.geneColors[gene];

                        if (geneProps) {
                          addGene(gene, geneProps.color, geneProps.size);
                        } else {
                          addGene(gene);
                        }
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
