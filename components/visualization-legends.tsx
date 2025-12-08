"use client";

import type { StandardizedDataset } from "@/lib/StandardizedDataset";

import React, { useEffect, useMemo } from "react";
import { X } from "lucide-react";

import { useVisualizationStore } from "@/lib/stores/visualizationStore";
import { useDatasetStore } from "@/lib/stores/datasetStore";
import { GeneScalebar } from "@/components/gene-scalebar";

export const VisualizationLegends: React.FC = () => {
  const {
    selectedGene,
    setSelectedGene,
    selectedColumn,
    selectedCelltypes,
    toggleCelltype,
    mode,
    geneScaleMin,
    geneScaleMax,
    setGeneScaleMin,
    setGeneScaleMax,
    numericalScaleMin,
    numericalScaleMax,
    setNumericalScaleMin,
    setNumericalScaleMax,
  } = useVisualizationStore();

  const { getCurrentDataset } = useDatasetStore();
  const dataset = getCurrentDataset();

  // Get the palette from the dataset's cluster data
  const colorPalette = useMemo(() => {
    if (!dataset || !selectedColumn) return {};

    // Check if dataset has clusters property (StandardizedDataset)
    if (!("clusters" in dataset)) return {};

    // Type-narrow to StandardizedDataset
    const standardizedDataset = dataset as StandardizedDataset;
    const selectedCluster = standardizedDataset.clusters?.find(
      (c) => c.column === selectedColumn,
    );

    return selectedCluster?.palette || {};
  }, [dataset, selectedColumn]);

  // Check if selected column is numerical
  const isNumericalColumn = useMemo(() => {
    if (!dataset || !selectedColumn) return false;
    if (!("clusters" in dataset)) return false;

    const standardizedDataset = dataset as StandardizedDataset;
    const selectedCluster = standardizedDataset.clusters?.find(
      (c) => c.column === selectedColumn,
    );

    return selectedCluster?.type === "numerical";
  }, [dataset, selectedColumn]);

  // Debug: Log color palette whenever it changes
  useEffect(() => {
    console.log("=== VisualizationLegends Debug ===");
    console.log("dataset:", dataset);
    console.log("dataset.type:", dataset?.type);
    console.log(
      "has clusters property:",
      dataset ? "clusters" in dataset : false,
    );
    console.log(
      "dataset.clusters:",
      dataset && "clusters" in dataset
        ? (dataset as StandardizedDataset).clusters
        : "N/A",
    );
    console.log("selectedColumn:", selectedColumn);

    if (dataset && "clusters" in dataset) {
      const standardizedDataset = dataset as StandardizedDataset;
      const selectedCluster = standardizedDataset.clusters?.find(
        (c) => c.column === selectedColumn,
      );

      console.log("selectedCluster:", selectedCluster);
      console.log("selectedCluster.palette:", selectedCluster?.palette);
    }

    console.log("colorPalette:", colorPalette);
    console.log("selectedCelltypes:", Array.from(selectedCelltypes));

    // Log colors for each selected celltype
    Array.from(selectedCelltypes).forEach((celltype) => {
      const color = colorPalette[celltype];

      console.log(`Color for "${celltype}":`, color || "NOT FOUND");
    });
  }, [dataset, colorPalette, selectedCelltypes, selectedColumn]);

  // Don't render if nothing is selected
  const hasGene = !!selectedGene;
  const hasCelltypes = selectedCelltypes.size > 0;
  const showScalebar =
    (mode.includes("gene") && selectedGene) ||
    (mode.includes("celltype") &&
      !selectedGene &&
      selectedColumn &&
      isNumericalColumn);

  if (!hasGene && !hasCelltypes && !showScalebar) {
    return null;
  }

  return (
    <div className="fixed right-6 top-24 z-10 flex flex-col items-end gap-4 max-w-xs">
      {/* Selected Gene Badge */}
      {hasGene && (
        <div className="flex flex-col items-end gap-2">
          <div className="text-xs text-white/70 font-medium">Selected Gene</div>
          <div
            className="group flex items-center gap-2 px-4 py-2 rounded-full bg-blue-500/70 hover:bg-blue-500 transition-colors cursor-pointer"
            onClick={() => setSelectedGene(null)}
          >
            <span className="text-xs font-medium text-white">
              {selectedGene}
            </span>
            <X className="w-2 h-2 text-white/70 group-hover:text-white" />
          </div>
        </div>
      )}

      {/* Gene/Numerical Scale Bar */}
      {showScalebar && (
        <div className="flex flex-col items-end">
          {mode.includes("gene") && selectedGene ? (
            <GeneScalebar
              maxValue={geneScaleMax}
              minValue={geneScaleMin}
              onMaxChange={setGeneScaleMax}
              onMinChange={setGeneScaleMin}
            />
          ) : (
            <GeneScalebar
              maxValue={numericalScaleMax}
              minValue={numericalScaleMin}
              onMaxChange={setNumericalScaleMax}
              onMinChange={setNumericalScaleMin}
            />
          )}
        </div>
      )}

      {/* Selected Celltypes Badges */}
      {hasCelltypes && (
        <div className="flex flex-col items-end gap-2">
          <div
            className="text-xs text-white/70 font-medium cursor-pointer hover:text-white transition-colors"
            onClick={() => {
              // Clear all selected celltypes
              Array.from(selectedCelltypes).forEach((celltype) => {
                toggleCelltype(celltype);
              });
            }}
          >
            {selectedColumn || "Celltypes"} ({selectedCelltypes.size}) - Clear
            All
          </div>
          <div className="flex flex-col items-end gap-2 max-h-96 overflow-y-auto">
            {Array.from(selectedCelltypes).map((celltype) => {
              const color = colorPalette[celltype] || "#888888";

              return (
                <div
                  key={celltype}
                  className="group flex items-center gap-2 px-4 py-2 rounded-full transition-colors cursor-pointer"
                  style={{
                    backgroundColor: `${color}`, // 70% opacity (B3 in hex)
                  }}
                  onClick={() => toggleCelltype(celltype)}
                  onMouseEnter={(e) => {
                    e.currentTarget.style.backgroundColor = color;
                  }}
                  onMouseLeave={(e) => {
                    e.currentTarget.style.backgroundColor = `${color}B3`;
                  }}
                >
                  <span className="text-xs font-medium text-black">
                    {celltype}
                  </span>
                  <X className="w-2 h-2 text-black/70 group-hover:text-black" />
                </div>
              );
            })}
          </div>
        </div>
      )}
    </div>
  );
};
