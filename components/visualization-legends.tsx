"use client";

import React, { useCallback, useEffect, useMemo, useState } from "react";
import { X } from "lucide-react";
import { Popover, PopoverTrigger, PopoverContent } from "@heroui/popover";
import { useVisualizationStore } from "@/lib/stores/visualizationStore";
import { useDatasetStore } from "@/lib/stores/datasetStore";
import type { StandardizedDataset } from "@/lib/StandardizedDataset";
import { GeneScalebar } from "@/components/gene-scalebar";
import {
  ColorPicker,
  ColorPickerEyeDropper,
  ColorPickerHue,
  ColorPickerOutput,
  ColorPickerSelection,
} from "@/components/ui/shadcn-io/color-picker";

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
    colorPalette: storeColorPalette,
    setColorPalette,
  } = useVisualizationStore();

  const { getCurrentDataset } = useDatasetStore();
  const dataset = getCurrentDataset();
  const [openPopoverCelltype, setOpenPopoverCelltype] = useState<string | null>(
    null
  );

  // Sync palette with selected cluster + local overrides
  useEffect(() => {
    if (!dataset || !selectedColumn) {
      setColorPalette({});
      return;
    }

    if (!("clusters" in dataset)) {
      setColorPalette({});
      return;
    }

    const standardizedDataset = dataset as StandardizedDataset;
    const selectedCluster = standardizedDataset.clusters?.find(
      (c) => c.column === selectedColumn
    );

    if (!selectedCluster) {
      setColorPalette({});
      return;
    }

    let palette = { ...(selectedCluster.palette || {}) };
    const storageKey = `cluster_palette_${standardizedDataset.id}_${selectedColumn}`;

    try {
      const stored = localStorage.getItem(storageKey);
      if (stored) {
        const overrides = JSON.parse(stored);
        palette = { ...palette, ...overrides };
      }
    } catch (error) {
      console.warn(
        "[VisualizationLegends] Failed to load palette overrides:",
        error
      );
    }

    selectedCluster.palette = palette;
    setColorPalette(palette);
  }, [dataset, selectedColumn, setColorPalette]);

  const handleClusterColorChange = useCallback(
    (celltype: string, color: string) => {
      if (!dataset || !selectedColumn) return;
      if (!("clusters" in dataset)) return;

      const standardizedDataset = dataset as StandardizedDataset;
      const selectedCluster = standardizedDataset.clusters?.find(
        (c) => c.column === selectedColumn
      );

      if (!selectedCluster) return;

      const updatedPalette = {
        ...(selectedCluster.palette || {}),
        [celltype]: color,
      };

      selectedCluster.palette = updatedPalette;
      setColorPalette(updatedPalette);

      const storageKey = `cluster_palette_${standardizedDataset.id}_${selectedColumn}`;
      try {
        localStorage.setItem(storageKey, JSON.stringify(updatedPalette));
      } catch (error) {
        console.warn(
          "[VisualizationLegends] Failed to persist palette overrides:",
          error
        );
      }
    },
    [dataset, selectedColumn, setColorPalette]
  );

  // Check if selected column is numerical
  const isNumericalColumn = useMemo(() => {
    if (!dataset || !selectedColumn) return false;
    if (!("clusters" in dataset)) return false;

    const standardizedDataset = dataset as StandardizedDataset;
    const selectedCluster = standardizedDataset.clusters?.find(
      (c) => c.column === selectedColumn
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
      dataset ? "clusters" in dataset : false
    );
    console.log(
      "dataset.clusters:",
      dataset && "clusters" in dataset
        ? (dataset as StandardizedDataset).clusters
        : "N/A"
    );
    console.log("selectedColumn:", selectedColumn);

    if (dataset && "clusters" in dataset) {
      const standardizedDataset = dataset as StandardizedDataset;
      const selectedCluster = standardizedDataset.clusters?.find(
        (c) => c.column === selectedColumn
      );
      console.log("selectedCluster:", selectedCluster);
      console.log("selectedCluster.palette:", selectedCluster?.palette);
    }

    console.log("colorPalette:", storeColorPalette);
    console.log("selectedCelltypes:", Array.from(selectedCelltypes));

    // Log colors for each selected celltype
    Array.from(selectedCelltypes).forEach((celltype) => {
      const color = storeColorPalette[celltype];
      console.log(`Color for "${celltype}":`, color || "NOT FOUND");
    });
  }, [dataset, storeColorPalette, selectedCelltypes, selectedColumn]);

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
              minValue={geneScaleMin}
              maxValue={geneScaleMax}
              onMinChange={setGeneScaleMin}
              onMaxChange={setGeneScaleMax}
            />
          ) : (
            <GeneScalebar
              minValue={numericalScaleMin}
              maxValue={numericalScaleMax}
              onMinChange={setNumericalScaleMin}
              onMaxChange={setNumericalScaleMax}
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
              const color = storeColorPalette[celltype] || "#888888";

              return (
                <Popover
                  key={celltype}
                  isOpen={openPopoverCelltype === celltype}
                  onOpenChange={(open) => {
                    if (!open) setOpenPopoverCelltype(null);
                  }}
                  placement="left"
                >
                  <PopoverTrigger>
                    <div
                      className="group flex items-center gap-2 px-4 py-2 rounded-full transition-all cursor-pointer"
                      style={{
                        backgroundColor: color,
                      }}
                      onClick={() => setOpenPopoverCelltype(celltype)}
                    >
                      <span className="text-xs font-medium text-black">
                        {celltype}
                      </span>
                      <X
                        className="w-2 h-2 text-black/70 group-hover:text-black"
                        onClick={(e) => {
                          e.stopPropagation();
                          toggleCelltype(celltype);
                        }}
                      />
                    </div>
                  </PopoverTrigger>
                  <PopoverContent className="!bg-[rgba(0,0,0,0.4)] backdrop-blur-[50px] border-white/20 p-3 w-64">
                    <div className="space-y-3">
                      <div className="flex items-center justify-between">
                        <h4 className="text-sm font-semibold">{celltype}</h4>
                        <button
                          type="button"
                          className="text-xs text-white/70 hover:text-white transition-colors"
                          onClick={() => setOpenPopoverCelltype(null)}
                        >
                          Close
                        </button>
                      </div>
                      <ColorPicker
                        key={celltype}
                        defaultValue={color}
                        onChange={(value) =>
                          handleClusterColorChange(celltype, value)
                        }
                        className="rounded-md p-2"
                      >
                        <ColorPickerSelection />
                        <div className="flex items-center gap-2 mt-2">
                          <ColorPickerEyeDropper />
                          <div className="grid w-full gap-1">
                            <ColorPickerHue />
                          </div>
                        </div>
                        <div className="mt-2 text-center">
                          <ColorPickerOutput />
                        </div>
                      </ColorPicker>
                    </div>
                  </PopoverContent>
                </Popover>
              );
            })}
          </div>
        </div>
      )}
    </div>
  );
};
