"use client";

import type { StandardizedDataset } from "@/lib/StandardizedDataset";

import React, { useCallback, useEffect, useMemo, useState } from "react";
import { X } from "lucide-react";
import { Popover, PopoverTrigger, PopoverContent } from "@heroui/popover";

import {
  usePanelVisualizationStore,
  usePanelDatasetStore,
} from "@/lib/hooks/usePanelStores";
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
    clusterVersion,
  } = usePanelVisualizationStore();

  const { getCurrentDataset } = usePanelDatasetStore();
  const dataset = getCurrentDataset();
  const [openPopoverCelltype, setOpenPopoverCelltype] = useState<string | null>(
    null,
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
      (c) => c.column === selectedColumn,
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
        error,
      );
    }

    selectedCluster.palette = palette;
    setColorPalette(palette);
  }, [dataset, selectedColumn, setColorPalette, clusterVersion]);

  const handleClusterColorChange = useCallback(
    (celltype: string, color: string) => {
      if (!dataset || !selectedColumn) return;
      if (!("clusters" in dataset)) return;

      const standardizedDataset = dataset as StandardizedDataset;
      const selectedCluster = standardizedDataset.clusters?.find(
        (c) => c.column === selectedColumn,
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
          error,
        );
      }
    },
    [dataset, selectedColumn, setColorPalette],
  );

  // Check if selected column is numerical
  const isNumericalColumn = useMemo(() => {
    if (!dataset || !selectedColumn) return false;
    if (!("clusters" in dataset)) return false;

    const standardizedDataset = dataset as StandardizedDataset;
    const selectedCluster = standardizedDataset.clusters?.find(
      (c) => c.column === selectedColumn,
    );

    return selectedCluster?.type === "numerical";
  }, [dataset, selectedColumn, clusterVersion]);

  // Debug logging removed for performance

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
    <div className="absolute right-6 top-24 z-10 flex flex-col items-end gap-4 max-w-xs">
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
          <div className="flex flex-col items-end gap-2 max-h-[calc(100vh-10rem)] overflow-y-auto">
            {Array.from(selectedCelltypes).map((celltype) => {
              const color = storeColorPalette[celltype] || "#888888";

              return (
                <Popover
                  key={celltype}
                  isOpen={openPopoverCelltype === celltype}
                  onOpenChange={(open) => {
                    if (!open) setOpenPopoverCelltype(null);
                  }}
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
                          className="text-xs text-white/70 hover:text-white transition-colors"
                          type="button"
                          onClick={() => setOpenPopoverCelltype(null)}
                        >
                          Close
                        </button>
                      </div>
                      <ColorPicker
                        key={celltype}
                        className="rounded-md p-2"
                        value={color}
                        onChange={(value) =>
                          handleClusterColorChange(celltype, value)
                        }
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
