"use client";

import type { StandardizedDataset } from "@/lib/StandardizedDataset";

import React, { useCallback, useEffect, useMemo, useState } from "react";
import { X, Globe, Eye, EyeOff } from "lucide-react";
import { Popover, PopoverTrigger, PopoverContent } from "@heroui/popover";
import { generateColorPalette } from "@/lib/utils/color-palette";

import {
  usePanelVisualizationStore,
  usePanelDatasetStore,
} from "@/lib/hooks/usePanelStores";
import { GeneScalebar } from "@/components/gene-scalebar";
import { ChannelScalebar } from "@/components/channel-scalebar";
import { getEffectiveColumnType } from "@/lib/utils/column-type-utils";
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
    geneEverywhere,
    setGeneEverywhere,
    hiddenCelltypes,
    toggleCelltypeVisibility,
    soloCelltype,
    toggleCelltype,
    setDegTarget,
    setDegTargetAuto,
    setDegPanelOpen,
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
    columnTypeOverrides,
    coexpressEnabled,
    selectedGene2,
    setSelectedGene2,
    coexpressSwapped,
    gene2ScaleMin,
    gene2ScaleMax,
    setGene2ScaleMin,
    setGene2ScaleMax,
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

    // If palette is missing (e.g. column was loaded as numerical), generate one
    let basePalette = selectedCluster.palette;

    if (
      !basePalette ||
      Object.keys(basePalette).length === 0
    ) {
      const uniqueValues = selectedCluster.uniqueValues
        ?? [...new Set(selectedCluster.values?.map(String) ?? [])];

      basePalette = generateColorPalette(uniqueValues);
    }

    let palette = { ...basePalette };
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

    return (
      getEffectiveColumnType(
        selectedColumn,
        dataset as StandardizedDataset,
        columnTypeOverrides,
      ) === "numerical"
    );
  }, [dataset, selectedColumn, clusterVersion, columnTypeOverrides]);

  // Whether this dataset supports DEG (drives the right-click-to-rank action)
  const hasDeStats = useMemo(() => {
    if (!dataset || !("clusters" in dataset)) return false;
    const ds = dataset as StandardizedDataset;

    return !!ds.deStats || (ds.availableDeStatsColumns?.length ?? 0) > 0;
  }, [dataset]);

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
      {hasGene && (() => {
        const coexpressActive = coexpressEnabled && selectedGene && selectedGene2;
        const gene1Color = coexpressActive
          ? (coexpressSwapped ? "#FF00FF" : "#00FF00")
          : null;
        const gene2Color = coexpressActive
          ? (coexpressSwapped ? "#00FF00" : "#FF00FF")
          : null;
        return (
        <div className="flex flex-col items-end gap-2">
          <div className="text-xs text-white/70 font-medium">
            {coexpressActive ? "Co-expression" : "Selected Gene"}
          </div>
          <div className="flex items-center gap-2">
            {/* Show-everywhere toggle — only meaningful while celltypes are
                selected (otherwise the gene already shows on every cell). */}
            {hasCelltypes && (
              <button
                className={`flex items-center justify-center w-9 h-9 rounded-full transition-colors ${
                  geneEverywhere
                    ? "bg-blue-500 text-white"
                    : "bg-white/10 text-white/60 hover:bg-white/20"
                }`}
                title={
                  geneEverywhere
                    ? "Showing gene on all cells — click to limit to selected celltypes"
                    : "Show gene expression on all cells (ignore celltype selection)"
                }
                type="button"
                onClick={() => setGeneEverywhere(!geneEverywhere)}
              >
                <Globe className="w-4 h-4" />
              </button>
            )}
            <div
              className="group flex items-center gap-2 px-4 py-2 rounded-full hover:opacity-100 transition-colors cursor-pointer"
              style={{
                backgroundColor: coexpressActive
                  ? `${gene1Color}b3` // 70% alpha
                  : "rgb(59 130 246 / 0.7)", // blue-500/70
              }}
              onClick={() => setSelectedGene(null)}
            >
              <span
                className="text-xs font-medium"
                style={{ color: coexpressActive ? "#000000" : "#ffffff" }}
              >
                {selectedGene}
              </span>
              <X className="w-2 h-2 text-black/60 group-hover:text-black" />
            </div>
          </div>
          {coexpressActive && (
            <div
              className="group flex items-center gap-2 px-4 py-2 rounded-full hover:opacity-100 transition-colors cursor-pointer"
              style={{ backgroundColor: `${gene2Color}b3` }}
              onClick={() => setSelectedGene2(null)}
            >
              <span className="text-xs font-medium text-black">
                {selectedGene2}
              </span>
              <X className="w-2 h-2 text-black/60 group-hover:text-black" />
            </div>
          )}
        </div>
        );
      })()}

      {/* Gene/Numerical Scale Bar */}
      {showScalebar && (
        <div className="flex flex-col items-end">
          {coexpressEnabled && selectedGene && selectedGene2 ? (
            <div className="flex items-start gap-3">
              <ChannelScalebar
                color={coexpressSwapped ? "#FF00FF" : "#00FF00"}
                label={selectedGene}
                maxValue={geneScaleMax}
                minValue={geneScaleMin}
                onMaxChange={setGeneScaleMax}
                onMinChange={setGeneScaleMin}
              />
              <ChannelScalebar
                color={coexpressSwapped ? "#00FF00" : "#FF00FF"}
                label={selectedGene2}
                maxValue={gene2ScaleMax}
                minValue={gene2ScaleMin}
                onMaxChange={setGene2ScaleMax}
                onMinChange={setGene2ScaleMin}
              />
            </div>
          ) : mode.includes("gene") && selectedGene ? (
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
              const hidden = hiddenCelltypes.has(celltype);

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
                        opacity: hidden ? 0.4 : 1,
                      }}
                      title={
                        hasDeStats
                          ? "Right-click to rank DEGs for this celltype"
                          : undefined
                      }
                      onClick={() => setOpenPopoverCelltype(celltype)}
                      onContextMenu={(e) => {
                        // Right-click → make this the DEG target and open
                        // the DEG panel. No-op (native menu) when the
                        // dataset has no DEG support.
                        if (!hasDeStats) return;
                        e.preventDefault();
                        setDegTarget(celltype);
                        // Manual mode so the right-clicked target sticks.
                        setDegTargetAuto(false);
                        setDegPanelOpen(true);
                      }}
                    >
                      <button
                        className="text-black/60 hover:text-black"
                        title={
                          hidden
                            ? "Hidden — click to show, ⌘-click to solo"
                            : "Visible — click to hide, ⌘-click to solo"
                        }
                        type="button"
                        onClick={(e) => {
                          e.stopPropagation();
                          if (e.metaKey || e.ctrlKey || e.altKey) {
                            soloCelltype(celltype);
                          } else {
                            toggleCelltypeVisibility(celltype);
                          }
                        }}
                      >
                        {hidden ? (
                          <EyeOff className="w-3 h-3" />
                        ) : (
                          <Eye className="w-3 h-3" />
                        )}
                      </button>
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
