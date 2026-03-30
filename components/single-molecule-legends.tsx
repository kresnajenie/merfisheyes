"use client";

import React, { useEffect, useState } from "react";
import { X } from "lucide-react";
import { Checkbox } from "@heroui/checkbox";
import { Switch } from "@heroui/switch";
import { Tooltip } from "@heroui/tooltip";
import { Popover, PopoverTrigger, PopoverContent } from "@heroui/popover";
import { Button } from "@heroui/button";
import { Slider } from "@heroui/react";

import {
  usePanelSingleMoleculeVisualizationStore,
  usePanelSingleMoleculeStore,
} from "@/lib/hooks/usePanelStores";
import { VISUALIZATION_CONFIG } from "@/lib/config/visualization.config";
import type { MoleculeShape } from "@/lib/stores/createSingleMoleculeVisualizationStore";
import {
  ColorPicker,
  ColorPickerEyeDropper,
  ColorPickerHue,
  ColorPickerOutput,
  ColorPickerSelection,
} from "@/components/ui/shadcn-io/color-picker";

export const SingleMoleculeLegends: React.FC = () => {
  const {
    selectedGenes,
    selectedGenesLegend,
    geneDataCache,
    removeGene,
    toggleGeneVisibility,
    setGeneColor,
    setGeneLocalScale,
    setGeneShowAssigned,
    setGeneShowUnassigned,
    setGeneAssignedShape,
    setGeneUnassignedShape,
    setGeneUnassignedColor,
    setGeneUnassignedLocalScale,
  } = usePanelSingleMoleculeVisualizationStore();
  const dataset = usePanelSingleMoleculeStore((state) => {
    const id = state.currentDatasetId;

    return id ? state.datasets.get(id) : null;
  });

  const [openPopoverGene, setOpenPopoverGene] = useState<string | null>(null);
  const [popoverMode, setPopoverMode] = useState<"assigned" | "unassigned">(
    "assigned",
  );

  // Debug logging
  useEffect(() => {
    console.log("[SingleMoleculeLegends] Component mounted/updated");
    console.log(
      "[SingleMoleculeLegends] selectedGenesLegend size:",
      selectedGenesLegend.size,
    );
    console.log(
      "[SingleMoleculeLegends] selectedGenes size:",
      selectedGenes.size,
    );
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

  console.log(
    "[SingleMoleculeLegends] Rendering",
    legendGenesArray.length,
    "gene badges",
  );

  const handleRightClick = (e: React.MouseEvent, gene: string) => {
    e.preventDefault();
    setOpenPopoverGene(gene);
  };

  const handleColorChange = (gene: string, color: string) => {
    if (dataset) {
      // Update visualization store immediately
      setGeneColor(gene, color);

      // Update dataset localStorage
      dataset.geneColors[gene] = {
        ...dataset.geneColors[gene],
        color: color,
      };
      const storageKey = `sm_gene_colors_${dataset.id}`;

      try {
        localStorage.setItem(storageKey, JSON.stringify(dataset.geneColors));
      } catch (error) {
        console.warn(
          `[SingleMoleculeLegends] Failed to update localStorage:`,
          error,
        );
      }
    }
  };

  const handleScaleChange = (gene: string, scale: number) => {
    if (dataset) {
      // Update visualization store immediately
      setGeneLocalScale(gene, scale);

      // Update dataset localStorage
      dataset.geneColors[gene] = {
        ...dataset.geneColors[gene],
        size: scale,
      };
      const storageKey = `sm_gene_colors_${dataset.id}`;

      try {
        localStorage.setItem(storageKey, JSON.stringify(dataset.geneColors));
      } catch (error) {
        console.warn(
          `[SingleMoleculeLegends] Failed to update localStorage:`,
          error,
        );
      }
    }
  };

  const handleUnassignedColorChange = (gene: string, color: string) => {
    setGeneUnassignedColor(gene, color);
  };

  return (
    <div className="absolute right-6 top-24 z-10 flex flex-col items-end gap-4 max-w-xs">
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
        <div className="flex flex-col items-end gap-2 max-h-[calc(100vh-10rem)] overflow-y-auto">
          {legendGenesArray.map(({ gene, geneViz, isVisible }) => {
            if (!geneViz) return null;

            console.log(
              "[SingleMoleculeLegends] Rendering badge for gene:",
              gene,
              "color:",
              geneViz.color,
              "visible:",
              isVisible,
            );

            return (
              <Popover
                key={gene}
                isOpen={openPopoverGene === gene}
                placement="left"
                onOpenChange={(open) => {
                  if (!open) setOpenPopoverGene(null);
                }}
              >
                <Tooltip
                  content="Click to change color and size"
                  delay={0}
                  placement="left"
                >
                  <PopoverTrigger>
                    <div
                      className="group flex items-center gap-2 px-4 py-2 rounded-full transition-all cursor-pointer"
                      style={{
                        backgroundColor: isVisible
                          ? geneViz.color
                          : `${geneViz.color}80`, // 50% opacity when hidden
                      }}
                      onClick={() => setOpenPopoverGene(gene)}
                    >
                      {/* Checkbox for visibility toggle */}
                      <Checkbox
                        classNames={{
                          wrapper: "bg-white/20",
                        }}
                        isSelected={isVisible}
                        size="sm"
                        onValueChange={() => {
                          console.log(
                            "[SingleMoleculeLegends] Toggling visibility for gene:",
                            gene,
                          );
                          toggleGeneVisibility(gene);
                        }}
                      />

                      {/* Gene name with strikethrough when hidden */}
                      <span
                        className={`text-xs font-medium text-black ${!isVisible ? "line-through" : ""}`}
                      >
                        {gene}
                      </span>

                      {/* X button to remove from legend */}
                      <X
                        className="w-3 h-3 text-black/70 group-hover:text-black cursor-pointer transition-colors"
                        onClick={(e) => {
                          e.stopPropagation();
                          console.log(
                            "[SingleMoleculeLegends] Removing gene from legend:",
                            gene,
                          );
                          removeGene(gene);
                        }}
                      />
                    </div>
                  </PopoverTrigger>
                </Tooltip>
                <PopoverContent className="!bg-[rgba(0,0,0,0.4)] backdrop-blur-[50px] border-white/20 p-3 w-64">
                  <div className="space-y-3">
                    {/* Header */}
                    <div className="flex items-center justify-between">
                      <h4 className="text-sm font-semibold">{gene}</h4>
                      <Button
                        isIconOnly
                        className="h-6 w-6 min-w-0"
                        size="sm"
                        variant="light"
                        onPress={() => setOpenPopoverGene(null)}
                      >
                        <X className="w-3 h-3" />
                      </Button>
                    </div>

                    {/* Mode Toggle: Assigned / Unassigned */}
                    {dataset?.hasUnassigned && (
                      <div className="flex gap-1 items-center">
                        <div className="flex gap-1 flex-1">
                          <Button
                            className="flex-1 text-xs"
                            color={
                              popoverMode === "assigned"
                                ? "primary"
                                : "default"
                            }
                            size="sm"
                            variant={
                              popoverMode === "assigned" ? "solid" : "bordered"
                            }
                            onPress={() => setPopoverMode("assigned")}
                          >
                            Assigned
                          </Button>
                          <Button
                            className="flex-1 text-xs"
                            color={
                              popoverMode === "unassigned"
                                ? "primary"
                                : "default"
                            }
                            size="sm"
                            variant={
                              popoverMode === "unassigned"
                                ? "solid"
                                : "bordered"
                            }
                            onPress={() => setPopoverMode("unassigned")}
                          >
                            Unassigned
                          </Button>
                        </div>
                        {/* Per-gene visibility toggle for current mode */}
                        <Tooltip
                          content={
                            popoverMode === "assigned"
                              ? geneViz.showAssigned
                                ? "Hide assigned"
                                : "Show assigned"
                              : geneViz.showUnassigned
                                ? "Hide unassigned"
                                : "Show unassigned"
                          }
                          placement="top"
                        >
                          <div>
                            <Switch
                              isSelected={
                                popoverMode === "assigned"
                                  ? geneViz.showAssigned
                                  : geneViz.showUnassigned
                              }
                              size="sm"
                              onValueChange={(val) =>
                                popoverMode === "assigned"
                                  ? setGeneShowAssigned(gene, val)
                                  : setGeneShowUnassigned(gene, val)
                              }
                            />
                          </div>
                        </Tooltip>
                      </div>
                    )}

                    {/* Color Picker — switches between assigned/unassigned */}
                    <ColorPicker
                      key={`${gene}-${popoverMode}`}
                      className="rounded-md p-2"
                      defaultValue={
                        popoverMode === "unassigned"
                          ? geneViz.unassignedColor
                          : geneViz.color
                      }
                      onChange={(color) =>
                        popoverMode === "unassigned"
                          ? handleUnassignedColorChange(gene, color)
                          : handleColorChange(gene, color)
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

                    {/* Size Slider — switches between assigned/unassigned */}
                    <div>
                      <label className="text-xs font-medium mb-1 block">
                        Size:{" "}
                        {popoverMode === "unassigned"
                          ? geneViz.unassignedLocalScale.toFixed(1)
                          : geneViz.localScale.toFixed(1)}
                        x
                      </label>
                      <Slider
                        aria-label="Local dot size"
                        className="w-full"
                        maxValue={
                          VISUALIZATION_CONFIG.SINGLE_MOLECULE_LOCAL_SCALE_MAX
                        }
                        minValue={
                          VISUALIZATION_CONFIG.SINGLE_MOLECULE_LOCAL_SCALE_MIN
                        }
                        size="sm"
                        step={
                          VISUALIZATION_CONFIG.SINGLE_MOLECULE_LOCAL_SCALE_STEP
                        }
                        value={
                          popoverMode === "unassigned"
                            ? geneViz.unassignedLocalScale
                            : geneViz.localScale
                        }
                        onChange={(value) =>
                          popoverMode === "unassigned"
                            ? setGeneUnassignedLocalScale(
                                gene,
                                value as number,
                              )
                            : handleScaleChange(gene, value as number)
                        }
                      />
                      <div className="flex justify-between text-[10px] text-default-400 mt-1">
                        <span>
                          {VISUALIZATION_CONFIG.SINGLE_MOLECULE_LOCAL_SCALE_MIN}
                          x
                        </span>
                        <span>
                          {VISUALIZATION_CONFIG.SINGLE_MOLECULE_LOCAL_SCALE_MAX}
                          x
                        </span>
                      </div>
                    </div>

                    {/* Shape Selector — switches between assigned/unassigned */}
                    <div>
                      <label className="text-xs font-medium mb-1 block">
                        Shape
                      </label>
                      <div className="flex gap-1">
                        {(["circle", "square"] as MoleculeShape[]).map(
                          (shape) => {
                            const currentShape =
                              popoverMode === "unassigned"
                                ? geneViz.unassignedShape
                                : geneViz.assignedShape;

                            return (
                              <Button
                                key={shape}
                                className="flex-1 capitalize text-xs"
                                color={
                                  currentShape === shape
                                    ? "primary"
                                    : "default"
                                }
                                size="sm"
                                variant={
                                  currentShape === shape
                                    ? "solid"
                                    : "bordered"
                                }
                                onPress={() =>
                                  popoverMode === "unassigned"
                                    ? setGeneUnassignedShape(gene, shape)
                                    : setGeneAssignedShape(gene, shape)
                                }
                              >
                                {shape}
                              </Button>
                            );
                          },
                        )}
                      </div>
                    </div>
                  </div>
                </PopoverContent>
              </Popover>
            );
          })}
        </div>
      </div>
    </div>
  );
};
