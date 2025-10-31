"use client";

import React, { useEffect, useState } from "react";
import { X } from "lucide-react";
import { Checkbox } from "@heroui/checkbox";
import { Tooltip } from "@heroui/tooltip";
import { Popover, PopoverTrigger, PopoverContent } from "@heroui/popover";
import { Button } from "@heroui/button";
import { Slider } from "@heroui/slider";
import { useSingleMoleculeVisualizationStore } from "@/lib/stores/singleMoleculeVisualizationStore";
import { useSingleMoleculeStore } from "@/lib/stores/singleMoleculeStore";
import { VISUALIZATION_CONFIG } from "@/lib/config/visualization.config";
import {
  ColorPicker,
  ColorPickerAlpha,
  ColorPickerEyeDropper,
  ColorPickerFormat,
  ColorPickerHue,
  ColorPickerOutput,
  ColorPickerSelection,
} from '@/components/ui/shadcn-io/color-picker';

export const SingleMoleculeLegends: React.FC = () => {
  const { selectedGenes, selectedGenesLegend, geneDataCache, removeGene, toggleGeneVisibility, setGeneColor, setGeneLocalScale } = useSingleMoleculeVisualizationStore();
  const dataset = useSingleMoleculeStore((state) => {
    const id = state.currentDatasetId;
    return id ? state.datasets.get(id) : null;
  });

  const [openPopoverGene, setOpenPopoverGene] = useState<string | null>(null);
  const [tempColors, setTempColors] = useState<Record<string, string>>({});
  const [tempScales, setTempScales] = useState<Record<string, number>>({});

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

  const handleRightClick = (e: React.MouseEvent, gene: string, geneViz: any) => {
    e.preventDefault();
    setOpenPopoverGene(gene);
    setTempColors(prev => ({ ...prev, [gene]: geneViz.color }));
    setTempScales(prev => ({ ...prev, [gene]: geneViz.localScale }));
  };

  const handleSaveChanges = (gene: string) => {
    if (dataset && tempColors[gene] && tempScales[gene] !== undefined) {
      // Update visualization store
      setGeneColor(gene, tempColors[gene]);
      setGeneLocalScale(gene, tempScales[gene]);

      // Update dataset localStorage
      dataset.geneColors[gene] = {
        color: tempColors[gene],
        size: tempScales[gene],
      };
      const storageKey = `sm_gene_colors_${dataset.id}`;
      try {
        localStorage.setItem(storageKey, JSON.stringify(dataset.geneColors));
        console.log(`[SingleMoleculeLegends] Updated gene colors in localStorage for ${gene}`);
      } catch (error) {
        console.warn(`[SingleMoleculeLegends] Failed to update localStorage:`, error);
      }

      setOpenPopoverGene(null);
    }
  };

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
              <Popover
                key={gene}
                isOpen={openPopoverGene === gene}
                onOpenChange={(open) => {
                  if (!open) setOpenPopoverGene(null);
                }}
                placement="left"
              >
                <Tooltip
                  content="Right click to change color and size"
                  placement="left"
                  delay={0}
                >
                  <PopoverTrigger>
                    <div
                      className="group flex items-center gap-2 px-4 py-2 rounded-full transition-all cursor-pointer"
                      style={{
                        backgroundColor: isVisible ? geneViz.color : `${geneViz.color}80`, // 50% opacity when hidden
                      }}
                      onContextMenu={(e) => handleRightClick(e, gene, geneViz)}
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
                        onClick={(e) => {
                          e.stopPropagation();
                          console.log("[SingleMoleculeLegends] Removing gene from legend:", gene);
                          removeGene(gene);
                        }}
                      />
                    </div>
                  </PopoverTrigger>
                </Tooltip>
                <PopoverContent className="p-4 w-80">
                  <div className="space-y-4">
                    <h4 className="font-semibold">Edit: {gene}</h4>

                    {/* Color Picker */}
                    <div>
                      <label className="text-sm font-medium mb-2 block">Color</label>
                      <ColorPicker
                        value={tempColors[gene] || geneViz.color}
                        onChange={(color) => setTempColors(prev => ({ ...prev, [gene]: color }))}
                        className="rounded-md border bg-background p-3 shadow-sm"
                      >
                        <ColorPickerSelection />
                        <div className="flex items-center gap-3 mt-3">
                          <ColorPickerEyeDropper />
                          <div className="grid w-full gap-1">
                            <ColorPickerHue />
                            <ColorPickerAlpha />
                          </div>
                        </div>
                        <div className="flex items-center gap-2 mt-3">
                          <ColorPickerOutput />
                          <ColorPickerFormat />
                        </div>
                      </ColorPicker>
                    </div>

                    {/* Local Dot Size Slider */}
                    <div>
                      <label className="text-sm font-medium mb-2 block">
                        Local Dot Size: {(tempScales[gene] || geneViz.localScale).toFixed(1)}x
                      </label>
                      <Slider
                        value={tempScales[gene] || geneViz.localScale}
                        onChange={(value) => setTempScales(prev => ({ ...prev, [gene]: value as number }))}
                        minValue={VISUALIZATION_CONFIG.SINGLE_MOLECULE_LOCAL_SCALE_MIN}
                        maxValue={VISUALIZATION_CONFIG.SINGLE_MOLECULE_LOCAL_SCALE_MAX}
                        step={VISUALIZATION_CONFIG.SINGLE_MOLECULE_LOCAL_SCALE_STEP}
                        aria-label="Local dot size"
                        className="w-full"
                        size="sm"
                      />
                      <div className="flex justify-between text-xs text-default-400 mt-1">
                        <span>{VISUALIZATION_CONFIG.SINGLE_MOLECULE_LOCAL_SCALE_MIN}x</span>
                        <span>{VISUALIZATION_CONFIG.SINGLE_MOLECULE_LOCAL_SCALE_MAX}x</span>
                      </div>
                    </div>

                    {/* Buttons */}
                    <div className="flex justify-end gap-2">
                      <Button
                        size="sm"
                        variant="light"
                        onPress={() => setOpenPopoverGene(null)}
                      >
                        Cancel
                      </Button>
                      <Button
                        size="sm"
                        color="primary"
                        onPress={() => handleSaveChanges(gene)}
                      >
                        Save
                      </Button>
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
