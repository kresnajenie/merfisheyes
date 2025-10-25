"use client";

import type { VisualizationMode } from "@/lib/stores/visualizationStore";

import { Button } from "@heroui/button";
import { Slider } from "@heroui/slider";
import { Tooltip } from "@heroui/tooltip";
import { useState } from "react";

import { VisualizationPanel } from "./visualization-panel";

import { useVisualizationStore } from "@/lib/stores/visualizationStore";
import { glassButton } from "@/components/primitives";

export function VisualizationControls() {
  const { panelMode, setPanelMode, sizeScale, setSizeScale } =
    useVisualizationStore();
  const [isPanelOpen, setIsPanelOpen] = useState(false);

  const handleModeChange = (newMode: VisualizationMode) => {
    if (panelMode === newMode) {
      // Toggle panel if clicking the same mode
      setIsPanelOpen(!isPanelOpen);
    } else {
      // Switch panel mode and open panel
      setPanelMode(newMode);
      setIsPanelOpen(true);
    }
  };

  const buttonBaseClass = "w-14 h-14 min-w-0 rounded-full font-medium text-xs";

  return (
    <div className="fixed top-20 left-4 z-50 flex flex-col gap-2">
      {/* Celltype Button */}
      <Button
        className={`${buttonBaseClass} ${isPanelOpen && panelMode === "celltype" ? "" : glassButton()}`}
        color={isPanelOpen && panelMode === "celltype" ? "primary" : "default"}
        variant={isPanelOpen && panelMode === "celltype" ? "shadow" : "light"}
        onPress={() => handleModeChange("celltype")}
      >
        Celltype
      </Button>

      {/* Gene Button */}
      <Button
        className={`${buttonBaseClass} ${isPanelOpen && panelMode === "gene" ? "" : glassButton()}`}
        color={isPanelOpen && panelMode === "gene" ? "primary" : "default"}
        variant={isPanelOpen && panelMode === "gene" ? "shadow" : "light"}
        onPress={() => handleModeChange("gene")}
      >
        Gene
      </Button>

      {/* Dot Size Slider */}
      <Tooltip content="Change dotsize" placement="right">
        <div
          className={`w-14 h-32 rounded-full border-2 border-default-200 p-2 flex flex-col items-center justify-center ${glassButton()}`}
        >
          <Slider
            aria-label="Dot size"
            className="h-full"
            maxValue={3}
            minValue={0.1}
            orientation="vertical"
            size="sm"
            step={0.1}
            value={sizeScale}
            onChange={(value) => setSizeScale(value as number)}
          />
        </div>
      </Tooltip>

      {/* Visualization Panel */}
      {isPanelOpen && (
        <VisualizationPanel
          mode={panelMode}
          onClose={() => setIsPanelOpen(false)}
        />
      )}
    </div>
  );
}
