"use client";

import { Button } from "@heroui/button";
import { Slider } from "@heroui/slider";
import { useState } from "react";
import { useVisualizationStore } from "@/lib/stores/visualizationStore";
import { VisualizationPanel } from "./visualization-panel";
import type { VisualizationMode } from "@/lib/stores/visualizationStore";

export function VisualizationControls() {
  const { mode, setMode, sizeScale, setSizeScale } = useVisualizationStore();
  const [isPanelOpen, setIsPanelOpen] = useState(false);

  const handleModeChange = (newMode: VisualizationMode) => {
    if (mode === newMode) {
      // Toggle panel if clicking the same mode
      setIsPanelOpen(!isPanelOpen);
    } else {
      // Switch mode and open panel
      setMode(newMode);
      setIsPanelOpen(true);
    }
  };

  const buttonBaseClass = "w-14 h-14 min-w-0 rounded-full font-medium text-xs";

  return (
    <div className="fixed top-20 left-4 z-50 flex flex-col gap-2">
      {/* Celltype Button */}
      <Button
        className={buttonBaseClass}
        color={isPanelOpen && mode === "celltype" ? "primary" : "default"}
        variant={isPanelOpen && mode === "celltype" ? "shadow" : "flat"}
        onPress={() => handleModeChange("celltype")}
      >
        Celltype
      </Button>

      {/* Gene Button */}
      <Button
        className={buttonBaseClass}
        color={isPanelOpen && mode === "gene" ? "primary" : "default"}
        variant={isPanelOpen && mode === "gene" ? "shadow" : "flat"}
        onPress={() => handleModeChange("gene")}
      >
        Gene
      </Button>

      {/* Dot Size Slider */}
      <div className="w-14 h-32 bg-content1 rounded-full border-2 border-default-200 shadow-medium p-2 flex flex-col items-center justify-center">
        <Slider
          orientation="vertical"
          size="sm"
          step={0.1}
          minValue={0.1}
          maxValue={3}
          value={sizeScale}
          onChange={(value) => setSizeScale(value as number)}
          className="h-full"
          aria-label="Dot size"
        />
      </div>

      {/* Visualization Panel */}
      {isPanelOpen && (
        <VisualizationPanel mode={mode} onClose={() => setIsPanelOpen(false)} />
      )}
    </div>
  );
}
