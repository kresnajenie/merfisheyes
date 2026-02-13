"use client";

import type { VisualizationMode } from "@/lib/stores/createVisualizationStore";

import { Button } from "@heroui/button";
import { Slider } from "@heroui/react";
import { Tooltip } from "@heroui/tooltip";
import { useState, useRef } from "react";

import { VisualizationPanel } from "./visualization-panel";

import { usePanelVisualizationStore, usePanelId } from "@/lib/hooks/usePanelStores";
import { useSplitScreenStore } from "@/lib/stores/splitScreenStore";
import { glassButton } from "@/components/primitives";
import { VISUALIZATION_CONFIG } from "@/lib/config/visualization.config";

export function VisualizationControls() {
  const { panelMode, setPanelMode, sizeScale, setSizeScale } =
    usePanelVisualizationStore();
  const { isSplitMode, enableSplit } = useSplitScreenStore();
  const panelId = usePanelId();
  const [isPanelOpen, setIsPanelOpen] = useState(false);
  const controlsRef = useRef<HTMLDivElement>(null);

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
    <div
      ref={controlsRef}
      className="absolute top-28 left-4 z-50 flex flex-col gap-2"
    >
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
            maxValue={VISUALIZATION_CONFIG.SINGLE_CELL_SIZE_SCALE_MAX}
            minValue={VISUALIZATION_CONFIG.SINGLE_CELL_SIZE_SCALE_MIN}
            orientation="vertical"
            size="sm"
            step={VISUALIZATION_CONFIG.SINGLE_CELL_SIZE_SCALE_STEP}
            value={sizeScale}
            onChange={(value) => setSizeScale(value as number)}
          />
        </div>
      </Tooltip>

      {/* Split Screen Button — only show on left panel (no panelId) when not already in split mode */}
      {!isSplitMode && !panelId && (
        <Tooltip content="Split screen" placement="right">
          <Button
            className={`${buttonBaseClass} ${glassButton()}`}
            color="default"
            variant="light"
            onPress={enableSplit}
          >
            <svg
              className="w-5 h-5"
              fill="none"
              stroke="currentColor"
              strokeWidth={1.5}
              viewBox="0 0 24 24"
            >
              <path
                d="M9 4.5v15m6-15v15M4.5 19.5h15a1.5 1.5 0 001.5-1.5V6a1.5 1.5 0 00-1.5-1.5h-15A1.5 1.5 0 003 6v12a1.5 1.5 0 001.5 1.5z"
                strokeLinecap="round"
                strokeLinejoin="round"
              />
            </svg>
          </Button>
        </Tooltip>
      )}

      {/* Visualization Panel */}
      {isPanelOpen && (
        <VisualizationPanel
          controlsRef={controlsRef}
          mode={panelMode}
          onClose={() => setIsPanelOpen(false)}
        />
      )}
    </div>
  );
}
