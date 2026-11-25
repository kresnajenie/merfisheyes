"use client";

import { Slider } from "@heroui/react";
import { Tooltip } from "@heroui/tooltip";

import { glassButton } from "@/components/primitives";
import { usePanelVisualizationStore } from "@/lib/hooks/usePanelStores";

/**
 * The cells' opacity, for this panel: every cell, spot or bin, over its own
 * colouring. Below 100% what lies under the cells (an H&E under Visium HD
 * bins, which tile it) shows through. Sits beside the dot-size slider and
 * looks like the molecules' one.
 */
export function CellOpacitySlider() {
  const opacity = usePanelVisualizationStore((s) => s.alphaScale);
  const setOpacity = usePanelVisualizationStore((s) => s.setAlphaScale);

  return (
    <Tooltip
      content={`Cell opacity: ${Math.round(opacity * 100)}%`}
      placement="right"
    >
      <div
        className={`w-14 h-32 rounded-full border-2 border-default-200 p-2 flex flex-col items-center justify-center ${glassButton()}`}
        data-testid="cell-opacity"
      >
        <Slider
          aria-label="Cell opacity"
          className="h-full"
          maxValue={1}
          minValue={0.05}
          orientation="vertical"
          size="sm"
          step={0.05}
          value={opacity}
          onChange={(value) => setOpacity(value as number)}
        />
      </div>
    </Tooltip>
  );
}
