"use client";

import React, { useState } from "react";

import { NumberScrubber } from "@/components/ui/number-scrubber";
import { VISUALIZATION_CONFIG } from "@/lib/config/visualization.config";
import { colormapCss } from "@/lib/utils/colormaps";
import { usePanelVisualizationStore } from "@/lib/hooks/usePanelStores";
import { ColormapPicker } from "./colormap-picker";

interface GeneScalebarProps {
  minValue?: number;
  maxValue?: number;
  onMinChange?: (value: number) => void;
  onMaxChange?: (value: number) => void;
}

export const GeneScalebar: React.FC<GeneScalebarProps> = ({
  minValue = VISUALIZATION_CONFIG.SCALE_BAR_DEFAULT_MIN,
  maxValue = VISUALIZATION_CONFIG.SCALE_BAR_DEFAULT_MAX,
  onMinChange,
  onMaxChange,
}) => {
  // Dynamic step size: configurable percentage of max value
  const dynamicStep = Math.max(
    VISUALIZATION_CONFIG.SCALE_BAR_MIN_STEP,
    maxValue * VISUALIZATION_CONFIG.SCALE_BAR_STEP_PERCENTAGE,
  );

  const colormap = usePanelVisualizationStore((s) => s.colormap);
  const setColormap = usePanelVisualizationStore((s) => s.setColormap);
  const [pickerPos, setPickerPos] = useState<{ x: number; y: number } | null>(
    null,
  );

  return (
    <div className="flex flex-col items-center gap-2">
      {/* Top number scrubber */}
      <NumberScrubber
        className="text-white drop-shadow-lg"
        decimals={VISUALIZATION_CONFIG.SCALE_BAR_DECIMALS}
        min={minValue}
        step={dynamicStep}
        value={maxValue}
        onChange={onMaxChange}
      />

      {/* Gradient bar — click to change colormap */}
      <button
        aria-label="Change colormap"
        className="w-8 h-32 rounded-lg shadow-lg cursor-pointer hover:ring-2 hover:ring-white/40 transition-shadow"
        style={{ background: colormapCss(colormap, "to bottom") }}
        title="Click to change colormap"
        type="button"
        onClick={(e) => {
          e.stopPropagation();
          setPickerPos({ x: e.clientX, y: e.clientY });
        }}
      />

      {pickerPos && (
        <ColormapPicker
          current={colormap}
          x={pickerPos.x}
          y={pickerPos.y}
          onClose={() => setPickerPos(null)}
          onSelect={(name) => {
            setColormap(name);
            setPickerPos(null);
          }}
        />
      )}

      {/* Bottom number scrubber */}
      <NumberScrubber
        className="text-white drop-shadow-lg"
        decimals={VISUALIZATION_CONFIG.SCALE_BAR_DECIMALS}
        max={maxValue}
        min={0}
        step={dynamicStep}
        value={minValue}
        onChange={onMinChange}
      />
    </div>
  );
};
