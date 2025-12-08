"use client";

import React from "react";

import { NumberScrubber } from "@/components/ui/number-scrubber";
import { VISUALIZATION_CONFIG } from "@/lib/config/visualization.config";

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

      {/* Gradient bar */}
      <div
        className="w-8 h-32 rounded-lg shadow-lg"
        style={{
          background:
            "linear-gradient(to bottom, #ff0000 0%, #ffffff 50%, #0000ff 100%)",
        }}
      />

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
