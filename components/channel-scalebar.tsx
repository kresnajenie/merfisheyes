"use client";

import React from "react";

import { NumberScrubber } from "@/components/ui/number-scrubber";
import { VISUALIZATION_CONFIG } from "@/lib/config/visualization.config";

interface ChannelScalebarProps {
  color: string; // "#00FF00" or "#FF00FF" — gradient top
  label?: string;
  minValue: number;
  maxValue: number;
  onMinChange: (v: number) => void;
  onMaxChange: (v: number) => void;
}

/**
 * Minimal scale bar that fades from black (bottom) → channel color (top).
 * Used in two-gene coexpression mode where each gene drives one channel.
 */
export const ChannelScalebar: React.FC<ChannelScalebarProps> = ({
  color,
  label,
  minValue,
  maxValue,
  onMinChange,
  onMaxChange,
}) => {
  const dynamicStep = Math.max(
    VISUALIZATION_CONFIG.SCALE_BAR_MIN_STEP,
    maxValue * VISUALIZATION_CONFIG.SCALE_BAR_STEP_PERCENTAGE,
  );

  return (
    <div className="flex flex-col items-center gap-1">
      <NumberScrubber
        className="text-white drop-shadow-lg text-xs"
        decimals={VISUALIZATION_CONFIG.SCALE_BAR_DECIMALS}
        min={minValue}
        step={dynamicStep}
        value={maxValue}
        onChange={onMaxChange}
      />
      <div
        className="w-6 h-24 rounded-md shadow-lg"
        style={{ background: `linear-gradient(to top, #000000, ${color})` }}
        title={label}
      />
      <NumberScrubber
        className="text-white drop-shadow-lg text-xs"
        decimals={VISUALIZATION_CONFIG.SCALE_BAR_DECIMALS}
        max={maxValue}
        min={0}
        step={dynamicStep}
        value={minValue}
        onChange={onMinChange}
      />
      {label && (
        <span className="text-[10px] text-white/80 drop-shadow font-medium max-w-[80px] truncate">
          {label}
        </span>
      )}
    </div>
  );
};
