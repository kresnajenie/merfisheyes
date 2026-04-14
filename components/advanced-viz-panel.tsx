"use client";

import { Button } from "@heroui/button";
import { Slider } from "@heroui/react";
import { useRef, useEffect } from "react";

import { usePanelVisualizationStore } from "@/lib/hooks/usePanelStores";
import { glassButton } from "@/components/primitives";
import { VISUALIZATION_CONFIG } from "@/lib/config/visualization.config";

interface AdvancedVizPanelProps {
  onClose: () => void;
  controlsRef?: React.RefObject<HTMLDivElement>;
}

export function AdvancedVizPanel({ onClose, controlsRef }: AdvancedVizPanelProps) {
  const panelRef = useRef<HTMLDivElement>(null);
  const {
    selectedSizeMultiplier,
    greyedOutSizeMultiplier,
    greyedOutAlpha,
    expressionAlphaMin,
    expressionAlphaMax,
    pointSizeMultiplierMin,
    pointSizeMultiplierMax,
    targetPx,
    setAdvancedViz,
  } = usePanelVisualizationStore();

  // Handle click outside to close panel
  useEffect(() => {
    const handleClickOutside = (event: MouseEvent) => {
      const target = event.target as Node;

      if (panelRef.current && panelRef.current.contains(target)) return;
      if (controlsRef?.current && controlsRef.current.contains(target)) return;

      onClose();
    };

    document.addEventListener("mousedown", handleClickOutside);
    return () => document.removeEventListener("mousedown", handleClickOutside);
  }, [onClose, controlsRef]);

  const sliderRow = (
    label: string,
    key: string,
    value: number,
    min: number,
    max: number,
    step: number,
  ) => (
    <div>
      <div className="flex justify-between text-xs mb-1">
        <span className="text-default-500">{label}</span>
        <span className="text-default-400">{value.toFixed(2)}</span>
      </div>
      <Slider
        aria-label={label}
        className="w-full"
        maxValue={max}
        minValue={min}
        size="sm"
        step={step}
        value={value}
        onChange={(v) => setAdvancedViz(key, v as number)}
      />
    </div>
  );

  return (
    <div
      ref={panelRef}
      className={`absolute top-0 left-16 z-50 w-[280px] border-2 border-white/20 rounded-3xl shadow-lg ${glassButton()}`}
    >
      <div className="p-4 space-y-4">
        {/* Header */}
        <div className="flex justify-between items-center">
          <span className="text-sm font-medium">Advanced Settings</span>
          <Button
            className="h-6 w-6 min-w-0"
            size="sm"
            variant="light"
            onPress={onClose}
          >
            <svg className="w-3 h-3" fill="none" stroke="currentColor" strokeWidth={2} viewBox="0 0 24 24">
              <path d="M6 18L18 6M6 6l12 12" strokeLinecap="round" />
            </svg>
          </Button>
        </div>

        {/* Base Dot Size */}
        <div>
          <span className="text-xs font-medium text-primary">Base Dot Size</span>
          {sliderRow(
            "Initial size (px)",
            "targetPx",
            targetPx,
            VISUALIZATION_CONFIG.TARGET_PX_MIN as number,
            VISUALIZATION_CONFIG.TARGET_PX_MAX as number,
            VISUALIZATION_CONFIG.TARGET_PX_STEP as number,
          )}
        </div>

        {/* Celltype Mode */}
        <div>
          <span className="text-xs font-medium text-primary">Celltype Mode</span>
          {sliderRow("Selected size", "selectedSizeMultiplier", selectedSizeMultiplier, 0.1, 10.0, 0.1)}
          {sliderRow("Unselected size", "greyedOutSizeMultiplier", greyedOutSizeMultiplier, 0.01, 5.0, 0.01)}
          {sliderRow("Unselected alpha", "greyedOutAlpha", greyedOutAlpha, 0.0, 1.0, 0.05)}
        </div>

        {/* Gene Expression Mode */}
        <div>
          <span className="text-xs font-medium text-primary">Gene Expression</span>
          {sliderRow("Size min", "pointSizeMultiplierMin", pointSizeMultiplierMin, 0.1, 5.0, 0.1)}
          {sliderRow("Size max", "pointSizeMultiplierMax", pointSizeMultiplierMax, 0.1, 10.0, 0.1)}
          {sliderRow("Alpha min", "expressionAlphaMin", expressionAlphaMin, 0.0, 1.0, 0.05)}
          {sliderRow("Alpha max", "expressionAlphaMax", expressionAlphaMax, 0.0, 1.0, 0.05)}
        </div>

        {/* Reset */}
        <Button
          className="w-full"
          color="danger"
          size="sm"
          variant="flat"
          onPress={() => {
            setAdvancedViz("selectedSizeMultiplier", VISUALIZATION_CONFIG.SELECTED_SIZE_MULTIPLIER as number);
            setAdvancedViz("greyedOutSizeMultiplier", VISUALIZATION_CONFIG.GREYED_OUT_SIZE_MULTIPLIER as number);
            setAdvancedViz("greyedOutAlpha", VISUALIZATION_CONFIG.GREYED_OUT_ALPHA as number);
            setAdvancedViz("expressionAlphaMin", VISUALIZATION_CONFIG.EXPRESSION_ALPHA_MIN as number);
            setAdvancedViz("expressionAlphaMax", VISUALIZATION_CONFIG.EXPRESSION_ALPHA_MAX as number);
            setAdvancedViz("pointSizeMultiplierMin", VISUALIZATION_CONFIG.POINT_SIZE_MULTIPLIER_MIN as number);
            setAdvancedViz("pointSizeMultiplierMax", VISUALIZATION_CONFIG.POINT_SIZE_MULTIPLIER_MAX as number);
            setAdvancedViz("targetPx", VISUALIZATION_CONFIG.TARGET_PX_DEFAULT as number);
          }}
        >
          Reset to Defaults
        </Button>
      </div>
    </div>
  );
}
