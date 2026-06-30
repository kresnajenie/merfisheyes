"use client";

import { Button } from "@heroui/button";
import { Slider } from "@heroui/react";
import { useRef, useEffect, useState } from "react";

import { useSliderRange } from "./slider-range-popover";
import { TourSection } from "./tour-section";

import { usePanelVisualizationStore } from "@/lib/hooks/usePanelStores";
import { glassPanel } from "@/components/primitives";
import { VISUALIZATION_CONFIG } from "@/lib/config/visualization.config";

interface AdvancedVizPanelProps {
  onClose: () => void;
  controlsRef?: React.RefObject<HTMLDivElement>;
}

export function AdvancedVizPanel({
  onClose,
  controlsRef,
}: AdvancedVizPanelProps) {
  const panelRef = useRef<HTMLDivElement>(null);
  const R = VISUALIZATION_CONFIG.ADVANCED_SLIDER_RANGES;
  const clearSliderRanges = usePanelVisualizationStore(
    (s) => s.clearSliderRanges,
  );
  const {
    selectedSizeMultiplier,
    greyedOutSizeMultiplier,
    greyedOutAlpha,
    expressionAlphaMin,
    expressionAlphaMax,
    pointSizeMultiplierMin,
    pointSizeMultiplierMax,
    targetPx,
    orbitSpeed,
    orbitYBob,
    viewMode,
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

  // Local state mirrors store values for responsive slider display
  const [local, setLocal] = useState({
    selectedSizeMultiplier,
    greyedOutSizeMultiplier,
    greyedOutAlpha,
    expressionAlphaMin,
    expressionAlphaMax,
    pointSizeMultiplierMin,
    pointSizeMultiplierMax,
    targetPx,
    orbitSpeed,
    orbitYBob,
  });

  // Sync local state when store changes externally (e.g. URL restore, reset)
  useEffect(() => {
    setLocal({
      selectedSizeMultiplier,
      greyedOutSizeMultiplier,
      greyedOutAlpha,
      expressionAlphaMin,
      expressionAlphaMax,
      pointSizeMultiplierMin,
      pointSizeMultiplierMax,
      targetPx,
      orbitSpeed,
      orbitYBob,
    });
  }, [
    selectedSizeMultiplier,
    greyedOutSizeMultiplier,
    greyedOutAlpha,
    expressionAlphaMin,
    expressionAlphaMax,
    pointSizeMultiplierMin,
    pointSizeMultiplierMax,
    targetPx,
    orbitSpeed,
    orbitYBob,
  ]);

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
        onChange={(v) => setLocal((prev) => ({ ...prev, [key]: v as number }))}
        onChangeEnd={(v) => setAdvancedViz(key, v as number)}
      />
    </div>
  );

  const sizeSliderRow = (
    label: string,
    key:
      | "selectedSizeMultiplier"
      | "greyedOutSizeMultiplier"
      | "pointSizeMultiplierMin"
      | "pointSizeMultiplierMax",
    value: number,
    step: number,
  ) => (
    <SizeSliderRow
      key={key}
      label={label}
      rangeKey={key}
      step={step}
      value={value}
      onChange={(v) => setLocal((prev) => ({ ...prev, [key]: v }))}
      onChangeEnd={(v) => setAdvancedViz(key, v)}
    />
  );

  return (
    <div
      ref={panelRef}
      className={`absolute top-0 left-16 z-[var(--z-panel)] w-[280px] ${glassPanel()}`}
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
            <svg
              className="w-3 h-3"
              fill="none"
              stroke="currentColor"
              strokeWidth={2}
              viewBox="0 0 24 24"
            >
              <path d="M6 18L18 6M6 6l12 12" strokeLinecap="round" />
            </svg>
          </Button>
        </div>

        {/* Celltype Mode */}
        <div>
          <span className="text-xs font-medium text-primary">
            Celltype Mode
          </span>
          {sizeSliderRow(
            "Selected size",
            "selectedSizeMultiplier",
            local.selectedSizeMultiplier,
            R.selectedSizeMultiplier.step,
          )}
          {sizeSliderRow(
            "Unselected size",
            "greyedOutSizeMultiplier",
            local.greyedOutSizeMultiplier,
            R.greyedOutSizeMultiplier.step,
          )}
          {sliderRow(
            "Unselected alpha",
            "greyedOutAlpha",
            local.greyedOutAlpha,
            R.greyedOutAlpha.min,
            R.greyedOutAlpha.max,
            R.greyedOutAlpha.step,
          )}
        </div>

        {/* Gene Expression Mode */}
        <div>
          <span className="text-xs font-medium text-primary">
            Gene Expression
          </span>
          {sizeSliderRow(
            "Size min",
            "pointSizeMultiplierMin",
            local.pointSizeMultiplierMin,
            R.pointSizeMultiplierMin.step,
          )}
          {sizeSliderRow(
            "Size max",
            "pointSizeMultiplierMax",
            local.pointSizeMultiplierMax,
            R.pointSizeMultiplierMax.step,
          )}
          {sliderRow(
            "Alpha min",
            "expressionAlphaMin",
            local.expressionAlphaMin,
            R.expressionAlphaMin.min,
            R.expressionAlphaMin.max,
            R.expressionAlphaMin.step,
          )}
          {sliderRow(
            "Alpha max",
            "expressionAlphaMax",
            local.expressionAlphaMax,
            R.expressionAlphaMax.min,
            R.expressionAlphaMax.max,
            R.expressionAlphaMax.step,
          )}
        </div>

        {/* 3D Orbit (O key) — only relevant in 3D view */}
        {viewMode === "3D" && (
          <div>
            <span className="text-xs font-medium text-primary">
              3D Orbit (O)
            </span>
            {sliderRow("Speed", "orbitSpeed", local.orbitSpeed, 0.25, 4, 0.05)}
            {sliderRow(
              "Vertical bob",
              "orbitYBob",
              local.orbitYBob,
              0,
              1,
              0.05,
            )}
          </div>
        )}

        {/* Tour */}
        <TourSection />

        {/* Reset */}
        <Button
          className="w-full"
          color="danger"
          size="sm"
          variant="flat"
          onPress={() => {
            setAdvancedViz(
              "selectedSizeMultiplier",
              VISUALIZATION_CONFIG.SELECTED_SIZE_MULTIPLIER as number,
            );
            setAdvancedViz(
              "greyedOutSizeMultiplier",
              VISUALIZATION_CONFIG.GREYED_OUT_SIZE_MULTIPLIER as number,
            );
            setAdvancedViz(
              "greyedOutAlpha",
              VISUALIZATION_CONFIG.GREYED_OUT_ALPHA as number,
            );
            setAdvancedViz(
              "expressionAlphaMin",
              VISUALIZATION_CONFIG.EXPRESSION_ALPHA_MIN as number,
            );
            setAdvancedViz(
              "expressionAlphaMax",
              VISUALIZATION_CONFIG.EXPRESSION_ALPHA_MAX as number,
            );
            setAdvancedViz(
              "pointSizeMultiplierMin",
              VISUALIZATION_CONFIG.POINT_SIZE_MULTIPLIER_MIN as number,
            );
            setAdvancedViz(
              "pointSizeMultiplierMax",
              VISUALIZATION_CONFIG.POINT_SIZE_MULTIPLIER_MAX as number,
            );
            setAdvancedViz(
              "targetPx",
              VISUALIZATION_CONFIG.TARGET_PX_DEFAULT as number,
            );
            clearSliderRanges([
              "selectedSizeMultiplier",
              "greyedOutSizeMultiplier",
              "pointSizeMultiplierMin",
              "pointSizeMultiplierMax",
            ]);
          }}
        >
          Reset to Defaults
        </Button>
      </div>
    </div>
  );
}

function SizeSliderRow({
  label,
  rangeKey,
  value,
  step,
  onChange,
  onChangeEnd,
}: {
  label: string;
  rangeKey:
    | "selectedSizeMultiplier"
    | "greyedOutSizeMultiplier"
    | "pointSizeMultiplierMin"
    | "pointSizeMultiplierMax";
  value: number;
  step: number;
  onChange: (v: number) => void;
  onChangeEnd: (v: number) => void;
}) {
  const R = VISUALIZATION_CONFIG.ADVANCED_SLIDER_RANGES;
  const { min, max, onContextMenu, popover } = useSliderRange(
    rangeKey,
    R[rangeKey].min,
    R[rangeKey].max,
    value,
  );

  return (
    <div onContextMenu={onContextMenu}>
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
        onChange={(v) => onChange(v as number)}
        onChangeEnd={(v) => onChangeEnd(v as number)}
      />
      {popover}
    </div>
  );
}
