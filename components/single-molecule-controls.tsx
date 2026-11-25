"use client";

import { Button } from "@heroui/button";
import { Slider } from "@heroui/react";
import { Tooltip } from "@heroui/tooltip";
import { useState, useRef, useEffect } from "react";

import { CameraPanel } from "./camera-panel";
import { SingleMoleculeGenePicker } from "./single-molecule-gene-picker";
import { useSliderRangeLocal } from "./slider-range-popover";
import { MoleculeOpacitySlider } from "./molecule-opacity-slider";
import { FadeScroll } from "./fade-scroll";

import {
  usePanelSingleMoleculeVisualizationStore,
  usePanelId,
} from "@/lib/hooks/usePanelStores";
import { useSplitScreenStore } from "@/lib/stores/splitScreenStore";
import { ViewerTopBar } from "@/components/viewer-top-bar";
import { glassButton, glassPanel } from "@/components/primitives";
import { VISUALIZATION_CONFIG } from "@/lib/config/visualization.config";

const BTN = "w-14 h-14 min-w-0 rounded-full font-medium text-xs";

export function SingleMoleculeControls() {
  const { isSplitMode, enableSplit } = useSplitScreenStore();
  const panelId = usePanelId();
  const [isPanelOpen, setIsPanelOpen] = useState(false);

  const {
    globalScale,
    setGlobalScale,
    viewMode,
    setViewMode,
    sceneRotation,
    setSceneRotation,
    flipX,
    setFlipX,
    flipY,
    setFlipY,
    exportBoxEnabled,
    exportBoxWidthMm,
    exportBoxHeightMm,
    setExportBoxEnabled,
    setExportBoxWidthMm,
    setExportBoxHeightMm,
    setExportBoxCenterPx,
  } = usePanelSingleMoleculeVisualizationStore();

  const [isCameraOpen, setIsCameraOpen] = useState(false);
  const controlsRef = useRef<HTMLDivElement>(null);
  const topRightRef = useRef<HTMLDivElement>(null);

  // Right-click the dot-size slider to edit its range (local, like the SC one).
  const dotRange = useSliderRangeLocal(
    VISUALIZATION_CONFIG.SINGLE_MOLECULE_GLOBAL_SCALE_MIN,
    VISUALIZATION_CONFIG.SINGLE_MOLECULE_GLOBAL_SCALE_MAX,
    globalScale,
  );

  // Share: copy the link to this view and show it with its QR code

  // Track shift key for 45° snap on rotation slider
  useEffect(() => {
    const handleKeyDown = (e: KeyboardEvent) => {
      if (e.key === "Shift") (window as any).__shiftHeld = true;
    };
    const handleKeyUp = (e: KeyboardEvent) => {
      if (e.key === "Shift") (window as any).__shiftHeld = false;
    };
    window.addEventListener("keydown", handleKeyDown);
    window.addEventListener("keyup", handleKeyUp);
    return () => {
      window.removeEventListener("keydown", handleKeyDown);
      window.removeEventListener("keyup", handleKeyUp);
    };
  }, []);

  return (
    <>
      {/* Left rail — data controls */}
      <div
        ref={controlsRef}
        className="absolute top-28 left-4 z-50 flex flex-col gap-2"
        data-ui-overlay
      >
        {/* The buttons scroll when the window is too short for them; the gene
            panel sits beside the column, outside what it clips */}
        <FadeScroll
          className="-mx-2 -my-1 flex max-h-[calc(100dvh-13rem)] flex-col gap-2 px-2 py-1"
          data-testid="control-rail"
        >
          {/* Gene Selection Button */}
          <Button
            className={`${BTN} ${isPanelOpen ? "" : glassButton()}`}
            color={isPanelOpen ? "primary" : "default"}
            variant={isPanelOpen ? "shadow" : "light"}
            onPress={() => setIsPanelOpen(!isPanelOpen)}
          >
            Genes
          </Button>

          {/* Dot Size Slider — right-click to change the range */}
          <Tooltip content="Dot size (right-click for range)" placement="right">
            <div
              className={`w-14 h-32 rounded-full border-2 border-default-200 p-2 flex flex-col items-center justify-center ${glassButton()}`}
              onContextMenu={dotRange.onContextMenu}
            >
              <Slider
                aria-label="Dot size"
                className="h-full"
                maxValue={dotRange.max}
                minValue={dotRange.min}
                orientation="vertical"
                size="sm"
                step={VISUALIZATION_CONFIG.SINGLE_MOLECULE_GLOBAL_SCALE_STEP}
                value={globalScale}
                onChange={(value) => setGlobalScale(value as number)}
              />
            </div>
          </Tooltip>
          {dotRange.popover}

          {/* Opacity — dense spots build up below 100% */}
          <MoleculeOpacitySlider />

          {/* Split Screen Button */}
          {!isSplitMode && !panelId && (
            <Tooltip content="Split screen" placement="right">
              <Button
                aria-label="Split screen"
                className={`${BTN} ${glassButton()}`}
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
        </FadeScroll>

        {/* Gene Selection Panel */}
        {isPanelOpen && (
          <div
            className={`absolute top-0 left-16 z-[var(--z-panel)] w-[320px] ${glassPanel()}`}
          >
            <div className="p-4">
              <SingleMoleculeGenePicker />
            </div>
          </div>
        )}
      </div>

      {/* Top-right cluster — view & session controls */}
      <ViewerTopBar
        barRef={topRightRef}
        camera={{
          open: isCameraOpen,
          onToggle: () => {
            setIsCameraOpen(!isCameraOpen);
            if (!isCameraOpen) setIsPanelOpen(false);
          },
        }}
      >
        {isCameraOpen && (
          <CameraPanel
            canExportBox={true}
            controlsRef={topRightRef}
            exportBox={{
              enabled: exportBoxEnabled,
              widthMm: exportBoxWidthMm,
              heightMm: exportBoxHeightMm,
              setEnabled: setExportBoxEnabled,
              setWidthMm: setExportBoxWidthMm,
              setHeightMm: setExportBoxHeightMm,
              resetCenter: () => setExportBoxCenterPx(null),
            }}
            flipX={flipX}
            flipY={flipY}
            placement="top-right"
            viewerKind="sm"
            sceneRotation={sceneRotation}
            setFlipX={setFlipX}
            setFlipY={setFlipY}
            setSceneRotation={setSceneRotation}
            setViewMode={setViewMode}
            viewMode={viewMode}
            onClose={() => setIsCameraOpen(false)}
          />
        )}
      </ViewerTopBar>
    </>
  );
}
