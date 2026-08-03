"use client";

import { Button } from "@heroui/button";
import { Slider } from "@heroui/react";
import { Tooltip } from "@heroui/tooltip";
import { useState, useRef, useEffect } from "react";
import { toast } from "react-toastify";

import { CameraPanel } from "./camera-panel";
import { SingleMoleculeGenePicker } from "./single-molecule-gene-picker";
import { useSliderRangeLocal } from "./slider-range-popover";

import {
  usePanelSingleMoleculeVisualizationStore,
  usePanelId,
} from "@/lib/hooks/usePanelStores";
import { useSplitScreenStore } from "@/lib/stores/splitScreenStore";
import { glassButton, glassPanel } from "@/components/primitives";
import { VISUALIZATION_CONFIG } from "@/lib/config/visualization.config";

const BTN = "w-14 h-14 min-w-0 rounded-full font-medium text-xs";

export function SingleMoleculeControls() {
  const { isSplitMode, enableSplit } = useSplitScreenStore();
  const setHideUi = useSplitScreenStore((s) => s.setHideUi);
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

  const handleShare = async () => {
    try {
      await navigator.clipboard.writeText(window.location.href);
      toast.success("Link copied");
    } catch {
      toast.error("Couldn't copy link");
    }
  };

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

        {/* Split Screen Button */}
        {!isSplitMode && !panelId && (
          <Tooltip content="Split screen" placement="right">
            <Button
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
      <div
        ref={topRightRef}
        data-ui-overlay
        className="absolute top-4 right-4 z-[var(--z-rail)] flex flex-row gap-2"
      >
        {/* Camera Button */}
        <Tooltip content="Camera controls" placement="bottom">
          <Button
            className={`${BTN} ${isCameraOpen ? "" : glassButton()}`}
            color={isCameraOpen ? "primary" : "default"}
            variant={isCameraOpen ? "shadow" : "light"}
            onPress={() => {
              setIsCameraOpen(!isCameraOpen);
              if (!isCameraOpen) setIsPanelOpen(false);
            }}
          >
            <svg className="w-5 h-5" fill="none" stroke="currentColor" strokeWidth={1.5} viewBox="0 0 24 24">
              <path d="M6.827 6.175A2.31 2.31 0 015.186 7.23c-.38.054-.757.112-1.134.175C2.999 7.58 2.25 8.507 2.25 9.574V18a2.25 2.25 0 002.25 2.25h15A2.25 2.25 0 0021.75 18V9.574c0-1.067-.75-1.994-1.802-2.169a47.865 47.865 0 00-1.134-.175 2.31 2.31 0 01-1.64-1.055l-.822-1.316a2.192 2.192 0 00-1.736-1.039 48.774 48.774 0 00-5.232 0 2.192 2.192 0 00-1.736 1.039l-.821 1.316z" strokeLinecap="round" strokeLinejoin="round" />
              <path d="M16.5 12.75a4.5 4.5 0 11-9 0 4.5 4.5 0 019 0z" strokeLinecap="round" strokeLinejoin="round" />
            </svg>
          </Button>
        </Tooltip>

        {/* Hide UI Button — strips overlays for screenshotting (H key) */}
        <Tooltip content="Hide UI for screenshot (H)" placement="bottom">
          <Button
            className={`${BTN} ${glassButton()}`}
            color="default"
            variant="light"
            onPress={() => setHideUi(true)}
          >
            <svg className="w-5 h-5" fill="none" stroke="currentColor" strokeWidth={1.5} viewBox="0 0 24 24">
              <path d="M3.98 8.223A10.477 10.477 0 001.934 12C3.226 16.338 7.244 19.5 12 19.5c.993 0 1.953-.138 2.863-.395M6.228 6.228A10.45 10.45 0 0112 4.5c4.756 0 8.773 3.162 10.065 7.498a10.523 10.523 0 01-4.293 5.774M6.228 6.228L3 3m3.228 3.228l3.65 3.65m7.894 7.894L21 21m-3.228-3.228l-3.65-3.65m0 0a3 3 0 10-4.243-4.243m4.242 4.242L9.88 9.88" strokeLinecap="round" strokeLinejoin="round" />
            </svg>
          </Button>
        </Tooltip>

        {/* Share Button — copies the current view URL to the clipboard */}
        <Tooltip content="Copy link to this view" placement="bottom">
          <Button
            className={`${BTN} ${glassButton()}`}
            color="default"
            variant="light"
            onPress={handleShare}
          >
            <svg className="w-5 h-5" fill="none" stroke="currentColor" strokeWidth={1.5} viewBox="0 0 24 24">
              <path d="M7.217 10.907a2.25 2.25 0 100 2.186m0-2.186c.18.324.283.696.283 1.093s-.103.77-.283 1.093m0-2.186l9.566-5.314m-9.566 7.5l9.566 5.314m0 0a2.25 2.25 0 103.935 2.186 2.25 2.25 0 00-3.935-2.186zm0-12.814a2.25 2.25 0 103.933-2.185 2.25 2.25 0 00-3.933 2.185z" strokeLinecap="round" strokeLinejoin="round" />
            </svg>
          </Button>
        </Tooltip>

        {/* Camera Panel */}
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
      </div>
    </>
  );
}
