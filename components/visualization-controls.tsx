"use client";

import type { VisualizationMode } from "@/lib/stores/createVisualizationStore";

import { Button } from "@heroui/button";
import { Slider } from "@heroui/react";
import { Tooltip } from "@heroui/tooltip";
import { useState, useRef, useEffect, useCallback, useMemo } from "react";
import type { StandardizedDataset } from "@/lib/StandardizedDataset";
import { getClusterValue } from "@/lib/StandardizedDataset";
import { getEffectiveColumnType } from "@/lib/utils/column-type-utils";

import { VisualizationPanel } from "./visualization-panel";
import { AdvancedVizPanel } from "./advanced-viz-panel";

import {
  usePanelVisualizationStore,
  usePanelDatasetStore,
  usePanelId,
} from "@/lib/hooks/usePanelStores";
import { useSplitScreenStore } from "@/lib/stores/splitScreenStore";
import { useVisualizationStore } from "@/lib/stores/visualizationStore";
import { glassButton } from "@/components/primitives";
import { VISUALIZATION_CONFIG } from "@/lib/config/visualization.config";

export function VisualizationControls() {
  const {
    panelMode, setPanelMode, sizeScale, setSizeScale, viewMode, setViewMode,
    celltypePlayback, setCelltypePlayback, celltypePlaybackInterval,
    celltypePlaybackSequence,
    selectedColumn, selectedCelltypes, setCelltypes, clusterVersion, columnTypeOverrides,
    sceneRotation, setSceneRotation, flipX, setFlipX, flipY, setFlipY,
  } = usePanelVisualizationStore();
  const { isSplitMode, enableSplit } = useSplitScreenStore();
  const panelId = usePanelId();

  const dataset = usePanelDatasetStore((s) => {
    const id = s.currentDatasetId;
    const ds = id ? s.datasets.get(id) : null;
    return ds && "spatial" in ds ? (ds as StandardizedDataset) : null;
  });
  const is3DDataset = dataset?.spatial?.dimensions === 3;

  // Celltype playback timer — lives here because this component doesn't unmount
  const playIndexRef = useRef(0);

  const celltypeItems = useMemo(() => {
    if (!dataset || !selectedColumn) return [];
    const cluster = dataset.clusters?.find((c) => c.column === selectedColumn);
    if (!cluster) return [];
    if (getEffectiveColumnType(selectedColumn, dataset, columnTypeOverrides) === "numerical") return [];
    const vals = cluster.uniqueValues ?? [...new Set(cluster.values.map(String))].sort();
    return vals;
  }, [dataset, selectedColumn, clusterVersion, columnTypeOverrides]);

  // Build the effective playback list: custom sequence (filtered to valid) or all items
  const playbackList = useMemo(() => {
    if (celltypePlaybackSequence.length === 0) return celltypeItems;

    const validSet = new Set(celltypeItems);
    const filtered = celltypePlaybackSequence.filter((name) => {
      if (validSet.has(name)) return true;
      console.warn(`Playback: skipping unknown celltype "${name}"`);
      return false;
    });

    return filtered.length > 0 ? filtered : celltypeItems;
  }, [celltypeItems, celltypePlaybackSequence]);

  const playbackListRef = useRef(playbackList);
  playbackListRef.current = playbackList;

  useEffect(() => {
    if (!celltypePlayback || playbackList.length === 0) return;

    playIndexRef.current = 0;
    setCelltypes(new Set([playbackList[0]]));

    const timer = setInterval(() => {
      const items = playbackListRef.current;
      if (items.length === 0) return;
      playIndexRef.current = (playIndexRef.current + 1) % items.length;
      setCelltypes(new Set([items[playIndexRef.current]]));
    }, celltypePlaybackInterval * 1000);

    return () => clearInterval(timer);
  }, [celltypePlayback, celltypePlaybackInterval, playbackList.length, setCelltypes]);
  const [isPanelOpen, setIsPanelOpen] = useState(false);
  const [isAdvancedOpen, setIsAdvancedOpen] = useState(false);
  const controlsRef = useRef<HTMLDivElement>(null);

  const handleModeChange = (newMode: VisualizationMode) => {
    setIsAdvancedOpen(false);
    if (panelMode === newMode) {
      setIsPanelOpen(!isPanelOpen);
    } else {
      setPanelMode(newMode);
      setIsPanelOpen(true);
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

  const buttonBaseClass = "w-14 h-14 min-w-0 rounded-full font-medium text-xs";

  // Read playback from global store (works for both left and right panels)
  const globalPlayback = useVisualizationStore((s) => s.celltypePlayback);
  const globalSetPlayback = useVisualizationStore((s) => s.setCelltypePlayback);

  // During playback, show only a floating pause button (left panel) or hide entirely (right panel)
  if (celltypePlayback || globalPlayback) {
    // Right panel: hide entirely during playback
    if (panelId) return null;

    // Left panel: show floating pause button
    return (
      <div className="absolute top-4 left-4 z-50">
        <Tooltip content="Stop playback" placement="right">
          <Button
            className="w-12 h-12 min-w-0 rounded-full"
            color="danger"
            variant="shadow"
            onPress={() => {
              setCelltypePlayback(false);
              globalSetPlayback(false);
            }}
          >
            <svg className="w-5 h-5" fill="currentColor" viewBox="0 0 24 24">
              <rect x="6" y="4" width="4" height="16" />
              <rect x="14" y="4" width="4" height="16" />
            </svg>
          </Button>
        </Tooltip>
      </div>
    );
  }

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

      {/* Rotation Slider */}
      <Tooltip content={`Rotate: ${sceneRotation}°`} placement="right">
        <div
          className={`w-14 h-32 rounded-full border-2 border-default-200 p-2 flex flex-col items-center justify-center ${glassButton()}`}
        >
          <span className="text-[8px] text-default-400 mb-1">↻</span>
          <Slider
            aria-label="Rotation"
            className="h-full"
            maxValue={360}
            minValue={0}
            orientation="vertical"
            size="sm"
            step={1}
            value={sceneRotation}
            onChange={(value) => {
              let v = value as number;
              if ((window as any).__shiftHeld) {
                v = Math.round(v / 45) * 45;
              }
              setSceneRotation(v);
            }}
          />
        </div>
      </Tooltip>

      {/* Flip Buttons */}
      <div className="flex flex-col gap-1">
        <Tooltip content="Flip horizontal" placement="right">
          <Button
            className={`w-14 h-8 min-w-0 rounded-full text-xs ${flipX ? "" : glassButton()}`}
            color={flipX ? "primary" : "default"}
            variant={flipX ? "solid" : "light"}
            onPress={() => setFlipX(!flipX)}
          >
            ⇔
          </Button>
        </Tooltip>
        <Tooltip content="Flip vertical" placement="right">
          <Button
            className={`w-14 h-8 min-w-0 rounded-full text-xs ${flipY ? "" : glassButton()}`}
            color={flipY ? "primary" : "default"}
            variant={flipY ? "solid" : "light"}
            onPress={() => setFlipY(!flipY)}
          >
            ⇕
          </Button>
        </Tooltip>
      </div>

      {/* 2D/3D View Toggle — only for 3D datasets */}
      {is3DDataset && (
        <Tooltip content={viewMode === "2D" ? "Switch to 3D" : "Switch to 2D"} placement="right">
          <Button
            className={`${buttonBaseClass} ${glassButton()}`}
            color="default"
            variant="light"
            onPress={() => setViewMode(viewMode === "2D" ? "3D" : "2D")}
          >
            {viewMode}
          </Button>
        </Tooltip>
      )}

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

      {/* Advanced Settings Button */}
      <Tooltip content="Advanced settings" placement="right">
        <Button
          className={`${buttonBaseClass} ${isAdvancedOpen ? "" : glassButton()}`}
          color={isAdvancedOpen ? "primary" : "default"}
          variant={isAdvancedOpen ? "shadow" : "light"}
          onPress={() => {
            setIsAdvancedOpen((prev) => {
              if (!prev) setIsPanelOpen(false);
              return !prev;
            });
          }}
        >
          <svg className="w-5 h-5" fill="none" stroke="currentColor" strokeWidth={1.5} viewBox="0 0 24 24">
            <path d="M9.594 3.94c.09-.542.56-.94 1.11-.94h2.593c.55 0 1.02.398 1.11.94l.213 1.281c.063.374.313.686.645.87.074.04.147.083.22.127.325.196.72.257 1.075.124l1.217-.456a1.125 1.125 0 011.37.49l1.296 2.247a1.125 1.125 0 01-.26 1.431l-1.003.827c-.293.241-.438.613-.431.992a6.759 6.759 0 010 .255c-.007.378.138.75.43.991l1.004.827c.424.35.534.954.26 1.43l-1.298 2.247a1.125 1.125 0 01-1.369.491l-1.217-.456c-.355-.133-.75-.072-1.076.124a6.57 6.57 0 01-.22.128c-.331.183-.581.495-.644.869l-.213 1.28c-.09.543-.56.941-1.11.941h-2.594c-.55 0-1.02-.398-1.11-.94l-.213-1.281c-.062-.374-.312-.686-.644-.87a6.52 6.52 0 01-.22-.127c-.325-.196-.72-.257-1.076-.124l-1.217.456a1.125 1.125 0 01-1.369-.49l-1.297-2.247a1.125 1.125 0 01.26-1.431l1.004-.827c.292-.24.437-.613.43-.991a6.932 6.932 0 010-.255c.007-.378-.138-.75-.43-.992l-1.004-.827a1.125 1.125 0 01-.26-1.43l1.297-2.247a1.125 1.125 0 011.37-.491l1.216.456c.356.133.751.072 1.076-.124.072-.044.146-.087.22-.128.332-.183.582-.495.644-.869l.214-1.281z" strokeLinecap="round" />
            <path d="M15 12a3 3 0 11-6 0 3 3 0 016 0z" strokeLinecap="round" />
          </svg>
        </Button>
      </Tooltip>

      {/* Visualization Panel */}
      {isPanelOpen && (
        <VisualizationPanel
          controlsRef={controlsRef}
          mode={panelMode}
          onClose={() => setIsPanelOpen(false)}
        />
      )}

      {/* Advanced Visualization Panel */}
      {isAdvancedOpen && (
        <AdvancedVizPanel
          controlsRef={controlsRef}
          onClose={() => setIsAdvancedOpen(false)}
        />
      )}
    </div>
  );
}
