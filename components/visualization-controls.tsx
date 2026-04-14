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
