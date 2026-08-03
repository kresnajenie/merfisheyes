"use client";

import type { VisualizationMode } from "@/lib/stores/createVisualizationStore";
import type { StandardizedDataset } from "@/lib/StandardizedDataset";

import { Button } from "@heroui/button";
import { Slider } from "@heroui/react";
import { Tooltip } from "@heroui/tooltip";
import { useState, useRef, useEffect, useMemo } from "react";
import { toast } from "react-toastify";

import { VisualizationPanel } from "./visualization-panel";
import { AdvancedVizPanel } from "./advanced-viz-panel";
import { CameraPanel } from "./camera-panel";
import { ScSmModeToggle, type ScSmMode } from "./sc-sm-mode-toggle";
import { useSliderRange } from "./slider-range-popover";
import { PlotPanel } from "./plot-panel";
import { DegPanel } from "./deg-panel";

import { getEffectiveColumnType } from "@/lib/utils/column-type-utils";
import {
  ensureDeStatsForColumn,
  isDeStatsInFlight,
} from "@/lib/utils/de-stats";
import {
  usePanelVisualizationStore,
  usePanelDatasetStore,
  usePanelId,
  usePanelSingleMoleculeVisualizationStore,
} from "@/lib/hooks/usePanelStores";
import { useSplitScreenStore } from "@/lib/stores/splitScreenStore";
import { useVisualizationStore } from "@/lib/stores/visualizationStore";
import { useSingleMoleculeVisualizationStore } from "@/lib/stores/singleMoleculeVisualizationStore";
import { useSingleMoleculeStore } from "@/lib/stores/singleMoleculeStore";
import { glassButton } from "@/components/primitives";
import { VISUALIZATION_CONFIG } from "@/lib/config/visualization.config";

export function VisualizationControls() {
  const {
    panelMode,
    setPanelMode,
    sizeScale,
    setSizeScale,
    viewMode,
    setViewMode,
    celltypePlayback,
    setCelltypePlayback,
    celltypePlaybackInterval,
    celltypePlaybackSequence,
    selectedColumn,
    selectedCelltypes,
    setCelltypes,
    clusterVersion,
    columnTypeOverrides,
    sceneRotation,
    setSceneRotation,
    flipX,
    setFlipX,
    flipY,
    setFlipY,
    plotPanelOpen,
    setPlotPanelOpen,
    exportBoxEnabled,
    exportBoxWidthMm,
    exportBoxHeightMm,
    setExportBoxEnabled,
    setExportBoxWidthMm,
    setExportBoxHeightMm,
    setExportBoxCenterPx,
    resetCamera,
    targetFilterEnabled,
    setTargetFilterEnabled,
    targetFilterRadius,
    setTargetFilterRadius,
    scLayerVisible,
    setScLayerVisible,
  } = usePanelVisualizationStore();
  // SM layer visibility reads the global SM store to match three-scene, which
  // renders the SM overlay from the global (non-panel) store.
  const smLayerVisible = useSingleMoleculeVisualizationStore(
    (s) => s.smLayerVisible,
  );
  const setSmLayerVisible = useSingleMoleculeVisualizationStore(
    (s) => s.setSmLayerVisible,
  );
  // An SM overlay is only present when a mapping.json with linkColumn
  // "__all__" has loaded an SM dataset into the global SM store. The SC/SM
  // layer toggles + SC↔SM mode toggles only make sense in that case.
  const hasSmOverlay = useSingleMoleculeStore((s) => s.currentDatasetId != null);
  const { isSplitMode, enableSplit } = useSplitScreenStore();
  const panelId = usePanelId();

  const dataset = usePanelDatasetStore((s) => {
    const id = s.currentDatasetId;
    const ds = id ? s.datasets.get(id) : null;

    return ds && "spatial" in ds ? (ds as StandardizedDataset) : null;
  });
  // is3DDataset check moved to CameraPanel

  // Celltype playback timer — lives here because this component doesn't unmount
  const playIndexRef = useRef(0);

  const celltypeItems = useMemo(() => {
    if (!dataset || !selectedColumn) return [];
    const cluster = dataset.clusters?.find((c) => c.column === selectedColumn);

    if (!cluster) return [];
    if (
      getEffectiveColumnType(selectedColumn, dataset, columnTypeOverrides) ===
      "numerical"
    )
      return [];
    const vals =
      cluster.uniqueValues ?? [...new Set(cluster.values.map(String))].sort();

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
  }, [
    celltypePlayback,
    celltypePlaybackInterval,
    playbackList.length,
    setCelltypes,
  ]);
  const [isPanelOpen, setIsPanelOpen] = useState(false);
  const [isAdvancedOpen, setIsAdvancedOpen] = useState(false);
  const [isCameraOpen, setIsCameraOpen] = useState(false);
  // Shared SC/SM mode for the Gene panel tab toggle AND the size-slider toggle.
  // TODO: gate visibility on an SM-overlay-loaded signal.
  const [scSmMode, setScSmMode] = useState<ScSmMode>("sc");
  // DEG panel open state lives in the store so it can be opened from
  // elsewhere (e.g. right-clicking a celltype badge in the legends).
  const isDegOpen = usePanelVisualizationStore((s) => s.degPanelOpen);
  const setIsDegOpen = usePanelVisualizationStore((s) => s.setDegPanelOpen);
  const setHideUi = useSplitScreenStore((s) => s.setHideUi);
  const controlsRef = useRef<HTMLDivElement>(null);
  // Top-right cluster (Camera / Hide UI / Settings / Share). Its own ref so the
  // Camera and Advanced flyouts anchor here and their click-outside ignores it.
  const topRightRef = useRef<HTMLDivElement>(null);

  const handleShare = async () => {
    try {
      await navigator.clipboard.writeText(window.location.href);
      toast.success("Link copied");
    } catch {
      toast.error("Couldn't copy link");
    }
  };

  const hasDeStats =
    !!dataset?.deStats || (dataset?.availableDeStatsColumns?.length ?? 0) > 0;

  // Keep the DEG target/reference in sync with the celltype selection.
  // Runs even when the panel is closed so opening it lands on a sensible
  // default. In Auto mode the target tracks the most-recently-selected
  // celltype and the reference (when its own Auto is on) tracks the 2nd
  // most-recently-selected one. In Manual mode each is only corrected when
  // it becomes invalid for the active column.
  const degTarget = usePanelVisualizationStore((s) => s.degTarget);
  const setDegTarget = usePanelVisualizationStore((s) => s.setDegTarget);
  const degReference = usePanelVisualizationStore((s) => s.degReference);
  const setDegReference = usePanelVisualizationStore((s) => s.setDegReference);
  const degTargetAuto = usePanelVisualizationStore((s) => s.degTargetAuto);
  const degReferenceAuto = usePanelVisualizationStore(
    (s) => s.degReferenceAuto,
  );
  const deStatsVersion = usePanelVisualizationStore((s) => s.deStatsVersion);

  useEffect(() => {
    if (!dataset) return;
    const activeCol =
      selectedColumn && dataset.deStatsByColumn.has(selectedColumn)
        ? selectedColumn
        : (dataset.deStats?.column ?? null);
    const deStats = activeCol
      ? (dataset.deStatsByColumn.get(activeCol) ?? null)
      : null;

    if (!deStats || deStats.celltypes.length === 0) return;

    const valid = new Set(deStats.celltypes);
    // Selected celltypes that exist in this column's deStats, in selection
    // order — the last entry is the most recently selected.
    const picks = [...selectedCelltypes].filter((ct) => valid.has(ct));
    const mostRecent = picks[picks.length - 1];
    const secondRecent = picks[picks.length - 2];

    // Target.
    let effectiveTarget: string;

    if (degTargetAuto) {
      effectiveTarget = mostRecent ?? deStats.celltypes[0];
    } else {
      effectiveTarget =
        degTarget && valid.has(degTarget) ? degTarget : deStats.celltypes[0];
    }
    if (effectiveTarget !== degTarget) setDegTarget(effectiveTarget);

    // Reference (null = vs Rest).
    if (degReferenceAuto) {
      const ref =
        secondRecent && secondRecent !== effectiveTarget ? secondRecent : null;

      if (ref !== degReference) setDegReference(ref);
    } else if (
      degReference &&
      (!valid.has(degReference) || degReference === effectiveTarget)
    ) {
      setDegReference(null);
    }
  }, [
    dataset,
    selectedColumn,
    deStatsVersion,
    selectedCelltypes,
    degTarget,
    setDegTarget,
    degReference,
    setDegReference,
    degTargetAuto,
    degReferenceAuto,
  ]);

  // Recompute deStats when the user changes the cluster column. Numerical
  // columns skip (the panel stays on its last categorical column). Cache hits
  // return synchronously inside ensureDeStatsForColumn. For chunked datasets
  // the adapter fetches the precomputed file from disk/S3 instead.
  const incrementDeStatsVersion = usePanelVisualizationStore(
    (s) => s.incrementDeStatsVersion,
  );

  useEffect(() => {
    if (!dataset || !selectedColumn) return;
    if (!hasDeStats) return; // dataset doesn't support DEG (e.g. Xenium)
    if (dataset.deStatsByColumn.has(selectedColumn)) return;
    if (isDeStatsInFlight(dataset.id, selectedColumn)) return;

    const canFetchFromAdapter =
      !!dataset.adapter?.loadDeStats &&
      dataset.availableDeStatsColumns?.includes(selectedColumn);

    if (!canFetchFromAdapter) {
      // Must recompute — requires the cluster to be loaded and matrix in memory.
      const cluster = dataset.clusters?.find(
        (c) => c.column === selectedColumn,
      );

      if (!cluster || cluster.type !== "categorical") return;
      if (!dataset.matrix) return;
    }

    const verb = canFetchFromAdapter ? "Loading" : "Computing";
    const toastId = toast.loading(`${verb} DEG for "${selectedColumn}"...`);
    let lastPct = -1;

    ensureDeStatsForColumn(dataset, selectedColumn, async (frac) => {
      const pct = Math.round(frac * 100);

      if (pct === lastPct) return;
      lastPct = pct;
      toast.update(toastId, {
        render: `${verb} DEG for "${selectedColumn}"... ${pct}%`,
      });
    })
      .then((result) => {
        if (result) {
          toast.update(toastId, {
            render: `DEG ready for "${selectedColumn}"`,
            type: "success",
            isLoading: false,
            autoClose: 1500,
          });
          incrementDeStatsVersion();
        } else {
          toast.dismiss(toastId);
        }
      })
      .catch((err) => {
        console.warn("[deStats] load/compute failed:", err);
        toast.update(toastId, {
          render: `DEG failed for "${selectedColumn}"`,
          type: "error",
          isLoading: false,
          autoClose: 3000,
        });
      });
  }, [
    dataset,
    selectedColumn,
    clusterVersion,
    incrementDeStatsVersion,
    hasDeStats,
  ]);

  const handleModeChange = (newMode: VisualizationMode) => {
    setIsAdvancedOpen(false);
    setIsCameraOpen(false);
    setIsDegOpen(false);
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
      <div className="absolute top-4 left-4 z-[var(--z-rail)]">
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
              <rect height="16" width="4" x="6" y="4" />
              <rect height="16" width="4" x="14" y="4" />
            </svg>
          </Button>
        </Tooltip>
      </div>
    );
  }

  return (
    <>
    <div
      ref={controlsRef}
      data-ui-overlay
      className="absolute top-28 left-4 z-[var(--z-rail)] flex flex-col gap-2"
    >
      {/* Celltype Button */}
      <Tooltip content="Color by cell type / cluster" placement="right">
        <Button
          className={`${buttonBaseClass} ${isPanelOpen && panelMode === "celltype" ? "" : glassButton()}`}
          color={isPanelOpen && panelMode === "celltype" ? "primary" : "default"}
          variant={isPanelOpen && panelMode === "celltype" ? "shadow" : "light"}
          onPress={() => handleModeChange("celltype")}
        >
          Celltype
        </Button>
      </Tooltip>

      {/* Gene Button */}
      <Tooltip content="Color by gene expression" placement="right">
        <Button
          className={`${buttonBaseClass} ${isPanelOpen && panelMode === "gene" ? "" : glassButton()}`}
          color={isPanelOpen && panelMode === "gene" ? "primary" : "default"}
          variant={isPanelOpen && panelMode === "gene" ? "shadow" : "light"}
          onPress={() => handleModeChange("gene")}
        >
          Gene
        </Button>
      </Tooltip>

      {/* DEG Button — only when deStats is available on the dataset */}
      {hasDeStats && (
        <Tooltip content="Differentially expressed genes" placement="right">
          <Button
            className={`${buttonBaseClass} ${isDegOpen ? "" : glassButton()}`}
            color={isDegOpen ? "primary" : "default"}
            data-testid="deg-button"
            variant={isDegOpen ? "shadow" : "light"}
            onPress={() => {
              const next = !isDegOpen;

              if (next) {
                setIsPanelOpen(false);
                setIsAdvancedOpen(false);
                setIsCameraOpen(false);
              }
              setIsDegOpen(next);
            }}
          >
            DEG
          </Button>
        </Tooltip>
      )}

      {/* Plot Panel Button */}
      <Tooltip content="Plot panel" placement="right">
        <Button
          className={`${buttonBaseClass} ${plotPanelOpen ? "" : glassButton()}`}
          color={plotPanelOpen ? "primary" : "default"}
          variant={plotPanelOpen ? "shadow" : "light"}
          onPress={() => setPlotPanelOpen(!plotPanelOpen)}
        >
          <svg
            className="w-5 h-5"
            fill="none"
            stroke="currentColor"
            strokeWidth={1.5}
            viewBox="0 0 24 24"
          >
            <path d="M3 3v18h18" strokeLinecap="round" strokeLinejoin="round" />
            <rect height="7" strokeLinejoin="round" width="3" x="6" y="11" />
            <rect height="11" strokeLinejoin="round" width="3" x="11" y="7" />
            <rect height="5" strokeLinejoin="round" width="3" x="16" y="13" />
          </svg>
        </Button>
      </Tooltip>

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

      {/* Dot Size Slider — switches between SC sizeScale and SM globalScale. */}
      <div className="flex flex-row items-center gap-1">
        <DotSizeSlider
          scSizeScale={sizeScale}
          scSmMode={hasSmOverlay ? scSmMode : "sc"}
          setScSizeScale={setSizeScale}
        />
        {/* SC↔SM toggle only when an SM overlay is present. */}
        {hasSmOverlay && (
          <ScSmModeToggle
            mode={scSmMode}
            variant="vertical"
            onChange={setScSmMode}
          />
        )}
      </div>

      {/* Layer visibility — show/hide the entire single-cell or single-molecule
          layer. Only meaningful when an SM overlay is present, so both toggles
          are hidden otherwise (a plain SC dataset has nothing to layer). */}
      {hasSmOverlay && (
        <>
          <LayerVisibilityButton
            buttonBaseClass={buttonBaseClass}
            kind="single cell"
            label="SC"
            visible={scLayerVisible}
            onToggle={() => setScLayerVisible(!scLayerVisible)}
          />
          <LayerVisibilityButton
            buttonBaseClass={buttonBaseClass}
            kind="single molecule"
            label="SM"
            visible={smLayerVisible}
            onToggle={() => setSmLayerVisible(!smLayerVisible)}
          />
        </>
      )}

      {/* Visualization Panel */}
      {isPanelOpen && (
        <VisualizationPanel
          controlsRef={controlsRef}
          mode={panelMode}
          scSmMode={hasSmOverlay ? scSmMode : undefined}
          setScSmMode={hasSmOverlay ? setScSmMode : undefined}
          onClose={() => setIsPanelOpen(false)}
        />
      )}

      {/* Plot Panel (floating, draggable, resizable) */}
      {plotPanelOpen && <PlotPanel />}

      {/* DEG Panel */}
      {isDegOpen && (
        <DegPanel
          controlsRef={controlsRef}
          onClose={() => setIsDegOpen(false)}
        />
      )}
    </div>

    {/* Top-right cluster — view & session controls (Camera / Hide UI /
        Settings / Share). Icon-only row; the Camera and Advanced flyouts open
        below-right, anchored to this container. */}
    <div
      ref={topRightRef}
      data-ui-overlay
      className="absolute top-4 right-4 z-[var(--z-rail)] flex flex-row gap-2"
    >
      {/* Camera Button */}
      <Tooltip content="Camera controls" placement="bottom">
        <Button
          className={`${buttonBaseClass} ${isCameraOpen ? "" : glassButton()}`}
          color={isCameraOpen ? "primary" : "default"}
          variant={isCameraOpen ? "shadow" : "light"}
          onPress={() => {
            setIsCameraOpen((prev) => {
              if (!prev) {
                setIsPanelOpen(false);
                setIsAdvancedOpen(false);
                setIsDegOpen(false);
              }

              return !prev;
            });
          }}
        >
          <svg
            className="w-5 h-5"
            fill="none"
            stroke="currentColor"
            strokeWidth={1.5}
            viewBox="0 0 24 24"
          >
            <path
              d="M6.827 6.175A2.31 2.31 0 015.186 7.23c-.38.054-.757.112-1.134.175C2.999 7.58 2.25 8.507 2.25 9.574V18a2.25 2.25 0 002.25 2.25h15A2.25 2.25 0 0021.75 18V9.574c0-1.067-.75-1.994-1.802-2.169a47.865 47.865 0 00-1.134-.175 2.31 2.31 0 01-1.64-1.055l-.822-1.316a2.192 2.192 0 00-1.736-1.039 48.774 48.774 0 00-5.232 0 2.192 2.192 0 00-1.736 1.039l-.821 1.316z"
              strokeLinecap="round"
              strokeLinejoin="round"
            />
            <path
              d="M16.5 12.75a4.5 4.5 0 11-9 0 4.5 4.5 0 019 0z"
              strokeLinecap="round"
              strokeLinejoin="round"
            />
          </svg>
        </Button>
      </Tooltip>

      {/* Hide UI Button — strips overlays for screenshotting (H key) */}
      <Tooltip content="Hide UI for screenshot (H)" placement="bottom">
        <Button
          className={`${buttonBaseClass} ${glassButton()}`}
          color="default"
          variant="light"
          onPress={() => setHideUi(true)}
        >
          <svg
            className="w-5 h-5"
            fill="none"
            stroke="currentColor"
            strokeWidth={1.5}
            viewBox="0 0 24 24"
          >
            <path
              d="M3.98 8.223A10.477 10.477 0 001.934 12C3.226 16.338 7.244 19.5 12 19.5c.993 0 1.953-.138 2.863-.395M6.228 6.228A10.45 10.45 0 0112 4.5c4.756 0 8.773 3.162 10.065 7.498a10.523 10.523 0 01-4.293 5.774M6.228 6.228L3 3m3.228 3.228l3.65 3.65m7.894 7.894L21 21m-3.228-3.228l-3.65-3.65m0 0a3 3 0 10-4.243-4.243m4.242 4.242L9.88 9.88"
              strokeLinecap="round"
              strokeLinejoin="round"
            />
          </svg>
        </Button>
      </Tooltip>

      {/* Advanced Settings Button */}
      <Tooltip content="Advanced settings" placement="bottom">
        <Button
          className={`${buttonBaseClass} ${isAdvancedOpen ? "" : glassButton()}`}
          color={isAdvancedOpen ? "primary" : "default"}
          variant={isAdvancedOpen ? "shadow" : "light"}
          onPress={() => {
            setIsAdvancedOpen((prev) => {
              if (!prev) {
                setIsPanelOpen(false);
                setIsCameraOpen(false);
                setIsDegOpen(false);
              }

              return !prev;
            });
          }}
        >
          <svg
            className="w-5 h-5"
            fill="none"
            stroke="currentColor"
            strokeWidth={1.5}
            viewBox="0 0 24 24"
          >
            <path
              d="M9.594 3.94c.09-.542.56-.94 1.11-.94h2.593c.55 0 1.02.398 1.11.94l.213 1.281c.063.374.313.686.645.87.074.04.147.083.22.127.325.196.72.257 1.075.124l1.217-.456a1.125 1.125 0 011.37.49l1.296 2.247a1.125 1.125 0 01-.26 1.431l-1.003.827c-.293.241-.438.613-.431.992a6.759 6.759 0 010 .255c-.007.378.138.75.43.991l1.004.827c.424.35.534.954.26 1.43l-1.298 2.247a1.125 1.125 0 01-1.369.491l-1.217-.456c-.355-.133-.75-.072-1.076.124a6.57 6.57 0 01-.22.128c-.331.183-.581.495-.644.869l-.213 1.28c-.09.543-.56.941-1.11.941h-2.594c-.55 0-1.02-.398-1.11-.94l-.213-1.281c-.062-.374-.312-.686-.644-.87a6.52 6.52 0 01-.22-.127c-.325-.196-.72-.257-1.076-.124l-1.217.456a1.125 1.125 0 01-1.369-.49l-1.297-2.247a1.125 1.125 0 01.26-1.431l1.004-.827c.292-.24.437-.613.43-.991a6.932 6.932 0 010-.255c.007-.378-.138-.75-.43-.992l-1.004-.827a1.125 1.125 0 01-.26-1.43l1.297-2.247a1.125 1.125 0 011.37-.491l1.216.456c.356.133.751.072 1.076-.124.072-.044.146-.087.22-.128.332-.183.582-.495.644-.869l.214-1.281z"
              strokeLinecap="round"
            />
            <path d="M15 12a3 3 0 11-6 0 3 3 0 016 0z" strokeLinecap="round" />
          </svg>
        </Button>
      </Tooltip>

      {/* Share Button — copies the current view URL to the clipboard */}
      <Tooltip content="Copy link to this view" placement="bottom">
        <Button
          className={`${buttonBaseClass} ${glassButton()}`}
          color="default"
          variant="light"
          onPress={handleShare}
        >
          <svg
            className="w-5 h-5"
            fill="none"
            stroke="currentColor"
            strokeWidth={1.5}
            viewBox="0 0 24 24"
          >
            <path
              d="M7.217 10.907a2.25 2.25 0 100 2.186m0-2.186c.18.324.283.696.283 1.093s-.103.77-.283 1.093m0-2.186l9.566-5.314m-9.566 7.5l9.566 5.314m0 0a2.25 2.25 0 103.935 2.186 2.25 2.25 0 00-3.935-2.186zm0-12.814a2.25 2.25 0 103.933-2.185 2.25 2.25 0 00-3.933 2.185z"
              strokeLinecap="round"
              strokeLinejoin="round"
            />
          </svg>
        </Button>
      </Tooltip>

      {/* Advanced Visualization Panel */}
      {isAdvancedOpen && (
        <AdvancedVizPanel
          controlsRef={topRightRef}
          placement="top-right"
          onClose={() => setIsAdvancedOpen(false)}
        />
      )}

      {/* Camera Panel */}
      {isCameraOpen && (
        <CameraPanel
          canExportBox={!!dataset && !dataset.metadata?.wasNormalized}
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
          is3DDataset={dataset?.spatial?.dimensions === 3}
          placement="top-right"
          sceneRotation={sceneRotation}
          setFlipX={setFlipX}
          setFlipY={setFlipY}
          setSceneRotation={setSceneRotation}
          setTargetFilterEnabled={setTargetFilterEnabled}
          setTargetFilterRadius={setTargetFilterRadius}
          setViewMode={setViewMode}
          targetFilterEnabled={targetFilterEnabled}
          targetFilterRadius={targetFilterRadius}
          viewMode={viewMode}
          onClose={() => setIsCameraOpen(false)}
          onResetCamera={resetCamera}
        />
      )}
    </div>
    </>
  );
}

function DotSizeSlider({
  scSizeScale,
  setScSizeScale,
  scSmMode,
}: {
  scSizeScale: number;
  setScSizeScale: (n: number) => void;
  scSmMode: ScSmMode;
}) {
  // SC: editable range via right-click popover. SM: static range from config.
  const sc = useSliderRange(
    "sizeScale",
    VISUALIZATION_CONFIG.SINGLE_CELL_SIZE_SCALE_MIN,
    VISUALIZATION_CONFIG.SINGLE_CELL_SIZE_SCALE_MAX,
    scSizeScale,
  );
  const smGlobalScale = usePanelSingleMoleculeVisualizationStore(
    (s) => s.globalScale,
  );
  const setSmGlobalScale = usePanelSingleMoleculeVisualizationStore(
    (s) => s.setGlobalScale,
  );

  const isSc = scSmMode === "sc";
  const value = isSc ? scSizeScale : smGlobalScale;
  const setValue = isSc ? setScSizeScale : setSmGlobalScale;
  const min = isSc
    ? sc.min
    : VISUALIZATION_CONFIG.SINGLE_MOLECULE_GLOBAL_SCALE_MIN;
  const max = isSc
    ? sc.max
    : VISUALIZATION_CONFIG.SINGLE_MOLECULE_GLOBAL_SCALE_MAX;
  const step = isSc
    ? VISUALIZATION_CONFIG.SINGLE_CELL_SIZE_SCALE_STEP
    : VISUALIZATION_CONFIG.SINGLE_MOLECULE_GLOBAL_SCALE_STEP;
  const borderColor = isSc
    ? "border-blue-400 shadow-[0_0_12px_rgba(59,130,246,0.35)]"
    : "border-purple-400 shadow-[0_0_12px_rgba(168,85,247,0.35)]";

  return (
    <Tooltip
      content={
        isSc
          ? "Change dotsize (right-click to edit range)"
          : "Change SM molecule size"
      }
      placement="right"
    >
      <div
        className={`w-14 h-32 rounded-full border-2 p-2 flex flex-col items-center justify-center transition-colors transition-shadow ${borderColor} ${glassButton()}`}
        onContextMenu={isSc ? sc.onContextMenu : undefined}
      >
        <Slider
          aria-label="Dot size"
          className="h-full"
          maxValue={max}
          minValue={min}
          orientation="vertical"
          size="sm"
          step={step}
          value={value}
          onChange={(v) => setValue(v as number)}
        />
        {isSc && sc.popover}
      </div>
    </Tooltip>
  );
}

function LayerVisibilityButton({
  label,
  kind,
  visible,
  onToggle,
  buttonBaseClass,
}: {
  label: string;
  kind: string;
  visible: boolean;
  onToggle: () => void;
  buttonBaseClass: string;
}) {
  return (
    <Tooltip
      content={`${visible ? "Hide" : "Show"} all ${kind}`}
      placement="right"
    >
      <Button
        className={`${buttonBaseClass} ${visible ? glassButton() : ""}`}
        color={visible ? "default" : "danger"}
        variant={visible ? "light" : "flat"}
        onPress={onToggle}
      >
        <div className="flex flex-col items-center gap-0.5">
          {visible ? (
            <svg
              className="w-4 h-4"
              fill="none"
              stroke="currentColor"
              strokeWidth={1.5}
              viewBox="0 0 24 24"
            >
              <path
                d="M2.036 12.322a1.012 1.012 0 010-.639C3.423 7.51 7.36 4.5 12 4.5c4.638 0 8.573 3.007 9.963 7.178.07.207.07.431 0 .639C20.577 16.49 16.64 19.5 12 19.5c-4.638 0-8.573-3.007-9.963-7.178z"
                strokeLinecap="round"
                strokeLinejoin="round"
              />
              <path
                d="M15 12a3 3 0 11-6 0 3 3 0 016 0z"
                strokeLinecap="round"
                strokeLinejoin="round"
              />
            </svg>
          ) : (
            <svg
              className="w-4 h-4"
              fill="none"
              stroke="currentColor"
              strokeWidth={1.5}
              viewBox="0 0 24 24"
            >
              <path
                d="M3.98 8.223A10.477 10.477 0 001.934 12C3.226 16.338 7.244 19.5 12 19.5c.993 0 1.953-.138 2.863-.395M6.228 6.228A10.45 10.45 0 0112 4.5c4.756 0 8.773 3.162 10.065 7.498a10.523 10.523 0 01-4.293 5.774M6.228 6.228L3 3m3.228 3.228l3.65 3.65m7.894 7.894L21 21m-3.228-3.228l-3.65-3.65m0 0a3 3 0 10-4.243-4.243m4.242 4.242L9.88 9.88"
                strokeLinecap="round"
                strokeLinejoin="round"
              />
            </svg>
          )}
          <span className="text-[10px] font-semibold leading-none">{label}</span>
        </div>
      </Button>
    </Tooltip>
  );
}

