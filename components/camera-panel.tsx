"use client";

import type { StandardizedDataset } from "@/lib/StandardizedDataset";

import { Button } from "@heroui/button";
import { Input } from "@heroui/input";
import { Select, SelectItem } from "@heroui/react";
import { Slider } from "@heroui/react";
import { useEffect, useMemo, useRef, useState } from "react";

import { glassPanel } from "@/components/primitives";
import {
  usePanelDatasetStore,
  usePanelVisualizationStore,
} from "@/lib/hooks/usePanelStores";
import { getEffectiveColumnType } from "@/lib/utils/column-type-utils";
import { loadClusterColumn } from "@/lib/utils/load-cluster-column";
import {
  buildSpatialBinGz,
  computeSampleCentroids,
  IDENTITY_TRANSFORM,
  type SampleTransform,
} from "@/lib/utils/sample-transforms";

const DEG_TO_RAD = Math.PI / 180;
const RAD_TO_DEG = 180 / Math.PI;

interface CameraPanelProps {
  onClose: () => void;
  controlsRef?: React.RefObject<HTMLDivElement>;
  sceneRotation: number;
  setSceneRotation: (degrees: number) => void;
  flipX: boolean;
  setFlipX: (flip: boolean) => void;
  flipY: boolean;
  setFlipY: (flip: boolean) => void;
  viewMode?: string;
  setViewMode?: (mode: any) => void;
  is3DDataset?: boolean;
  // Export box tool — provided by viewers that support exact-mm screenshots.
  // Section only renders when canExportBox is true and viewMode === "2D".
  canExportBox?: boolean;
  exportBox?: {
    enabled: boolean;
    widthMm: number;
    heightMm: number;
    setEnabled: (e: boolean) => void;
    setWidthMm: (mm: number) => void;
    setHeightMm: (mm: number) => void;
    resetCenter: () => void;
  };
  // Optional reset-view action. When provided, a "Reset View" button is shown.
  onResetCamera?: () => void;
}

export function CameraPanel({
  onClose,
  controlsRef,
  sceneRotation,
  setSceneRotation,
  flipX,
  setFlipX,
  flipY,
  setFlipY,
  viewMode,
  setViewMode,
  is3DDataset = false,
  canExportBox = false,
  exportBox,
  onResetCamera,
}: CameraPanelProps) {
  const panelRef = useRef<HTMLDivElement>(null);
  const [alignOpen, setAlignOpen] = useState(false);
  const [exportBoxOpen, setExportBoxOpen] = useState(false);

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

  return (
    <div
      ref={panelRef}
      className={`absolute top-0 left-16 z-[var(--z-panel)] w-[300px] ${glassPanel()}`}
    >
      <div className="p-4 space-y-4">
        {/* Header */}
        <div className="flex justify-between items-center">
          <span className="text-sm font-medium">Camera</span>
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

        {/* Reset View — restores the camera lookAt to the initial framing. */}
        {onResetCamera && (
          <div>
            <Button
              className="w-full text-xs"
              size="sm"
              variant="flat"
              onPress={onResetCamera}
            >
              Reset View (R)
            </Button>
            <p className="text-[10px] text-default-400 mt-1">
              Cmd/Ctrl+click a cell to center on it
            </p>
          </div>
        )}

        {/* Rotation */}
        <div>
          <div className="flex justify-between text-xs mb-1">
            <span className="text-default-500">Rotation</span>
            <span className="text-default-400">{sceneRotation}°</span>
          </div>
          <Slider
            aria-label="Rotation"
            className="w-full"
            maxValue={360}
            minValue={0}
            size="sm"
            step={1}
            value={sceneRotation}
            onChange={(v) => {
              let val = v as number;
              if ((window as any).__shiftHeld) {
                val = Math.round(val / 45) * 45;
              }
              setSceneRotation(val);
            }}
          />
          <p className="text-[10px] text-default-400 mt-1">Hold Shift for 45° snapping</p>
        </div>

        {/* Flip */}
        <div>
          <span className="text-xs text-default-500 mb-2 block">Flip</span>
          <div className="flex gap-2">
            <Button
              className="flex-1 text-xs"
              color={flipX ? "primary" : "default"}
              size="sm"
              variant={flipX ? "solid" : "bordered"}
              onPress={() => setFlipX(!flipX)}
            >
              ⇔ Horizontal
            </Button>
            <Button
              className="flex-1 text-xs"
              color={flipY ? "primary" : "default"}
              size="sm"
              variant={flipY ? "solid" : "bordered"}
              onPress={() => setFlipY(!flipY)}
            >
              ⇕ Vertical
            </Button>
          </div>
        </div>

        {/* 2D/3D Toggle — always available; 2D datasets render as a flat plane in 3D. */}
        {setViewMode && viewMode && (
          <div>
            <span className="text-xs text-default-500 mb-2 block">View Mode</span>
            <div className="flex gap-2">
              <Button
                className="flex-1 text-xs"
                color={viewMode === "2D" ? "primary" : "default"}
                size="sm"
                variant={viewMode === "2D" ? "solid" : "bordered"}
                onPress={() => setViewMode("2D")}
              >
                2D
              </Button>
              <Button
                className="flex-1 text-xs"
                color={viewMode === "3D" ? "primary" : "default"}
                size="sm"
                variant={viewMode === "3D" ? "solid" : "bordered"}
                onPress={() => setViewMode("3D")}
              >
                3D
              </Button>
            </div>
          </div>
        )}

        {/* Reset */}
        <Button
          className="w-full"
          color="danger"
          size="sm"
          variant="flat"
          onPress={() => {
            setSceneRotation(0);
            setFlipX(false);
            setFlipY(false);
          }}
        >
          Reset
        </Button>

        {/* Export box (exact-mm screenshot tool) */}
        {canExportBox && exportBox && viewMode === "2D" && (
          <>
            <button
              aria-expanded={exportBoxOpen}
              className="w-full flex items-center justify-between text-xs text-default-400 hover:text-default-200 pt-1"
              onClick={() => setExportBoxOpen(!exportBoxOpen)}
            >
              <span>Export box</span>
              <span
                className={`transition-transform ${exportBoxOpen ? "rotate-90" : ""}`}
              >
                ▸
              </span>
            </button>
            {exportBoxOpen && (
              <ExportBoxSection
                heightMm={exportBox.heightMm}
                enabled={exportBox.enabled}
                resetCenter={exportBox.resetCenter}
                setEnabled={exportBox.setEnabled}
                setHeightMm={exportBox.setHeightMm}
                setWidthMm={exportBox.setWidthMm}
                widthMm={exportBox.widthMm}
              />
            )}
          </>
        )}

        {/* Per-sample align (advanced, collapsed by default) */}
        <button
          aria-expanded={alignOpen}
          className="w-full flex items-center justify-between text-xs text-default-400 hover:text-default-200 pt-1"
          onClick={() => setAlignOpen(!alignOpen)}
        >
          <span>Per-sample align</span>
          <span className={`transition-transform ${alignOpen ? "rotate-90" : ""}`}>
            ▸
          </span>
        </button>
        {alignOpen && <SampleAlignSection />}
      </div>
    </div>
  );
}

function SampleAlignSection() {
  const { getCurrentDataset } = usePanelDatasetStore();
  const {
    transformColumn,
    setTransformColumn,
    activeSampleId,
    setActiveSampleId,
    sampleTransforms,
    setSampleTransform,
    clearSampleTransform,
    clearAllSampleTransforms,
    columnTypeOverrides,
    clusterVersion,
    incrementClusterVersion,
  } = usePanelVisualizationStore();

  const rawDataset = getCurrentDataset();
  const dataset =
    rawDataset && "clusters" in rawDataset
      ? (rawDataset as StandardizedDataset)
      : null;

  const [loading, setLoading] = useState(false);
  const [exporting, setExporting] = useState(false);

  const candidates = useMemo(() => {
    if (!dataset) return [] as string[];
    const all =
      dataset.allClusterColumnNames?.length > 0
        ? dataset.allClusterColumnNames
        : (dataset.clusters ?? []).map((c) => c.column);
    return all.filter(
      (name) =>
        getEffectiveColumnType(name, dataset, columnTypeOverrides) ===
        "categorical",
    );
  }, [dataset, columnTypeOverrides, clusterVersion]);

  useEffect(() => {
    if (transformColumn || !candidates.length) return;
    const pick = candidates.includes("_sample_id")
      ? "_sample_id"
      : candidates[0];
    setTransformColumn(pick);
  }, [candidates, transformColumn, setTransformColumn]);

  useEffect(() => {
    if (!dataset || !transformColumn) return;
    if (dataset.clusters?.some((c) => c.column === transformColumn)) return;
    let cancelled = false;
    setLoading(true);
    loadClusterColumn(dataset, transformColumn)
      .then((added) => {
        if (!cancelled && added) incrementClusterVersion();
      })
      .catch((e) => {
        // eslint-disable-next-line no-console
        console.warn(`Failed to load column ${transformColumn}`, e);
      })
      .finally(() => {
        if (!cancelled) setLoading(false);
      });
    return () => {
      cancelled = true;
    };
  }, [dataset, transformColumn, incrementClusterVersion]);

  const cluster = useMemo(() => {
    if (!dataset || !transformColumn) return null;
    return dataset.clusters?.find((c) => c.column === transformColumn) ?? null;
  }, [dataset, transformColumn, clusterVersion]);

  const samples = useMemo<string[]>(() => {
    if (!cluster) return [];
    const raw =
      cluster.uniqueValues && cluster.uniqueValues.length > 0
        ? cluster.uniqueValues.slice()
        : Array.from(new Set(cluster.values.map(String)));
    return raw.sort((a, b) =>
      a.localeCompare(b, undefined, { numeric: true }),
    );
  }, [cluster]);

  useEffect(() => {
    if (activeSampleId || !samples.length) return;
    setActiveSampleId(samples[0]);
  }, [samples, activeSampleId, setActiveSampleId]);

  const currentTransform: SampleTransform =
    (activeSampleId && sampleTransforms.get(activeSampleId)) ||
    IDENTITY_TRANSFORM;

  const onChangeTransform = (patch: Partial<SampleTransform>) => {
    if (!activeSampleId) return;
    setSampleTransform(activeSampleId, { ...currentTransform, ...patch });
  };

  const onExport = async () => {
    if (!dataset) return;
    setExporting(true);
    try {
      let centroids = null;
      if (cluster) {
        centroids = computeSampleCentroids(dataset.spatial, cluster);
      }
      const blob = await buildSpatialBinGz(
        dataset.spatial,
        cluster,
        centroids,
        sampleTransforms,
      );
      const url = URL.createObjectURL(blob);
      const a = document.createElement("a");
      a.href = url;
      a.download = "spatial_aligned.bin.gz";
      document.body.appendChild(a);
      a.click();
      a.remove();
      URL.revokeObjectURL(url);
    } catch (e) {
      // eslint-disable-next-line no-console
      console.error("Export failed", e);
    } finally {
      setExporting(false);
    }
  };

  const editedCount = sampleTransforms.size;

  if (!dataset) {
    return (
      <div className="text-xs text-default-400 py-2">Load a dataset first.</div>
    );
  }
  if (!candidates.length) {
    return (
      <div className="text-xs text-default-400 py-2">
        No categorical obs column available.
      </div>
    );
  }

  return (
    <div className="space-y-2 pt-1">
      <p className="text-[10px] text-default-400">
        Translate + rotate per sample (pivot = centroid). Values in raw coord
        units.
      </p>

      <Select
        label="Sample column"
        labelPlacement="outside"
        selectedKeys={transformColumn ? [transformColumn] : []}
        size="sm"
        onChange={(e) => setTransformColumn(e.target.value || null)}
      >
        {candidates.map((name) => (
          <SelectItem key={name}>{name}</SelectItem>
        ))}
      </Select>

      <Select
        isDisabled={!samples.length || loading}
        label={
          loading
            ? `Loading "${transformColumn}"…`
            : `Sample (${samples.length})`
        }
        labelPlacement="outside"
        selectedKeys={activeSampleId ? [activeSampleId] : []}
        size="sm"
        onChange={(e) => setActiveSampleId(e.target.value || null)}
      >
        {samples.map((s) => (
          <SelectItem
            key={s}
            endContent={
              sampleTransforms.has(s) ? (
                <span className="text-xs text-primary-400">●</span>
              ) : null
            }
          >
            {s}
          </SelectItem>
        ))}
      </Select>

      <div className="grid grid-cols-3 gap-1">
        <Input
          isDisabled={!activeSampleId}
          label="dx"
          labelPlacement="outside"
          size="sm"
          step="any"
          type="number"
          value={String(currentTransform.dx)}
          onValueChange={(v) =>
            onChangeTransform({ dx: parseFloatOr(v, 0) })
          }
        />
        <Input
          isDisabled={!activeSampleId}
          label="dy"
          labelPlacement="outside"
          size="sm"
          step="any"
          type="number"
          value={String(currentTransform.dy)}
          onValueChange={(v) =>
            onChangeTransform({ dy: parseFloatOr(v, 0) })
          }
        />
        <Input
          endContent={<span className="text-[10px] text-default-400">°</span>}
          isDisabled={!activeSampleId}
          label="rot"
          labelPlacement="outside"
          size="sm"
          step="any"
          type="number"
          value={String(currentTransform.theta * RAD_TO_DEG)}
          onValueChange={(v) =>
            onChangeTransform({ theta: parseFloatOr(v, 0) * DEG_TO_RAD })
          }
        />
      </div>

      <div className="flex items-center gap-1">
        <Button
          className="flex-1 text-xs"
          isDisabled={
            !activeSampleId || !sampleTransforms.has(activeSampleId)
          }
          size="sm"
          variant="flat"
          onPress={() => activeSampleId && clearSampleTransform(activeSampleId)}
        >
          Reset
        </Button>
        <Button
          className="flex-1 text-xs"
          color="danger"
          isDisabled={editedCount === 0}
          size="sm"
          variant="flat"
          onPress={clearAllSampleTransforms}
        >
          Reset all ({editedCount})
        </Button>
      </div>

      <Button
        className="w-full"
        color="primary"
        isDisabled={exporting}
        isLoading={exporting}
        size="sm"
        title="Bake transforms and download spatial_aligned.bin.gz"
        variant="solid"
        onPress={onExport}
      >
        Export spatial_aligned.bin.gz
      </Button>
    </div>
  );
}

function parseFloatOr(s: string, fallback: number): number {
  if (s === "" || s === "-" || s === "." || s === "-.") return fallback;
  const n = Number(s);
  return Number.isFinite(n) ? n : fallback;
}

interface ExportBoxSectionProps {
  enabled: boolean;
  widthMm: number;
  heightMm: number;
  setEnabled: (e: boolean) => void;
  setWidthMm: (mm: number) => void;
  setHeightMm: (mm: number) => void;
  resetCenter: () => void;
}

function ExportBoxSection({
  enabled,
  widthMm,
  heightMm,
  setEnabled,
  setWidthMm,
  setHeightMm,
  resetCenter,
}: ExportBoxSectionProps) {
  // Local string state so the user can type freely (e.g. "0." while heading
  // toward "0.5") without the store rounding/clamping every keystroke.
  const [wStr, setWStr] = useState(String(widthMm));
  const [hStr, setHStr] = useState(String(heightMm));
  useEffect(() => {
    setWStr(String(widthMm));
  }, [widthMm]);
  useEffect(() => {
    setHStr(String(heightMm));
  }, [heightMm]);

  const commitW = () => {
    const n = parseFloat(wStr);
    if (Number.isFinite(n) && n > 0) setWidthMm(n);
    else setWStr(String(widthMm));
  };
  const commitH = () => {
    const n = parseFloat(hStr);
    if (Number.isFinite(n) && n > 0) setHeightMm(n);
    else setHStr(String(heightMm));
  };

  return (
    <div className="space-y-2 text-xs">
      <label className="flex items-center gap-2 cursor-pointer">
        <input
          checked={enabled}
          className="accent-yellow-300"
          type="checkbox"
          onChange={(e) => setEnabled(e.target.checked)}
        />
        <span>Show export box</span>
      </label>

      <div className="grid grid-cols-2 gap-2">
        <div>
          <div className="text-default-400 mb-1">Width (mm)</div>
          <Input
            classNames={{ input: "text-xs" }}
            min={0.01}
            size="sm"
            step={0.1}
            type="number"
            value={wStr}
            onBlur={commitW}
            onChange={(e) => setWStr(e.target.value)}
            onKeyDown={(e) => {
              if (e.key === "Enter") commitW();
            }}
          />
        </div>
        <div>
          <div className="text-default-400 mb-1">Height (mm)</div>
          <Input
            classNames={{ input: "text-xs" }}
            min={0.01}
            size="sm"
            step={0.1}
            type="number"
            value={hStr}
            onBlur={commitH}
            onChange={(e) => setHStr(e.target.value)}
            onKeyDown={(e) => {
              if (e.key === "Enter") commitH();
            }}
          />
        </div>
      </div>

      <Button
        className="w-full"
        isDisabled={!enabled}
        size="sm"
        variant="flat"
        onPress={resetCenter}
      >
        Recenter
      </Button>

      <p className="text-[10px] text-default-400 leading-snug">
        Drag the yellow box to position. Use Save / Copy on the box toolbar to
        export the cropped region at the current zoom.
      </p>
    </div>
  );
}
