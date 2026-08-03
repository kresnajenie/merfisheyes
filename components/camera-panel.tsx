"use client";

import type { StandardizedDataset } from "@/lib/StandardizedDataset";

import { Button } from "@heroui/button";
import { Input } from "@heroui/input";
import { Select, SelectItem } from "@heroui/react";
import { Slider, Switch } from "@heroui/react";
import { useEffect, useMemo, useRef, useState } from "react";
import { useSession } from "next-auth/react";
import { toast } from "react-toastify";

import { glassPanel } from "@/components/primitives";
import { ThumbnailCaptureModal } from "@/components/thumbnail-capture-modal";
import { GeneMultiSelect } from "@/components/gene-multiselect";
import {
  usePanelDatasetStore,
  usePanelVisualizationStore,
  usePanelSingleMoleculeStore,
  usePanelSingleMoleculeVisualizationStore,
  usePanelId,
} from "@/lib/hooks/usePanelStores";
import { useViewerRegistrationStore } from "@/lib/stores/viewerRegistrationStore";
import { useSingleMoleculeStore } from "@/lib/stores/singleMoleculeStore";
import { useSingleMoleculeVisualizationStore } from "@/lib/stores/singleMoleculeVisualizationStore";
import {
  extractViewerConfig,
  type ViewerConfig,
} from "@/lib/utils/viewer-config";
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
// Per-sample translations are stored in raw coordinate units (microns); the UI
// shows them in millimetres. display_mm = microns * MICRONS_TO_MM.
const MICRONS_TO_MM = 1 / 1000;

interface CameraPanelProps {
  onClose: () => void;
  controlsRef?: React.RefObject<HTMLDivElement>;
  // Where the flyout anchors: "rail" opens to the right of the left rail;
  // "top-right" drops below the top-right control cluster.
  placement?: "rail" | "top-right";
  // Which viewer this panel serves. "sc" (default) uses the single-cell stores
  // for per-sample align + owner defaults; "sm" uses the single-molecule store
  // for a camera-only owner preset and hides the SC-only align section.
  viewerKind?: "sc" | "sm";
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
  // Optional distance-from-target filter controls. When provided AND
  // viewMode === "3D", a "Distance filter" section is shown.
  targetFilterEnabled?: boolean;
  setTargetFilterEnabled?: (b: boolean) => void;
  targetFilterRadius?: number;
  setTargetFilterRadius?: (r: number) => void;
}

export function CameraPanel({
  onClose,
  controlsRef,
  placement = "rail",
  viewerKind = "sc",
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
  targetFilterEnabled,
  setTargetFilterEnabled,
  targetFilterRadius,
  setTargetFilterRadius,
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

  const posClass =
    placement === "top-right"
      ? "top-14 right-0 max-h-[calc(100vh-6rem)] overflow-y-auto"
      : "top-0 left-16";

  return (
    <div
      ref={panelRef}
      className={`absolute ${posClass} z-[var(--z-panel)] w-[300px] ${glassPanel()}`}
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

        {/* Distance filter — 3D only. Fades points beyond the radius from the
            current orbit target (cmd+click recenter point). */}
        {viewMode === "3D" &&
          setTargetFilterEnabled &&
          setTargetFilterRadius &&
          targetFilterRadius !== undefined && (
            <div>
              <div className="flex items-center justify-between mb-2">
                <span className="text-xs text-default-500">
                  Distance filter
                </span>
                <Switch
                  isSelected={!!targetFilterEnabled}
                  size="sm"
                  onValueChange={setTargetFilterEnabled}
                />
              </div>
              {targetFilterEnabled && (
                <>
                  <div className="flex justify-between items-center text-xs mb-1">
                    <span className="text-default-500">Radius</span>
                    <div className="flex items-center gap-1">
                      <input
                        aria-label="Distance filter radius value"
                        className="w-16 px-1 py-0.5 text-xs rounded bg-default-100/60 border border-default-300/40 text-right text-default-700 outline-none focus:border-primary"
                        step={10}
                        type="number"
                        value={Math.round(targetFilterRadius)}
                        onChange={(e) => {
                          const n = parseFloat(e.target.value);

                          if (!Number.isNaN(n)) {
                            setTargetFilterRadius(Math.max(0, n));
                          }
                        }}
                      />
                      <span className="text-default-400">µm</span>
                    </div>
                  </div>
                  <Slider
                    aria-label="Distance filter radius"
                    className="w-full"
                    maxValue={30000}
                    minValue={10}
                    size="sm"
                    step={10}
                    value={Math.min(30000, Math.max(10, targetFilterRadius))}
                    onChange={(v) => setTargetFilterRadius(v as number)}
                  />
                  <p className="text-[10px] text-default-400 mt-1">
                    Radius from cmd+click target. 10% feather at the edge.
                  </p>
                </>
              )}
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

        {/* Per-sample align (advanced, collapsed by default) — single-cell only. */}
        {viewerKind === "sc" && (
          <>
            <button
              aria-expanded={alignOpen}
              className="w-full flex items-center justify-between text-xs text-default-400 hover:text-default-200 pt-1"
              onClick={() => setAlignOpen(!alignOpen)}
            >
              <span>Per-sample align</span>
              <span
                className={`transition-transform ${alignOpen ? "rotate-90" : ""}`}
              >
                ▸
              </span>
            </button>
            {alignOpen && <SampleAlignSection />}
          </>
        )}

        {viewerKind === "sm" ? (
          <SMOwnerDefaultsSection />
        ) : (
          <OwnerDefaultsSection />
        )}
      </div>
    </div>
  );
}

/**
 * Merge a partial into the dataset's saved viewerConfig and PATCH it, so a
 * single-field owner control (priority column, default genes) doesn't clobber
 * other saved fields. Updates the registration store on success.
 */
async function patchViewerConfig(
  dbId: string,
  patch: Partial<ViewerConfig>,
): Promise<boolean> {
  const current = (useViewerRegistrationStore.getState().viewerConfig ?? {
    version: 1,
  }) as ViewerConfig;
  const merged: ViewerConfig = { ...current, ...patch, version: 1 };

  try {
    const res = await fetch(`/api/ingest/mine/${dbId}`, {
      method: "PATCH",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ viewerConfig: merged }),
    });

    if (!res.ok) return false;
    useViewerRegistrationStore.getState().set({ viewerConfig: merged });

    return true;
  } catch {
    return false;
  }
}

/**
 * Owner-only "Save current view as default" control. Persists the current
 * rotation, flips, priority cluster column, per-sample transforms, and custom
 * colors to the dataset's viewerConfig (applied on future loads). Only shown to
 * the dataset owner (or an admin) on the primary panel.
 */
function OwnerDefaultsSection() {
  const { data: session } = useSession();
  const panelId = usePanelId();
  const { dbId, ownerId, adminOwned } = useViewerRegistrationStore();
  const { getCurrentDataset } = usePanelDatasetStore();
  const store = usePanelVisualizationStore();
  const [saving, setSaving] = useState(false);
  const [thumbOpen, setThumbOpen] = useState(false);
  const [overlayGenes, setOverlayGenes] = useState<string[]>([]);
  const [savingOverlay, setSavingOverlay] = useState(false);

  const dataset = getCurrentDataset() as StandardizedDataset | null;

  // SM overlay (combined SC+SM viewer) — from the GLOBAL SM stores, matching
  // how three-scene renders the overlay. Present only when a mapping.json
  // "__all__" overlay has loaded.
  const smOverlay = useSingleMoleculeStore((s) =>
    s.currentDatasetId ? s.datasets.get(s.currentDatasetId) ?? null : null,
  );
  const smSelectedGenes = useSingleMoleculeVisualizationStore(
    (s) => s.selectedGenes,
  );
  const overlayUniqueGenes = smOverlay?.uniqueGenes ?? [];

  // Seed the overlay picker from this SC dataset's saved overlay genes, else
  // the overlay's current selection, when the overlay loads.
  useEffect(() => {
    if (!smOverlay) return;
    const saved = (
      useViewerRegistrationStore.getState().viewerConfig as ViewerConfig | null
    )?.defaultGenes;

    setOverlayGenes(
      saved && saved.length > 0
        ? saved.filter((g) => overlayUniqueGenes.includes(g))
        : Array.from(smSelectedGenes.keys()),
    );
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [smOverlay?.id]);

  // Categorical columns selectable as the default (priority) cell-type column.
  const priorityCandidates = useMemo(() => {
    if (!dataset) return [] as string[];
    const all =
      dataset.allClusterColumnNames?.length > 0
        ? dataset.allClusterColumnNames
        : (dataset.clusters ?? []).map((c) => c.column);

    return all.filter(
      (name) =>
        getEffectiveColumnType(name, dataset, store.columnTypeOverrides) ===
        "categorical",
    );
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [dataset, store.columnTypeOverrides, store.clusterVersion]);

  const userId = session?.user?.id ?? null;
  const isAdmin =
    session?.user?.role === "ADMIN" || session?.user?.role === "SUPER_ADMIN";
  // Personal owner edits their own; any admin edits admin-owned (shared).
  const canSave =
    !panelId &&
    !!dbId &&
    !!dataset &&
    ((!!ownerId && ownerId === userId) || (adminOwned && isAdmin));

  if (!canSave) return null;

  const save = async () => {
    if (!dataset || !dbId) return;
    setSaving(true);
    try {
      const viewerConfig = extractViewerConfig(store as any, dataset);
      const res = await fetch(`/api/ingest/mine/${dbId}`, {
        method: "PATCH",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ viewerConfig }),
      });

      if (!res.ok) {
        toast.error("Couldn't save defaults.");

        return;
      }
      useViewerRegistrationStore.getState().set({ viewerConfig });
      toast.success("Saved as this dataset's defaults.");
    } catch {
      toast.error("Couldn't save defaults.");
    } finally {
      setSaving(false);
    }
  };

  return (
    <div className="border-t border-white/10 pt-3 space-y-2">
      {priorityCandidates.length > 0 && (
        <Select
          label="Default cell-type column"
          labelPlacement="outside"
          selectedKeys={store.selectedColumn ? [store.selectedColumn] : []}
          size="sm"
          onChange={(e) => {
            const col = e.target.value;

            if (!col || !dbId) return;
            store.setSelectedColumn(col, false);
            patchViewerConfig(dbId, { priorityColumn: col }).then((ok) =>
              ok
                ? toast.success("Default column saved.")
                : toast.error("Couldn't save default column."),
            );
          }}
        >
          {priorityCandidates.map((name) => (
            <SelectItem key={name}>{name}</SelectItem>
          ))}
        </Select>
      )}

      <Button
        className="w-full text-xs"
        color="primary"
        isLoading={saving}
        size="sm"
        variant="flat"
        onPress={save}
      >
        Save current view as default
      </Button>
      <p className="text-[10px] text-default-400 mt-1">
        Saves rotation, colors, alignment, and the{" "}
        <span className="text-default-300">
          {store.selectedColumn ?? "current"}
        </span>{" "}
        column as this dataset&apos;s defaults on load.
      </p>

      {smOverlay && overlayUniqueGenes.length > 0 && (
        <div className="pt-1 space-y-1">
          <GeneMultiSelect
            genes={overlayUniqueGenes}
            label="Default overlay genes"
            value={overlayGenes}
            onChange={setOverlayGenes}
          />
          <Button
            className="w-full text-xs"
            color="primary"
            isDisabled={overlayGenes.length === 0}
            isLoading={savingOverlay}
            size="sm"
            variant="flat"
            onPress={async () => {
              if (!dbId) return;
              setSavingOverlay(true);
              const ok = await patchViewerConfig(dbId, {
                defaultGenes: overlayGenes,
              });

              setSavingOverlay(false);
              if (ok) toast.success("Overlay genes saved.");
              else toast.error("Couldn't save overlay genes.");
            }}
          >
            Save overlay genes
          </Button>
        </div>
      )}

      <Button
        className="w-full text-xs"
        size="sm"
        variant="flat"
        onPress={() => setThumbOpen(true)}
      >
        Set thumbnail from screenshot
      </Button>

      {dbId && (
        <ThumbnailCaptureModal
          dbId={dbId}
          isOpen={thumbOpen}
          viewerKind="sc"
          onClose={() => setThumbOpen(false)}
        />
      )}
    </div>
  );
}

/**
 * Single-molecule owner preset: saves the camera (rotation, flips, 2D/3D) to
 * the dataset's viewerConfig, applied on load. Owner/admin only, primary panel.
 */
function SMOwnerDefaultsSection() {
  const { data: session } = useSession();
  const panelId = usePanelId();
  const { dbId, ownerId, adminOwned } = useViewerRegistrationStore();
  const store = usePanelSingleMoleculeVisualizationStore();
  const { getCurrentDataset } = usePanelSingleMoleculeStore();
  const [saving, setSaving] = useState(false);
  const [savingGenes, setSavingGenes] = useState(false);
  const [thumbOpen, setThumbOpen] = useState(false);
  const [defaultGenes, setDefaultGenes] = useState<string[]>([]);

  const smDataset = getCurrentDataset();
  const uniqueGenes = smDataset?.uniqueGenes ?? [];

  // Seed the picker from the saved defaults, else the current selection.
  useEffect(() => {
    const saved = (
      useViewerRegistrationStore.getState().viewerConfig as ViewerConfig | null
    )?.defaultGenes;

    setDefaultGenes(
      saved && saved.length > 0 ? saved : Array.from(store.selectedGenes.keys()),
    );
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  const userId = session?.user?.id ?? null;
  const isAdmin =
    session?.user?.role === "ADMIN" || session?.user?.role === "SUPER_ADMIN";
  const canSave =
    !panelId &&
    !!dbId &&
    ((!!ownerId && ownerId === userId) || (adminOwned && isAdmin));

  if (!canSave) return null;

  const save = async () => {
    if (!dbId) return;
    setSaving(true);
    // Merge so a camera save doesn't wipe saved default genes.
    const ok = await patchViewerConfig(dbId, {
      sceneRotation: store.sceneRotation,
      flipX: store.flipX,
      flipY: store.flipY,
      viewMode: store.viewMode,
    });

    setSaving(false);
    if (ok) toast.success("Saved as this dataset's defaults.");
    else toast.error("Couldn't save defaults.");
  };

  const saveGenes = async () => {
    if (!dbId) return;
    setSavingGenes(true);
    const ok = await patchViewerConfig(dbId, {
      defaultGenes: defaultGenes.slice(0, 20),
    });

    setSavingGenes(false);
    if (ok) toast.success("Default genes saved.");
    else toast.error("Couldn't save default genes.");
  };

  return (
    <div className="border-t border-white/10 pt-3 space-y-2">
      <Button
        className="w-full text-xs"
        color="primary"
        isLoading={saving}
        size="sm"
        variant="flat"
        onPress={save}
      >
        Save current view as default
      </Button>
      <p className="text-[10px] text-default-400 mt-1">
        Saves rotation, flips, and 2D/3D view as this dataset&apos;s defaults on
        load.
      </p>

      {uniqueGenes.length > 0 && (
        <div className="pt-1 space-y-1">
          <GeneMultiSelect
            genes={uniqueGenes}
            label="Default genes shown on load"
            value={defaultGenes}
            onChange={setDefaultGenes}
          />
          <Button
            className="w-full text-xs"
            color="primary"
            isDisabled={defaultGenes.length === 0}
            isLoading={savingGenes}
            size="sm"
            variant="flat"
            onPress={saveGenes}
          >
            Save default genes
          </Button>
        </div>
      )}

      <Button
        className="w-full text-xs"
        size="sm"
        variant="flat"
        onPress={() => setThumbOpen(true)}
      >
        Set thumbnail from screenshot
      </Button>

      {dbId && (
        <ThumbnailCaptureModal
          dbId={dbId}
          isOpen={thumbOpen}
          viewerKind="sm"
          onClose={() => setThumbOpen(false)}
        />
      )}
    </div>
  );
}

/** Trim float noise for display (e.g. 1.4999999 -> 1.5). */
function fmtNum(n: number): string {
  if (!Number.isFinite(n)) return "";

  return String(Math.round(n * 1e6) / 1e6);
}

/**
 * Numeric field that keeps a local text buffer so partial input ("-", "1.",
 * empty) doesn't get parsed-and-reset mid-typing. Commits only finite numbers.
 * `displayFactor` converts the stored (canonical) value to the shown value:
 * display = value * displayFactor; value = shown / displayFactor.
 */
function AlignNumberInput({
  label,
  endContent,
  value,
  displayFactor,
  disabled,
  onCommit,
}: {
  label: string;
  endContent?: React.ReactNode;
  value: number;
  displayFactor: number;
  disabled?: boolean;
  onCommit: (value: number) => void;
}) {
  const toDisplay = (v: number) => fmtNum(v * displayFactor);
  const [str, setStr] = useState(() => toDisplay(value));
  const [focused, setFocused] = useState(false);

  // Sync from the store when the field isn't being edited (e.g. sample switch).
  useEffect(() => {
    if (!focused) setStr(toDisplay(value));
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [value, displayFactor, focused]);

  return (
    <Input
      endContent={endContent}
      inputMode="decimal"
      isDisabled={disabled}
      label={label}
      labelPlacement="outside"
      size="sm"
      type="text"
      value={str}
      onBlur={() => {
        setFocused(false);
        const n = parseFloat(str);

        setStr(Number.isFinite(n) ? toDisplay(n / displayFactor) : toDisplay(value));
      }}
      onFocus={() => setFocused(true)}
      onValueChange={(v) => {
        setStr(v);
        const n = parseFloat(v);

        if (Number.isFinite(n)) onCommit(n / displayFactor);
      }}
    />
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
        Translate (mm) + rotate per sample (pivot = centroid).
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
        <AlignNumberInput
          disabled={!activeSampleId}
          displayFactor={MICRONS_TO_MM}
          endContent={<span className="text-[10px] text-default-400">mm</span>}
          label="dx"
          value={currentTransform.dx}
          onCommit={(dx) => onChangeTransform({ dx })}
        />
        <AlignNumberInput
          disabled={!activeSampleId}
          displayFactor={MICRONS_TO_MM}
          endContent={<span className="text-[10px] text-default-400">mm</span>}
          label="dy"
          value={currentTransform.dy}
          onCommit={(dy) => onChangeTransform({ dy })}
        />
        <AlignNumberInput
          disabled={!activeSampleId}
          displayFactor={RAD_TO_DEG}
          endContent={<span className="text-[10px] text-default-400">°</span>}
          label="rot"
          value={currentTransform.theta}
          onCommit={(theta) => onChangeTransform({ theta })}
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
