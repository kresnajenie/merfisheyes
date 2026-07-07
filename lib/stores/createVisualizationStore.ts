import { createStore } from "zustand";

import { VISUALIZATION_CONFIG } from "../config/visualization.config";
import type { SampleTransform } from "../utils/sample-transforms";

export type VisualizationMode = "celltype" | "gene";
export type CellViewMode = "2D" | "3D";

export interface VisualizationState {
  mode: VisualizationMode[];
  panelMode: VisualizationMode;
  viewMode: CellViewMode;
  selectedGene: string | null;
  geneScaleMin: number;
  geneScaleMax: number;
  selectedClusterColumn: string | null;
  selectedColumn: string | null;
  selectedCelltypes: Set<string>;
  // When true, a selected gene is rendered on every cell even while
  // celltypes are selected (bypasses the combined gene+celltype view).
  geneEverywhere: boolean;
  // Celltypes hidden from the 3D scene via the per-badge eye toggle. They
  // stay in selectedCelltypes (so DEG / plots use the full selection) but
  // are greyed out in the point cloud.
  hiddenCelltypes: Set<string>;
  numericalScaleMin: number;
  numericalScaleMax: number;
  celltypeSearchTerm: string;
  geneSearchTerm: string;
  selectedEmbedding: string | null;
  colorPalette: Record<string, string>;
  alphaScale: number;
  sizeScale: number;
  clusterVersion: number;
  deStatsVersion: number;
  // Bumped when the user requests a camera reset (button or shortcut).
  // ThreeScene watches this and snaps controls.target back to the initial lookAt.
  cameraResetSignal: number;
  // 3D orbit (O key). orbitSpeed is a multiplier on the base azimuth period
  // (1.0 = 4.5s/rev); orbitYBob is the vertical bob amplitude as a fraction
  // of the camera→target XZ radius (0 = flat ring, 0.3 = current default).
  orbitSpeed: number;
  orbitYBob: number;
  // Distance-from-target filter. When enabled, points fade out smoothly
  // beyond `targetFilterRadius * (dataset extent)` from the camera's orbit
  // target (cmd+click recenter point). The feather band is 10% of that.
  targetFilterEnabled: boolean;
  targetFilterRadius: number;
  // Master show/hide for the whole single-cell point cloud layer.
  scLayerVisible: boolean;
  degTarget: string | null;
  degReference: string | null; // null = vs Rest
  degTargetAuto: boolean; // target follows the most-recently-selected celltype
  degReferenceAuto: boolean; // reference follows the 2nd-most-recently-selected celltype
  degSearchTerm: string;
  degSortKey: "log2FC" | "meanIn" | "pctIn";
  degSortDesc: boolean;
  degPanelOpen: boolean;
  columnTypeOverrides: Record<string, "categorical" | "numerical">;
  celltypePlayback: boolean;
  celltypePlaybackInterval: number;
  celltypePlaybackSequence: string[];

  // Camera/scene transforms
  sceneRotation: number; // degrees
  flipX: boolean;
  flipY: boolean;

  // Advanced visualization settings
  selectedSizeMultiplier: number;
  greyedOutSizeMultiplier: number;
  greyedOutAlpha: number;
  expressionAlphaMin: number;
  expressionAlphaMax: number;
  pointSizeMultiplierMin: number;
  pointSizeMultiplierMax: number;
  targetPx: number;

  setMode: (mode: VisualizationMode[]) => void;
  setPanelMode: (mode: VisualizationMode) => void;
  setViewMode: (mode: CellViewMode) => void;
  setSelectedGene: (gene: string | null) => void;
  setGeneEverywhere: (everywhere: boolean) => void;
  toggleCelltypeVisibility: (celltype: string) => void;
  soloCelltype: (celltype: string) => void;
  setGeneScaleMin: (min: number) => void;
  setGeneScaleMax: (max: number) => void;
  setNumericalScaleMin: (min: number) => void;
  setNumericalScaleMax: (max: number) => void;
  toggleCelltype: (celltype: string) => void;
  setCelltypes: (celltypes: Set<string>) => void;
  setClusterColumn: (column: string | null) => void;
  setSelectedColumn: (column: string | null, isNumerical?: boolean) => void;
  setSelectedEmbedding: (embedding: string | null) => void;
  setColorPalette: (palette: Record<string, string>) => void;
  setAlphaScale: (alpha: number) => void;
  setSizeScale: (size: number) => void;
  setCelltypeSearchTerm: (value: string) => void;
  setGeneSearchTerm: (value: string) => void;
  incrementClusterVersion: () => void;
  incrementDeStatsVersion: () => void;
  resetCamera: () => void;
  setTargetFilterEnabled: (enabled: boolean) => void;
  setScLayerVisible: (visible: boolean) => void;
  setTargetFilterRadius: (r: number) => void;
  setDegTarget: (target: string | null) => void;
  setDegReference: (reference: string | null) => void;
  setDegTargetAuto: (auto: boolean) => void;
  setDegReferenceAuto: (auto: boolean) => void;
  setDegSearchTerm: (term: string) => void;
  setDegSortKey: (key: "log2FC" | "meanIn" | "pctIn") => void;
  setDegPanelOpen: (open: boolean) => void;
  setDegSortDesc: (desc: boolean) => void;
  toggleColumnType: (
    column: string,
    currentType: "categorical" | "numerical",
  ) => void;
  setColumnTypeOverrides: (
    overrides: Record<string, "categorical" | "numerical">,
  ) => void;
  setCelltypePlayback: (playing: boolean) => void;
  setCelltypePlaybackInterval: (interval: number) => void;
  setCelltypePlaybackSequence: (sequence: string[]) => void;
  setSceneRotation: (degrees: number) => void;
  setFlipX: (flip: boolean) => void;
  setFlipY: (flip: boolean) => void;
  setAdvancedViz: (key: string, value: number) => void;
  pinnedTooltipColumns: Set<string>;
  togglePinnedTooltipColumn: (column: string) => void;
  sliderRanges: Record<string, { min: number; max: number }>;
  setSliderRange: (key: string, min: number, max: number) => void;
  clearSliderRanges: (keys: string[]) => void;
  colormap: string;
  setColormap: (name: string) => void;
  plotPanelOpen: boolean;
  setPlotPanelOpen: (open: boolean) => void;
  // Plot-panel-only secondary grouping (e.g. compare expression across
  // treatments within each starred celltype). Does not affect the 3D scene.
  secondaryColumn: string | null;
  selectedSecondaryValues: Set<string>;
  secondaryPaletteOverrides: Record<string, string>;
  // User-customised ordering for the current secondary column's values. Empty
  // means "use the column's natural order". Resets when secondaryColumn changes.
  secondaryValueOrder: string[];
  setSecondaryColumn: (column: string | null) => void;
  setSelectedSecondaryValues: (values: Set<string>) => void;
  toggleSecondaryValue: (value: string) => void;
  setSecondaryPaletteOverride: (value: string, color: string) => void;
  setSecondaryValueOrder: (values: string[]) => void;

  // "Export box" overlay tool — draws an exact-mm rectangle in the 2D scene
  // for cropped PNG screenshots. Panel-local, not persisted.
  exportBoxEnabled: boolean;
  exportBoxWidthMm: number;
  exportBoxHeightMm: number;
  // Center position in CSS px relative to the scene container.
  // null = use viewport center on the next render.
  exportBoxCenterPx: { x: number; y: number } | null;
  setExportBoxEnabled: (enabled: boolean) => void;
  setExportBoxWidthMm: (mm: number) => void;
  setExportBoxHeightMm: (mm: number) => void;
  setExportBoxCenterPx: (pos: { x: number; y: number } | null) => void;

  // Two-gene coexpression mode.
  coexpressEnabled: boolean;
  selectedGene2: string | null;
  coexpressSwapped: boolean;
  gene2ScaleMin: number;
  gene2ScaleMax: number;
  setCoexpressEnabled: (on: boolean) => void;
  setSelectedGene2: (gene: string | null) => void;
  setCoexpressSwapped: (swapped: boolean) => void;
  setGene2ScaleMin: (min: number) => void;
  setGene2ScaleMax: (max: number) => void;

  // Per-sample 2D transforms (translate + rotate around centroid).
  transformColumn: string | null;
  activeSampleId: string | null;
  sampleTransforms: Map<string, SampleTransform>;
  transformVersion: number;
  setTransformColumn: (column: string | null) => void;
  setActiveSampleId: (id: string | null) => void;
  setSampleTransform: (sampleId: string, transform: SampleTransform) => void;
  clearSampleTransform: (sampleId: string) => void;
  clearAllSampleTransforms: () => void;

  reset: () => void;
}

const initialState = {
  mode: ["celltype"] as VisualizationMode[],
  panelMode: "celltype" as VisualizationMode,
  viewMode: "2D" as CellViewMode,
  selectedGene: null,
  geneEverywhere: false,
  hiddenCelltypes: new Set<string>(),
  geneScaleMin: VISUALIZATION_CONFIG.SCALE_BAR_DEFAULT_MIN,
  geneScaleMax: VISUALIZATION_CONFIG.SCALE_BAR_DEFAULT_MAX,
  selectedClusterColumn: null,
  selectedColumn: null,
  selectedCelltypes: new Set<string>(),
  numericalScaleMin: VISUALIZATION_CONFIG.SCALE_BAR_DEFAULT_MIN,
  numericalScaleMax: VISUALIZATION_CONFIG.SCALE_BAR_DEFAULT_MAX,
  celltypeSearchTerm: "",
  geneSearchTerm: "",
  selectedEmbedding: null,
  colorPalette: {},
  alphaScale: VISUALIZATION_CONFIG.POINT_BASE_ALPHA,
  sizeScale: 1.0,
  clusterVersion: 0,
  cameraResetSignal: 0,
  orbitSpeed: 1.0,
  orbitYBob: 0.3,
  targetFilterEnabled: false,
  scLayerVisible: true,
  // Distance filter radius is now absolute world units (microns for raw-coord
  // datasets), measured from controls.target.
  targetFilterRadius: 5000,
  deStatsVersion: 0,
  degTarget: null as string | null,
  degReference: null as string | null,
  degTargetAuto: true,
  degReferenceAuto: false,
  degSearchTerm: "",
  degSortKey: "log2FC" as "log2FC" | "meanIn" | "pctIn",
  degSortDesc: true,
  degPanelOpen: false,
  columnTypeOverrides: {} as Record<string, "categorical" | "numerical">,
  celltypePlayback: false,
  celltypePlaybackInterval: 1.0,
  celltypePlaybackSequence: [] as string[],
  sceneRotation: 0,
  flipX: false,
  flipY: false,
  selectedSizeMultiplier: VISUALIZATION_CONFIG.SELECTED_SIZE_MULTIPLIER as number,
  greyedOutSizeMultiplier: VISUALIZATION_CONFIG.GREYED_OUT_SIZE_MULTIPLIER as number,
  greyedOutAlpha: VISUALIZATION_CONFIG.GREYED_OUT_ALPHA as number,
  expressionAlphaMin: VISUALIZATION_CONFIG.EXPRESSION_ALPHA_MIN as number,
  expressionAlphaMax: VISUALIZATION_CONFIG.EXPRESSION_ALPHA_MAX as number,
  pointSizeMultiplierMin: VISUALIZATION_CONFIG.POINT_SIZE_MULTIPLIER_MIN as number,
  pointSizeMultiplierMax: VISUALIZATION_CONFIG.POINT_SIZE_MULTIPLIER_MAX as number,
  targetPx: VISUALIZATION_CONFIG.TARGET_PX_DEFAULT as number,
  pinnedTooltipColumns: new Set<string>(),
  sliderRanges: {} as Record<string, { min: number; max: number }>,
  colormap: "bwr",
  plotPanelOpen: false,
  secondaryColumn: null as string | null,
  selectedSecondaryValues: new Set<string>(),
  secondaryPaletteOverrides: {} as Record<string, string>,
  secondaryValueOrder: [] as string[],
  exportBoxEnabled: false,
  exportBoxWidthMm: 1,
  exportBoxHeightMm: 1,
  exportBoxCenterPx: null as { x: number; y: number } | null,
  transformColumn: null as string | null,
  activeSampleId: null as string | null,
  sampleTransforms: new Map<string, SampleTransform>(),
  transformVersion: 0,
  coexpressEnabled: false,
  selectedGene2: null as string | null,
  coexpressSwapped: false,
  gene2ScaleMin: VISUALIZATION_CONFIG.SCALE_BAR_DEFAULT_MIN,
  gene2ScaleMax: VISUALIZATION_CONFIG.SCALE_BAR_DEFAULT_MAX,
};

const updateModeArray = (
  currentMode: VisualizationMode[],
  add?: VisualizationMode,
  remove?: VisualizationMode,
): VisualizationMode[] => {
  let newMode = [...currentMode];

  if (remove && newMode.includes(remove)) {
    newMode = newMode.filter((m) => m !== remove);
  }

  if (add && !newMode.includes(add)) {
    newMode.push(add);
  }

  if (newMode.length === 0) {
    newMode = ["celltype"];
  }

  return newMode;
};

export function createVisualizationStoreInstance() {
  return createStore<VisualizationState>((set, get) => ({
    ...initialState,

    setMode: (mode) => {
      set({ mode });
    },

    setPanelMode: (mode) => {
      set({ panelMode: mode });
    },

    setViewMode: (mode) => {
      set({ viewMode: mode });
    },

    setSelectedGene: (gene) => {
      set((state) => {
        if (gene) {
          const newMode: VisualizationMode[] =
            state.selectedCelltypes.size > 0
              ? updateModeArray(state.mode, "gene")
              : (["gene"] as VisualizationMode[]);

          return { selectedGene: gene, mode: newMode };
        } else {
          const newMode: VisualizationMode[] = updateModeArray(
            state.mode,
            undefined,
            "gene",
          );

          return { selectedGene: gene, mode: newMode };
        }
      });
    },

    setGeneScaleMin: (min) => {
      set({ geneScaleMin: min });
    },

    setGeneScaleMax: (max) => {
      set({ geneScaleMax: max });
    },

    setNumericalScaleMin: (min) => {
      set({ numericalScaleMin: min });
    },

    setNumericalScaleMax: (max) => {
      set({ numericalScaleMax: max });
    },

    toggleCelltype: (celltype) => {
      set((state) => {
        const newCelltypes = new Set(state.selectedCelltypes);

        if (newCelltypes.has(celltype)) {
          newCelltypes.delete(celltype);
        } else {
          newCelltypes.add(celltype);
        }

        let newMode: VisualizationMode[] = [...state.mode];

        if (newCelltypes.size > 0) {
          if (!newMode.includes("celltype")) {
            newMode.push("celltype");
          }
        } else {
          if (!state.selectedGene) {
            newMode = newMode.filter((m) => m !== "celltype");
            if (newMode.length === 0) {
              newMode = ["celltype"];
            }
          }
        }

        // Toggling a celltype's selection resets its isolate-visibility.
        const newHidden = state.hiddenCelltypes.has(celltype)
          ? new Set([...state.hiddenCelltypes].filter((ct) => ct !== celltype))
          : state.hiddenCelltypes;

        return {
          selectedCelltypes: newCelltypes,
          mode: newMode,
          hiddenCelltypes: newHidden,
        };
      });
    },

    setCelltypes: (celltypes) => {
      set((state) => {
        let newMode: VisualizationMode[] = [...state.mode];

        if (celltypes.size > 0) {
          if (!newMode.includes("celltype")) {
            newMode.push("celltype");
          }
        } else {
          if (!state.selectedGene) {
            newMode = newMode.filter((m) => m !== "celltype");
            if (newMode.length === 0) {
              newMode = ["celltype"];
            }
          }
        }

        // Drop any hidden celltypes that are no longer selected.
        const newHidden = new Set(
          [...state.hiddenCelltypes].filter((ct) => celltypes.has(ct)),
        );

        return {
          selectedCelltypes: celltypes,
          mode: newMode,
          hiddenCelltypes: newHidden,
        };
      });
    },

    setClusterColumn: (column) => {
      set({ selectedClusterColumn: column });
    },

    setSelectedColumn: (column, isNumerical) => {
      set((state) => {
        const updates: Partial<VisualizationState> = {
          selectedColumn: column,
          selectedCelltypes: new Set<string>(),
          hiddenCelltypes: new Set<string>(),
          // Changing the primary column invalidates the grouping pairs.
          secondaryColumn: null,
          selectedSecondaryValues: new Set<string>(),
          secondaryPaletteOverrides: {},
          secondaryValueOrder: [],
        };

        if (isNumerical && column) {
          // Numerical column is mutually exclusive with gene expression colouring.
          updates.selectedGene = null;
          updates.mode = ["celltype"];
        } else {
          // Categorical column or cleared column: keep gene expression visible
          // if one is selected. selectedCelltypes was just cleared so combined
          // mode (`["gene","celltype"]`) doesn't apply here.
          updates.mode = state.selectedGene ? ["gene"] : ["celltype"];
        }

        return updates;
      });
    },

    setSelectedEmbedding: (embedding) => {
      set({ selectedEmbedding: embedding });
    },

    setColorPalette: (palette) => {
      set({ colorPalette: palette });
    },

    setAlphaScale: (alpha) => {
      set({ alphaScale: Math.max(0, Math.min(1, alpha)) });
    },

    setSizeScale: (size) => {
      set({ sizeScale: Math.max(0.1, size) });
    },

    setCelltypeSearchTerm: (value) => {
      set({ celltypeSearchTerm: value });
    },

    setGeneSearchTerm: (value) => {
      set({ geneSearchTerm: value });
    },

    incrementClusterVersion: () => {
      set((state) => ({ clusterVersion: state.clusterVersion + 1 }));
    },

    incrementDeStatsVersion: () => {
      set((state) => ({ deStatsVersion: state.deStatsVersion + 1 }));
    },

    resetCamera: () => {
      set((state) => ({ cameraResetSignal: state.cameraResetSignal + 1 }));
    },

    setTargetFilterEnabled: (enabled) => set({ targetFilterEnabled: enabled }),

    setScLayerVisible: (visible) => set({ scLayerVisible: visible }),
    setTargetFilterRadius: (r) => set({ targetFilterRadius: r }),

    setDegTarget: (target) => set({ degTarget: target }),
    setDegReference: (reference) => set({ degReference: reference }),
    setDegTargetAuto: (auto) => set({ degTargetAuto: auto }),
    setDegReferenceAuto: (auto) => set({ degReferenceAuto: auto }),
    setGeneEverywhere: (everywhere) => set({ geneEverywhere: everywhere }),
    toggleCelltypeVisibility: (celltype) =>
      set((state) => {
        const next = new Set(state.hiddenCelltypes);

        if (next.has(celltype)) next.delete(celltype);
        else next.add(celltype);

        return { hiddenCelltypes: next };
      }),
    soloCelltype: (celltype) =>
      set((state) => {
        const others = [...state.selectedCelltypes].filter(
          (ct) => ct !== celltype,
        );
        // Already soloed to this celltype → restore everything.
        const alreadySolo =
          others.length > 0 &&
          !state.hiddenCelltypes.has(celltype) &&
          others.every((ct) => state.hiddenCelltypes.has(ct));

        return {
          hiddenCelltypes: alreadySolo ? new Set<string>() : new Set(others),
        };
      }),
    setDegSearchTerm: (term) => set({ degSearchTerm: term }),
    setDegSortKey: (key) => set({ degSortKey: key }),
    setDegSortDesc: (desc) => set({ degSortDesc: desc }),
    setDegPanelOpen: (open) => set({ degPanelOpen: open }),

    toggleColumnType: (column, currentType) => {
      set((state) => {
        const newType: "categorical" | "numerical" =
          currentType === "categorical" ? "numerical" : "categorical";
        const newOverrides: Record<string, "categorical" | "numerical"> = {
          ...state.columnTypeOverrides,
          [column]: newType,
        };

        const updates: Partial<VisualizationState> = {
          columnTypeOverrides: newOverrides,
          clusterVersion: state.clusterVersion + 1,
        };

        // If toggling the currently selected column to numerical, clear gene/celltypes
        if (state.selectedColumn === column && newType === "numerical") {
          updates.selectedGene = null;
          updates.selectedCelltypes = new Set<string>();
          updates.hiddenCelltypes = new Set<string>();
          updates.mode = ["celltype"];
        }

        return updates;
      });
    },

    setColumnTypeOverrides: (overrides) => {
      set({ columnTypeOverrides: overrides });
    },

    setCelltypePlayback: (playing) => {
      set({ celltypePlayback: playing });
    },

    setCelltypePlaybackInterval: (interval) => {
      set({ celltypePlaybackInterval: interval });
    },

    setCelltypePlaybackSequence: (sequence) => {
      set({ celltypePlaybackSequence: sequence });
    },

    setSceneRotation: (degrees) => {
      set({ sceneRotation: degrees });
    },

    setFlipX: (flip) => {
      set({ flipX: flip });
    },

    setFlipY: (flip) => {
      set({ flipY: flip });
    },

    setAdvancedViz: (key, value) => {
      set({ [key]: value } as any);
    },

    togglePinnedTooltipColumn: (column) => {
      set((state) => {
        const next = new Set(state.pinnedTooltipColumns);
        if (next.has(column)) {
          next.delete(column);
        } else {
          next.add(column);
        }
        return { pinnedTooltipColumns: next };
      });
    },

    setSliderRange: (key, min, max) => {
      if (!(min < max)) return;
      set((state) => ({
        sliderRanges: { ...state.sliderRanges, [key]: { min, max } },
      }));
    },

    clearSliderRanges: (keys) => {
      set((state) => {
        const next = { ...state.sliderRanges };
        for (const k of keys) delete next[k];
        return { sliderRanges: next };
      });
    },

    setColormap: (name) => {
      set({ colormap: name });
    },

    setPlotPanelOpen: (open) => {
      set({ plotPanelOpen: open });
    },

    setSecondaryColumn: (column) => {
      set((state) => {
        // No-op when picking the primary column as secondary (UI also blocks
        // this, but defend against direct calls).
        if (column && column === state.selectedColumn) return {} as Partial<VisualizationState>;
        return {
          secondaryColumn: column,
          selectedSecondaryValues: new Set<string>(),
          secondaryPaletteOverrides: {},
          secondaryValueOrder: [],
        };
      });
    },

    setSelectedSecondaryValues: (values) => {
      set({ selectedSecondaryValues: new Set(values) });
    },

    toggleSecondaryValue: (value) => {
      set((state) => {
        const next = new Set(state.selectedSecondaryValues);
        if (next.has(value)) next.delete(value);
        else next.add(value);
        return { selectedSecondaryValues: next };
      });
    },

    setSecondaryPaletteOverride: (value, color) => {
      set((state) => ({
        secondaryPaletteOverrides: {
          ...state.secondaryPaletteOverrides,
          [value]: color,
        },
      }));
    },

    setSecondaryValueOrder: (values) => {
      set({ secondaryValueOrder: [...values] });
    },

    setExportBoxEnabled: (enabled) => set({ exportBoxEnabled: enabled }),
    setExportBoxWidthMm: (mm) =>
      set({ exportBoxWidthMm: Math.max(0.01, mm) }),
    setExportBoxHeightMm: (mm) =>
      set({ exportBoxHeightMm: Math.max(0.01, mm) }),
    setExportBoxCenterPx: (pos) => set({ exportBoxCenterPx: pos }),

    setCoexpressEnabled: (on) => {
      set({ coexpressEnabled: on });
    },
    setSelectedGene2: (gene) => {
      set({ selectedGene2: gene });
    },
    setCoexpressSwapped: (swapped) => {
      set({ coexpressSwapped: swapped });
    },
    setGene2ScaleMin: (min) => {
      set({ gene2ScaleMin: min });
    },
    setGene2ScaleMax: (max) => {
      set({ gene2ScaleMax: max });
    },

    setTransformColumn: (column) => {
      // Switching the sample column invalidates the active selection;
      // existing transforms are kept (keyed by sample value, not column).
      set({ transformColumn: column, activeSampleId: null });
    },

    setActiveSampleId: (id) => {
      set({ activeSampleId: id });
    },

    setSampleTransform: (sampleId, transform) => {
      set((state) => {
        const next = new Map(state.sampleTransforms);
        next.set(sampleId, transform);
        return {
          sampleTransforms: next,
          transformVersion: state.transformVersion + 1,
        };
      });
    },

    clearSampleTransform: (sampleId) => {
      set((state) => {
        if (!state.sampleTransforms.has(sampleId)) return {} as Partial<VisualizationState>;
        const next = new Map(state.sampleTransforms);
        next.delete(sampleId);
        return {
          sampleTransforms: next,
          transformVersion: state.transformVersion + 1,
        };
      });
    },

    clearAllSampleTransforms: () => {
      set((state) => ({
        sampleTransforms: new Map<string, SampleTransform>(),
        transformVersion: state.transformVersion + 1,
      }));
    },

    reset: () => {
      set(initialState);
    },
  }));
}

export type VisualizationStore = ReturnType<
  typeof createVisualizationStoreInstance
>;
