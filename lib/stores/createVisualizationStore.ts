import { createStore } from "zustand";

import { VISUALIZATION_CONFIG } from "../config/visualization.config";

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
  numericalScaleMin: number;
  numericalScaleMax: number;
  celltypeSearchTerm: string;
  geneSearchTerm: string;
  selectedEmbedding: string | null;
  colorPalette: Record<string, string>;
  alphaScale: number;
  sizeScale: number;
  clusterVersion: number;
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
  reset: () => void;
}

const initialState = {
  mode: ["celltype"] as VisualizationMode[],
  panelMode: "celltype" as VisualizationMode,
  viewMode: "2D" as CellViewMode,
  selectedGene: null,
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

        return { selectedCelltypes: newCelltypes, mode: newMode };
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

        return { selectedCelltypes: celltypes, mode: newMode };
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

    reset: () => {
      set(initialState);
    },
  }));
}

export type VisualizationStore = ReturnType<
  typeof createVisualizationStoreInstance
>;
