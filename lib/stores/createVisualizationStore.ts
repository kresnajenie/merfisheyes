import { createStore } from "zustand";

import { VISUALIZATION_CONFIG } from "../config/visualization.config";

export type VisualizationMode = "celltype" | "gene";

export interface VisualizationState {
  mode: VisualizationMode[];
  panelMode: VisualizationMode;
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

  setMode: (mode: VisualizationMode[]) => void;
  setPanelMode: (mode: VisualizationMode) => void;
  setSelectedGene: (gene: string | null) => void;
  setGeneScaleMin: (min: number) => void;
  setGeneScaleMax: (max: number) => void;
  setNumericalScaleMin: (min: number) => void;
  setNumericalScaleMax: (max: number) => void;
  toggleCelltype: (celltype: string) => void;
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
  reset: () => void;
}

const initialState = {
  mode: ["celltype"] as VisualizationMode[],
  panelMode: "celltype" as VisualizationMode,
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
          updates.selectedGene = null;
          updates.mode = ["celltype"];
        } else if (!column) {
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

    reset: () => {
      set(initialState);
    },
  }));
}

export type VisualizationStore = ReturnType<
  typeof createVisualizationStoreInstance
>;
