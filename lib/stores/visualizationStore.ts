import { create } from "zustand";

export type VisualizationMode = "celltype" | "gene";

interface VisualizationState {
  // Visualization mode
  mode: VisualizationMode;

  // Gene-specific settings
  selectedGene: string | null;

  // Celltype/cluster settings
  selectedClusterColumn: string | null;
  selectedColumn: string | null; // Currently active obs column for celltype mode
  selectedCelltypes: Set<string>; // Selected celltypes from the list

  // Embedding settings
  selectedEmbedding: string | null;

  // Visual properties
  colorPalette: Record<string, string>;
  alphaScale: number; // 0-1
  sizeScale: number; // multiplier

  // Actions
  setMode: (mode: VisualizationMode) => void;
  setSelectedGene: (gene: string | null) => void;
  toggleCelltype: (celltype: string) => void;
  setClusterColumn: (column: string | null) => void;
  setSelectedColumn: (column: string | null) => void;
  setSelectedEmbedding: (embedding: string | null) => void;
  setColorPalette: (palette: Record<string, string>) => void;
  setAlphaScale: (alpha: number) => void;
  setSizeScale: (size: number) => void;
  reset: () => void;
}

const initialState = {
  mode: "celltype" as VisualizationMode,
  selectedGene: null,
  selectedClusterColumn: null,
  selectedColumn: null,
  selectedCelltypes: new Set<string>(),
  selectedEmbedding: null,
  colorPalette: {},
  alphaScale: 1.0,
  sizeScale: 1.0,
};

export const useVisualizationStore = create<VisualizationState>((set) => ({
  ...initialState,

  setMode: (mode) => {
    console.log("Visualization mode changed to:", mode);
    set({ mode });
  },

  setSelectedGene: (gene) => {
    console.log("Selected gene:", gene);
    set({ selectedGene: gene });
  },

  toggleCelltype: (celltype) => {
    set((state) => {
      const newCelltypes = new Set(state.selectedCelltypes);
      if (newCelltypes.has(celltype)) {
        newCelltypes.delete(celltype);
      } else {
        newCelltypes.add(celltype);
      }
      console.log("Toggled celltype:", celltype, "Selected celltypes:", Array.from(newCelltypes));
      return { selectedCelltypes: newCelltypes };
    });
  },

  setClusterColumn: (column) => {
    console.log("Selected cluster column:", column);
    set({ selectedClusterColumn: column });
  },

  setSelectedColumn: (column) => {
    console.log("Selected column:", column);
    set({
      selectedColumn: column,
      selectedCelltypes: new Set<string>(),
      selectedGene: null,
    });
  },

  setSelectedEmbedding: (embedding) => {
    console.log("Selected embedding:", embedding);
    set({ selectedEmbedding: embedding });
  },

  setColorPalette: (palette) => {
    set({ colorPalette: palette });
  },

  setAlphaScale: (alpha) => {
    set({ alphaScale: Math.max(0, Math.min(1, alpha)) }); // Clamp 0-1
  },

  setSizeScale: (size) => {
    set({ sizeScale: Math.max(0.1, size) }); // Minimum 0.1
  },

  reset: () => {
    console.log("Resetting visualization settings");
    set(initialState);
  },
}));
