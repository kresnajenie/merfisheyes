import { create } from "zustand";

export type VisualizationMode = "celltype" | "gene";

interface VisualizationState {
  // Visualization mode (what's actually being rendered) - can be multiple modes
  mode: VisualizationMode[];

  // Panel mode (which panel is open)
  panelMode: VisualizationMode;

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
  setMode: (mode: VisualizationMode[]) => void;
  setPanelMode: (mode: VisualizationMode) => void;
  setSelectedGene: (gene: string | null) => void;
  toggleCelltype: (celltype: string) => void;
  setClusterColumn: (column: string | null) => void;
  setSelectedColumn: (column: string | null, isNumerical?: boolean) => void;
  setSelectedEmbedding: (embedding: string | null) => void;
  setColorPalette: (palette: Record<string, string>) => void;
  setAlphaScale: (alpha: number) => void;
  setSizeScale: (size: number) => void;
  reset: () => void;
}

const initialState = {
  mode: ["celltype"] as VisualizationMode[],
  panelMode: "celltype" as VisualizationMode,
  selectedGene: null,
  selectedClusterColumn: null,
  selectedColumn: null,
  selectedCelltypes: new Set<string>(),
  selectedEmbedding: null,
  colorPalette: {},
  alphaScale: 1.0,
  sizeScale: 1.0,
};

// Helper function to update mode array
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

  // Ensure mode array is never empty, default to celltype
  if (newMode.length === 0) {
    newMode = ["celltype"];
  }

  return newMode;
};

export const useVisualizationStore = create<VisualizationState>((set, get) => ({
  ...initialState,

  setMode: (mode) => {
    console.log("Visualization mode changed to:", mode);
    set({ mode });
  },

  setPanelMode: (mode) => {
    console.log("Panel mode changed to:", mode);
    set({ panelMode: mode });
  },

  setSelectedGene: (gene) => {
    console.log("Selected gene:", gene);
    set((state) => {
      // When selecting a gene, keep selectedColumn but update mode
      if (gene) {
        // Add gene to mode, keep celltype if celltypes are selected
        const newMode: VisualizationMode[] =
          state.selectedCelltypes.size > 0
            ? updateModeArray(state.mode, "gene")
            : (["gene"] as VisualizationMode[]);

        console.log("Gene selected - new mode:", newMode);
        return { selectedGene: gene, mode: newMode };
      } else {
        // When clearing gene, remove from mode
        const newMode: VisualizationMode[] = updateModeArray(
          state.mode,
          undefined,
          "gene",
        );

        console.log("Gene cleared - new mode:", newMode);
        return { selectedGene: gene, mode: newMode };
      }
    });
  },

  toggleCelltype: (celltype) => {
    set((state) => {
      const newCelltypes = new Set(state.selectedCelltypes);

      if (newCelltypes.has(celltype)) {
        newCelltypes.delete(celltype);
      } else {
        newCelltypes.add(celltype);
      }
      console.log(
        "Toggled celltype:",
        celltype,
        "Selected celltypes:",
        Array.from(newCelltypes),
      );

      // Update mode based on selections
      let newMode: VisualizationMode[] = [...state.mode];

      if (newCelltypes.size > 0) {
        // Add celltype to mode if celltypes are selected
        if (!newMode.includes("celltype")) {
          newMode.push("celltype");
        }
      } else {
        // Remove celltype from mode if no celltypes selected and gene is not selected
        if (!state.selectedGene) {
          newMode = newMode.filter((m) => m !== "celltype");
          if (newMode.length === 0) {
            newMode = ["celltype"];
          }
        }
      }

      console.log("After toggle - new mode:", newMode);
      return { selectedCelltypes: newCelltypes, mode: newMode };
    });
  },

  setClusterColumn: (column) => {
    console.log("Selected cluster column:", column);
    set({ selectedClusterColumn: column });
  },

  setSelectedColumn: (column, isNumerical) => {
    console.log("Selected column:", column, "isNumerical:", isNumerical);
    set((state) => {
      const updates: Partial<VisualizationState> = {
        selectedColumn: column,
        selectedCelltypes: new Set<string>(),
      };

      // If selecting a numerical column, clear gene and set mode to celltype only
      if (isNumerical && column) {
        console.log("Numerical column selected - clearing gene");
        updates.selectedGene = null;
        updates.mode = ["celltype"];
      } else if (!column) {
        // If clearing column, reset mode
        updates.mode = state.selectedGene ? ["gene"] : ["celltype"];
      }

      return updates;
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
