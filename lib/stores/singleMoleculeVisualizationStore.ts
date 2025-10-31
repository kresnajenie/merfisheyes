import { create } from "zustand";

// Generate a random bright color for dark background
function generateBrightColor(): string {
  const hue = Math.random() * 360;
  const saturation = 70 + Math.random() * 30; // 70-100%
  const lightness = 50 + Math.random() * 20; // 50-70%

  return `hsl(${hue}, ${saturation}%, ${lightness}%)`;
}

export interface GeneVisualization {
  gene: string;
  color: string;
  localScale: number; // Local scale multiplier for this gene
}

export type ViewMode = "2D" | "3D";

interface SingleMoleculeVisualizationState {
  // Selected genes with their visualization properties
  selectedGenes: Map<string, GeneVisualization>;

  // Global scale multiplier (affects all genes)
  globalScale: number;

  // Camera view mode
  viewMode: ViewMode;

  // Actions
  addGene: (gene: string, color?: string, localScale?: number) => void;
  removeGene: (gene: string) => void;
  setGeneColor: (gene: string, color: string) => void;
  setGeneLocalScale: (gene: string, scale: number) => void;
  setGlobalScale: (scale: number) => void;
  setViewMode: (mode: ViewMode) => void;
  clearGenes: () => void;
}

export const useSingleMoleculeVisualizationStore =
  create<SingleMoleculeVisualizationState>((set) => ({
    selectedGenes: new Map(),
    globalScale: 1.0,
    viewMode: "2D",

    addGene: (gene: string, color?: string, localScale?: number) =>
      set((state) => {
        const newSelectedGenes = new Map(state.selectedGenes);

        if (!newSelectedGenes.has(gene)) {
          newSelectedGenes.set(gene, {
            gene,
            color: color || generateBrightColor(),
            localScale: localScale || 1.0,
          });
        }

        return { selectedGenes: newSelectedGenes };
      }),

    removeGene: (gene: string) =>
      set((state) => {
        const newSelectedGenes = new Map(state.selectedGenes);

        newSelectedGenes.delete(gene);

        return { selectedGenes: newSelectedGenes };
      }),

    setGeneColor: (gene: string, color: string) =>
      set((state) => {
        const newSelectedGenes = new Map(state.selectedGenes);
        const geneViz = newSelectedGenes.get(gene);

        if (geneViz) {
          newSelectedGenes.set(gene, { ...geneViz, color });
        }

        return { selectedGenes: newSelectedGenes };
      }),

    setGeneLocalScale: (gene: string, scale: number) =>
      set((state) => {
        const newSelectedGenes = new Map(state.selectedGenes);
        const geneViz = newSelectedGenes.get(gene);

        if (geneViz) {
          newSelectedGenes.set(gene, { ...geneViz, localScale: scale });
        }

        return { selectedGenes: newSelectedGenes };
      }),

    setGlobalScale: (scale: number) => set({ globalScale: scale }),

    setViewMode: (mode: ViewMode) => set({ viewMode: mode }),

    clearGenes: () => set({ selectedGenes: new Map() }),
  }));
