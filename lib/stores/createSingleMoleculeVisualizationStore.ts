import { createStore } from "zustand";

function generateBrightColor(): string {
  const hue = Math.random() * 360;
  const saturation = 70 + Math.random() * 30;
  const lightness = 50 + Math.random() * 20;

  return `hsl(${hue}, ${saturation}%, ${lightness}%)`;
}

export interface GeneVisualization {
  gene: string;
  color: string;
  localScale: number;
}

export type ViewMode = "2D" | "3D";

export interface SingleMoleculeVisualizationState {
  selectedGenes: Map<string, GeneVisualization>;
  selectedGenesLegend: Set<string>;
  geneDataCache: Map<string, GeneVisualization>;
  globalScale: number;
  viewMode: ViewMode;

  addGene: (gene: string, color?: string, localScale?: number) => void;
  removeGene: (gene: string) => void;
  toggleGeneVisibility: (gene: string) => void;
  setGeneColor: (gene: string, color: string) => void;
  setGeneLocalScale: (gene: string, scale: number) => void;
  setGlobalScale: (scale: number) => void;
  setViewMode: (mode: ViewMode) => void;
  clearGenes: () => void;
  loadFromLocalStorage: (datasetId: string) => void;
  saveToLocalStorage: (datasetId: string) => void;
}

export function createSingleMoleculeVisualizationStoreInstance() {
  return createStore<SingleMoleculeVisualizationState>((set, get) => ({
    selectedGenes: new Map(),
    selectedGenesLegend: new Set(),
    geneDataCache: new Map(),
    globalScale: 1.0,
    viewMode: "2D",

    addGene: (gene: string, color?: string, localScale?: number) =>
      set((state) => {
        const newSelectedGenes = new Map(state.selectedGenes);
        const newSelectedGenesLegend = new Set(state.selectedGenesLegend);
        const newGeneDataCache = new Map(state.geneDataCache);

        const geneViz = {
          gene,
          color: color || generateBrightColor(),
          localScale: localScale || 1.0,
        };

        if (!newSelectedGenes.has(gene)) {
          newSelectedGenes.set(gene, geneViz);
        }

        newSelectedGenesLegend.add(gene);
        newGeneDataCache.set(gene, geneViz);

        return {
          selectedGenes: newSelectedGenes,
          selectedGenesLegend: newSelectedGenesLegend,
          geneDataCache: newGeneDataCache,
        };
      }),

    removeGene: (gene: string) =>
      set((state) => {
        const newSelectedGenes = new Map(state.selectedGenes);
        const newSelectedGenesLegend = new Set(state.selectedGenesLegend);
        const newGeneDataCache = new Map(state.geneDataCache);

        newSelectedGenes.delete(gene);
        newSelectedGenesLegend.delete(gene);
        newGeneDataCache.delete(gene);

        return {
          selectedGenes: newSelectedGenes,
          selectedGenesLegend: newSelectedGenesLegend,
          geneDataCache: newGeneDataCache,
        };
      }),

    toggleGeneVisibility: (gene: string) =>
      set((state) => {
        const newSelectedGenes = new Map(state.selectedGenes);

        if (newSelectedGenes.has(gene)) {
          newSelectedGenes.delete(gene);
        } else {
          const cachedData = state.geneDataCache.get(gene);

          if (cachedData) {
            newSelectedGenes.set(gene, cachedData);
          }
        }

        return { selectedGenes: newSelectedGenes };
      }),

    setGeneColor: (gene: string, color: string) =>
      set((state) => {
        const newSelectedGenes = new Map(state.selectedGenes);
        const newGeneDataCache = new Map(state.geneDataCache);
        const geneViz =
          newSelectedGenes.get(gene) || newGeneDataCache.get(gene);

        if (geneViz) {
          const updatedGeneViz = { ...geneViz, color };

          if (newSelectedGenes.has(gene)) {
            newSelectedGenes.set(gene, updatedGeneViz);
          }
          newGeneDataCache.set(gene, updatedGeneViz);
        }

        return {
          selectedGenes: newSelectedGenes,
          geneDataCache: newGeneDataCache,
        };
      }),

    setGeneLocalScale: (gene: string, scale: number) =>
      set((state) => {
        const newSelectedGenes = new Map(state.selectedGenes);
        const newGeneDataCache = new Map(state.geneDataCache);
        const geneViz =
          newSelectedGenes.get(gene) || newGeneDataCache.get(gene);

        if (geneViz) {
          const updatedGeneViz = { ...geneViz, localScale: scale };

          if (newSelectedGenes.has(gene)) {
            newSelectedGenes.set(gene, updatedGeneViz);
          }
          newGeneDataCache.set(gene, updatedGeneViz);
        }

        return {
          selectedGenes: newSelectedGenes,
          geneDataCache: newGeneDataCache,
        };
      }),

    setGlobalScale: (scale: number) => set({ globalScale: scale }),

    setViewMode: (mode: ViewMode) => set({ viewMode: mode }),

    clearGenes: () =>
      set({
        selectedGenes: new Map(),
        selectedGenesLegend: new Set(),
        geneDataCache: new Map(),
      }),

    loadFromLocalStorage: (datasetId: string) => {
      if (typeof window === "undefined") return;

      try {
        const storageKey = `sm_gene_visibility_${datasetId}`;
        const stored = localStorage.getItem(storageKey);

        if (stored) {
          const { visibleGenes, legendGenes, geneData } = JSON.parse(stored);

          const newSelectedGenes = new Map<string, GeneVisualization>();

          visibleGenes.forEach((gene: string) => {
            if (geneData[gene]) {
              newSelectedGenes.set(gene, geneData[gene]);
            }
          });

          const newSelectedGenesLegend = new Set<string>(legendGenes);

          const newGeneDataCache = new Map<string, GeneVisualization>();

          legendGenes.forEach((gene: string) => {
            if (geneData[gene]) {
              newGeneDataCache.set(gene, geneData[gene]);
            }
          });

          set({
            selectedGenes: newSelectedGenes,
            selectedGenesLegend: newSelectedGenesLegend,
            geneDataCache: newGeneDataCache,
          });
        }
      } catch (error) {
        console.warn(
          `[SingleMoleculeVisualizationStore] Failed to load from localStorage:`,
          error,
        );
      }
    },

    saveToLocalStorage: (datasetId: string) => {
      if (typeof window === "undefined") return;

      try {
        const state = get();
        const storageKey = `sm_gene_visibility_${datasetId}`;

        const visibleGenes = Array.from(state.selectedGenes.keys());
        const legendGenes = Array.from(state.selectedGenesLegend);

        const geneData: Record<string, GeneVisualization> = {};

        state.geneDataCache.forEach((geneViz, gene) => {
          geneData[gene] = geneViz;
        });

        const data = { visibleGenes, legendGenes, geneData };

        localStorage.setItem(storageKey, JSON.stringify(data));
      } catch (error) {
        console.warn(
          `[SingleMoleculeVisualizationStore] Failed to save to localStorage:`,
          error,
        );
      }
    },
  }));
}

export type SingleMoleculeVisualizationStore = ReturnType<
  typeof createSingleMoleculeVisualizationStoreInstance
>;
