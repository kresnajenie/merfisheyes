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
  // Selected genes with their visualization properties (visible in scene)
  selectedGenes: Map<string, GeneVisualization>;

  // Genes that appear in legend panel (may be hidden from scene)
  selectedGenesLegend: Set<string>;

  // Cache of gene visualization data for all legend genes (even if hidden)
  geneDataCache: Map<string, GeneVisualization>;

  // Global scale multiplier (affects all genes)
  globalScale: number;

  // Camera view mode
  viewMode: ViewMode;

  // Actions
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

export const useSingleMoleculeVisualizationStore =
  create<SingleMoleculeVisualizationState>((set, get) => ({
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

        // Also add to legend and cache
        newSelectedGenesLegend.add(gene);
        newGeneDataCache.set(gene, geneViz);

        return {
          selectedGenes: newSelectedGenes,
          selectedGenesLegend: newSelectedGenesLegend,
          geneDataCache: newGeneDataCache
        };
      }),

    removeGene: (gene: string) =>
      set((state) => {
        const newSelectedGenes = new Map(state.selectedGenes);
        const newSelectedGenesLegend = new Set(state.selectedGenesLegend);
        const newGeneDataCache = new Map(state.geneDataCache);

        // Remove from visibility, legend, and cache
        newSelectedGenes.delete(gene);
        newSelectedGenesLegend.delete(gene);
        newGeneDataCache.delete(gene);

        return {
          selectedGenes: newSelectedGenes,
          selectedGenesLegend: newSelectedGenesLegend,
          geneDataCache: newGeneDataCache
        };
      }),

    toggleGeneVisibility: (gene: string) =>
      set((state) => {
        const newSelectedGenes = new Map(state.selectedGenes);

        if (newSelectedGenes.has(gene)) {
          // Gene is visible, hide it (but keep in legend and cache)
          newSelectedGenes.delete(gene);
        } else {
          // Gene is hidden, show it (restore from cache)
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
        const geneViz = newSelectedGenes.get(gene) || newGeneDataCache.get(gene);

        if (geneViz) {
          const updatedGeneViz = { ...geneViz, color };
          // Update in both visible genes (if present) and cache
          if (newSelectedGenes.has(gene)) {
            newSelectedGenes.set(gene, updatedGeneViz);
          }
          newGeneDataCache.set(gene, updatedGeneViz);
        }

        return { selectedGenes: newSelectedGenes, geneDataCache: newGeneDataCache };
      }),

    setGeneLocalScale: (gene: string, scale: number) =>
      set((state) => {
        const newSelectedGenes = new Map(state.selectedGenes);
        const newGeneDataCache = new Map(state.geneDataCache);
        const geneViz = newSelectedGenes.get(gene) || newGeneDataCache.get(gene);

        if (geneViz) {
          const updatedGeneViz = { ...geneViz, localScale: scale };
          // Update in both visible genes (if present) and cache
          if (newSelectedGenes.has(gene)) {
            newSelectedGenes.set(gene, updatedGeneViz);
          }
          newGeneDataCache.set(gene, updatedGeneViz);
        }

        return { selectedGenes: newSelectedGenes, geneDataCache: newGeneDataCache };
      }),

    setGlobalScale: (scale: number) => set({ globalScale: scale }),

    setViewMode: (mode: ViewMode) => set({ viewMode: mode }),

    clearGenes: () => set({ selectedGenes: new Map(), selectedGenesLegend: new Set(), geneDataCache: new Map() }),

    loadFromLocalStorage: (datasetId: string) => {
      if (typeof window === "undefined") return;

      try {
        const storageKey = `sm_gene_visibility_${datasetId}`;
        const stored = localStorage.getItem(storageKey);

        if (stored) {
          const { visibleGenes, legendGenes, geneData } = JSON.parse(stored);

          // Reconstruct selectedGenes Map (visible genes)
          const newSelectedGenes = new Map<string, GeneVisualization>();
          visibleGenes.forEach((gene: string) => {
            if (geneData[gene]) {
              newSelectedGenes.set(gene, geneData[gene]);
            }
          });

          // Reconstruct selectedGenesLegend Set
          const newSelectedGenesLegend = new Set<string>(legendGenes);

          // Reconstruct geneDataCache Map (all legend genes)
          const newGeneDataCache = new Map<string, GeneVisualization>();
          legendGenes.forEach((gene: string) => {
            if (geneData[gene]) {
              newGeneDataCache.set(gene, geneData[gene]);
            }
          });

          set({
            selectedGenes: newSelectedGenes,
            selectedGenesLegend: newSelectedGenesLegend,
            geneDataCache: newGeneDataCache
          });

          console.log(`[SingleMoleculeVisualizationStore] Loaded visibility state from localStorage for dataset: ${datasetId}`);
        }
      } catch (error) {
        console.warn(`[SingleMoleculeVisualizationStore] Failed to load visibility state from localStorage:`, error);
      }
    },

    saveToLocalStorage: (datasetId: string) => {
      if (typeof window === "undefined") return;

      try {
        const state = get();
        const storageKey = `sm_gene_visibility_${datasetId}`;

        // Convert Map and Set to serializable format
        const visibleGenes = Array.from(state.selectedGenes.keys());
        const legendGenes = Array.from(state.selectedGenesLegend);

        // Store full gene data from cache (includes hidden genes)
        const geneData: Record<string, GeneVisualization> = {};
        state.geneDataCache.forEach((geneViz, gene) => {
          geneData[gene] = geneViz;
        });

        const data = {
          visibleGenes,
          legendGenes,
          geneData
        };

        localStorage.setItem(storageKey, JSON.stringify(data));
        console.log(`[SingleMoleculeVisualizationStore] Saved visibility state to localStorage for dataset: ${datasetId}`);
      } catch (error) {
        console.warn(`[SingleMoleculeVisualizationStore] Failed to save visibility state to localStorage:`, error);
      }
    }
  }));
