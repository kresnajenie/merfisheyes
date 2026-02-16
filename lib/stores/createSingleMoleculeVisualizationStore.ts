import { createStore } from "zustand";

import {
  getColorForSlot,
  findLowestAvailableSlot,
} from "@/lib/utils/gene-color-palette";

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
  geneColorSlots: Map<string, number>;
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
    geneColorSlots: new Map(),
    globalScale: 1.0,
    viewMode: "2D",

    addGene: (gene: string, color?: string, localScale?: number) =>
      set((state) => {
        const newSelectedGenes = new Map(state.selectedGenes);
        const newSelectedGenesLegend = new Set(state.selectedGenesLegend);
        const newGeneDataCache = new Map(state.geneDataCache);
        const newGeneColorSlots = new Map(state.geneColorSlots);

        // Assign color by lowest-available-slot if no explicit color provided
        let assignedColor = color;

        if (!assignedColor) {
          let slot = newGeneColorSlots.get(gene);

          if (slot === undefined) {
            const usedSlots = new Set(newGeneColorSlots.values());

            slot = findLowestAvailableSlot(usedSlots);
            newGeneColorSlots.set(gene, slot);
          }
          assignedColor = getColorForSlot(slot);
        } else if (!newGeneColorSlots.has(gene)) {
          // Even with explicit color (URL restore), assign a slot
          const usedSlots = new Set(newGeneColorSlots.values());
          const slot = findLowestAvailableSlot(usedSlots);

          newGeneColorSlots.set(gene, slot);
        }

        const geneViz = {
          gene,
          color: assignedColor,
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
          geneColorSlots: newGeneColorSlots,
        };
      }),

    removeGene: (gene: string) =>
      set((state) => {
        const newSelectedGenes = new Map(state.selectedGenes);
        const newSelectedGenesLegend = new Set(state.selectedGenesLegend);
        const newGeneDataCache = new Map(state.geneDataCache);
        const newGeneColorSlots = new Map(state.geneColorSlots);

        newSelectedGenes.delete(gene);
        newSelectedGenesLegend.delete(gene);
        newGeneDataCache.delete(gene);
        newGeneColorSlots.delete(gene);

        return {
          selectedGenes: newSelectedGenes,
          selectedGenesLegend: newSelectedGenesLegend,
          geneDataCache: newGeneDataCache,
          geneColorSlots: newGeneColorSlots,
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
        geneColorSlots: new Map(),
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

          // Rebuild color slots: assign slots in legend order
          const newGeneColorSlots = new Map<string, number>();

          legendGenes.forEach((gene: string, index: number) => {
            newGeneColorSlots.set(gene, index);
          });

          set({
            selectedGenes: newSelectedGenes,
            selectedGenesLegend: newSelectedGenesLegend,
            geneDataCache: newGeneDataCache,
            geneColorSlots: newGeneColorSlots,
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
