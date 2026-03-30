import { create } from "zustand";

import {
  getColorForSlot,
  findLowestAvailableSlot,
} from "@/lib/utils/gene-color-palette";

export type MoleculeShape = "circle" | "square";

export interface GeneVisualization {
  gene: string;
  color: string;
  localScale: number;
  showAssigned: boolean;
  showUnassigned: boolean;
  assignedShape: MoleculeShape;
  unassignedShape: MoleculeShape;
  unassignedColor: string;
  unassignedLocalScale: number;
}

export type ViewMode = "2D" | "3D";

interface SingleMoleculeVisualizationState {
  // Selected genes with their visualization properties (visible in scene)
  selectedGenes: Map<string, GeneVisualization>;

  // Genes that appear in legend panel (may be hidden from scene)
  selectedGenesLegend: Set<string>;

  // Cache of gene visualization data for all legend genes (even if hidden)
  geneDataCache: Map<string, GeneVisualization>;

  // Maps gene name to its color slot index
  geneColorSlots: Map<string, number>;

  // Global scale multiplier (affects all genes)
  globalScale: number;

  // Camera view mode
  viewMode: ViewMode;

  // Assigned/unassigned molecule visibility
  showAssigned: boolean;
  showUnassigned: boolean;

  // Actions
  addGene: (gene: string, color?: string, localScale?: number) => void;
  removeGene: (gene: string) => void;
  toggleGeneVisibility: (gene: string) => void;
  setGeneColor: (gene: string, color: string) => void;
  setGeneLocalScale: (gene: string, scale: number) => void;
  setGlobalScale: (scale: number) => void;
  setViewMode: (mode: ViewMode) => void;
  setShowAssigned: (show: boolean) => void;
  setShowUnassigned: (show: boolean) => void;
  setGeneShowAssigned: (gene: string, show: boolean) => void;
  setGeneShowUnassigned: (gene: string, show: boolean) => void;
  setGeneAssignedShape: (gene: string, shape: MoleculeShape) => void;
  setGeneUnassignedShape: (gene: string, shape: MoleculeShape) => void;
  setGeneUnassignedColor: (gene: string, color: string) => void;
  setGeneUnassignedLocalScale: (gene: string, scale: number) => void;
  clearGenes: () => void;
}

export const useSingleMoleculeVisualizationStore =
  create<SingleMoleculeVisualizationState>((set, get) => ({
    selectedGenes: new Map(),
    selectedGenesLegend: new Set(),
    geneDataCache: new Map(),
    geneColorSlots: new Map(),
    globalScale: 1.0,
    viewMode: "2D",
    showAssigned: true,
    showUnassigned: true,

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

        const geneViz: GeneVisualization = {
          gene,
          color: assignedColor,
          localScale: localScale || 1.0,
          showAssigned: true,
          showUnassigned: true,
          assignedShape: "circle",
          unassignedShape: "square",
          unassignedColor: assignedColor,
          unassignedLocalScale: localScale || 1.0,
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

        // Remove from visibility, legend, cache, and color slot
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
        const geneViz =
          newSelectedGenes.get(gene) || newGeneDataCache.get(gene);

        if (geneViz) {
          const updatedGeneViz = { ...geneViz, color };

          // Update in both visible genes (if present) and cache
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

          // Update in both visible genes (if present) and cache
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

    setShowAssigned: (show: boolean) =>
      set((state) => {
        const newSelectedGenes = new Map(state.selectedGenes);
        const newGeneDataCache = new Map(state.geneDataCache);

        for (const [key, geneViz] of newSelectedGenes) {
          newSelectedGenes.set(key, { ...geneViz, showAssigned: show });
        }
        for (const [key, geneViz] of newGeneDataCache) {
          newGeneDataCache.set(key, { ...geneViz, showAssigned: show });
        }

        return {
          showAssigned: show,
          selectedGenes: newSelectedGenes,
          geneDataCache: newGeneDataCache,
        };
      }),

    setShowUnassigned: (show: boolean) =>
      set((state) => {
        const newSelectedGenes = new Map(state.selectedGenes);
        const newGeneDataCache = new Map(state.geneDataCache);

        for (const [key, geneViz] of newSelectedGenes) {
          newSelectedGenes.set(key, { ...geneViz, showUnassigned: show });
        }
        for (const [key, geneViz] of newGeneDataCache) {
          newGeneDataCache.set(key, { ...geneViz, showUnassigned: show });
        }

        return {
          showUnassigned: show,
          selectedGenes: newSelectedGenes,
          geneDataCache: newGeneDataCache,
        };
      }),

    setGeneShowAssigned: (gene: string, show: boolean) =>
      set((state) => {
        const newSelectedGenes = new Map(state.selectedGenes);
        const newGeneDataCache = new Map(state.geneDataCache);
        const geneViz =
          newSelectedGenes.get(gene) || newGeneDataCache.get(gene);

        if (geneViz) {
          const updated = { ...geneViz, showAssigned: show };

          if (newSelectedGenes.has(gene)) newSelectedGenes.set(gene, updated);
          newGeneDataCache.set(gene, updated);
        }

        return { selectedGenes: newSelectedGenes, geneDataCache: newGeneDataCache };
      }),

    setGeneShowUnassigned: (gene: string, show: boolean) =>
      set((state) => {
        const newSelectedGenes = new Map(state.selectedGenes);
        const newGeneDataCache = new Map(state.geneDataCache);
        const geneViz =
          newSelectedGenes.get(gene) || newGeneDataCache.get(gene);

        if (geneViz) {
          const updated = { ...geneViz, showUnassigned: show };

          if (newSelectedGenes.has(gene)) newSelectedGenes.set(gene, updated);
          newGeneDataCache.set(gene, updated);
        }

        return { selectedGenes: newSelectedGenes, geneDataCache: newGeneDataCache };
      }),

    setGeneAssignedShape: (gene: string, shape: MoleculeShape) =>
      set((state) => {
        const newSelectedGenes = new Map(state.selectedGenes);
        const newGeneDataCache = new Map(state.geneDataCache);
        const geneViz =
          newSelectedGenes.get(gene) || newGeneDataCache.get(gene);

        if (geneViz) {
          const updated = { ...geneViz, assignedShape: shape };

          if (newSelectedGenes.has(gene)) newSelectedGenes.set(gene, updated);
          newGeneDataCache.set(gene, updated);
        }

        return { selectedGenes: newSelectedGenes, geneDataCache: newGeneDataCache };
      }),

    setGeneUnassignedShape: (gene: string, shape: MoleculeShape) =>
      set((state) => {
        const newSelectedGenes = new Map(state.selectedGenes);
        const newGeneDataCache = new Map(state.geneDataCache);
        const geneViz =
          newSelectedGenes.get(gene) || newGeneDataCache.get(gene);

        if (geneViz) {
          const updated = { ...geneViz, unassignedShape: shape };

          if (newSelectedGenes.has(gene)) newSelectedGenes.set(gene, updated);
          newGeneDataCache.set(gene, updated);
        }

        return { selectedGenes: newSelectedGenes, geneDataCache: newGeneDataCache };
      }),

    setGeneUnassignedColor: (gene: string, color: string) =>
      set((state) => {
        const newSelectedGenes = new Map(state.selectedGenes);
        const newGeneDataCache = new Map(state.geneDataCache);
        const geneViz =
          newSelectedGenes.get(gene) || newGeneDataCache.get(gene);

        if (geneViz) {
          const updated = { ...geneViz, unassignedColor: color };

          if (newSelectedGenes.has(gene)) newSelectedGenes.set(gene, updated);
          newGeneDataCache.set(gene, updated);
        }

        return { selectedGenes: newSelectedGenes, geneDataCache: newGeneDataCache };
      }),

    setGeneUnassignedLocalScale: (gene: string, scale: number) =>
      set((state) => {
        const newSelectedGenes = new Map(state.selectedGenes);
        const newGeneDataCache = new Map(state.geneDataCache);
        const geneViz =
          newSelectedGenes.get(gene) || newGeneDataCache.get(gene);

        if (geneViz) {
          const updated = { ...geneViz, unassignedLocalScale: scale };

          if (newSelectedGenes.has(gene)) newSelectedGenes.set(gene, updated);
          newGeneDataCache.set(gene, updated);
        }

        return { selectedGenes: newSelectedGenes, geneDataCache: newGeneDataCache };
      }),

    clearGenes: () =>
      set({
        selectedGenes: new Map(),
        selectedGenesLegend: new Set(),
        geneDataCache: new Map(),
        geneColorSlots: new Map(),
      }),

  }));
