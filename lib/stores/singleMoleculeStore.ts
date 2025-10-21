import { create } from "zustand";
import type { SingleMoleculeDataset } from "../SingleMoleculeDataset";

interface SingleMoleculeState {
  datasets: Map<string, SingleMoleculeDataset>;
  currentDatasetId: string | null;
  isLoading: boolean;
  error: string | null;

  // Actions
  addDataset: (dataset: SingleMoleculeDataset) => void;
  removeDataset: (id: string) => void;
  setCurrentDataset: (id: string | null) => void;
  getCurrentDataset: () => SingleMoleculeDataset | null;
  setLoading: (loading: boolean) => void;
  setError: (error: string | null) => void;
  clearError: () => void;
  reset: () => void;
}

export const useSingleMoleculeStore = create<SingleMoleculeState>(
  (set, get) => ({
    datasets: new Map(),
    currentDatasetId: null,
    isLoading: false,
    error: null,

    addDataset: (dataset: SingleMoleculeDataset) => {
      set((state) => {
        const newDatasets = new Map(state.datasets);
        newDatasets.set(dataset.id, dataset);

        console.log("Single molecule dataset added to store:", {
          id: dataset.id,
          name: dataset.name,
          type: dataset.type,
          summary: dataset.getSummary(),
        });

        console.log(dataset);

        return {
          datasets: newDatasets,
          currentDatasetId: dataset.id,
          error: null,
        };
      });
    },

    removeDataset: (id: string) => {
      set((state) => {
        const newDatasets = new Map(state.datasets);
        newDatasets.delete(id);

        return {
          datasets: newDatasets,
          currentDatasetId:
            state.currentDatasetId === id ? null : state.currentDatasetId,
        };
      });
    },

    setCurrentDataset: (id: string | null) => {
      set({ currentDatasetId: id });
    },

    getCurrentDataset: () => {
      const state = get();
      if (!state.currentDatasetId) return null;
      return state.datasets.get(state.currentDatasetId) || null;
    },

    setLoading: (loading: boolean) => {
      set({ isLoading: loading });
    },

    setError: (error: string | null) => {
      set({ error });
    },

    clearError: () => {
      set({ error: null });
    },

    reset: () => {
      set({
        datasets: new Map(),
        currentDatasetId: null,
        isLoading: false,
        error: null,
      });
    },
  })
);
