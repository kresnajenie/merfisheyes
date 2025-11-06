import type { StandardizedDataset } from "../StandardizedDataset";
import type { SingleMoleculeDataset } from "../SingleMoleculeDataset";

import { create } from "zustand";

type Dataset = StandardizedDataset | SingleMoleculeDataset;

interface DatasetState {
  datasets: Map<string, Dataset>;
  currentDatasetId: string | null;
  isLoading: boolean;
  error: string | null;

  // Actions
  addDataset: (dataset: Dataset) => void;
  removeDataset: (id: string) => void;
  setCurrentDataset: (id: string | null) => void;
  getCurrentDataset: () => Dataset | null;
  setLoading: (loading: boolean) => void;
  setError: (error: string | null) => void;
  clearError: () => void;
  reset: () => void;
}

export const useDatasetStore = create<DatasetState>((set, get) => ({
  datasets: new Map(),
  currentDatasetId: null,
  isLoading: false,
  error: null,

  addDataset: (dataset: Dataset) => {
    set((state) => {
      const newDatasets = new Map(state.datasets);

      newDatasets.set(dataset.id, dataset);

      if (typeof window !== "undefined") {
        try {
          window.localStorage.setItem("lastDatasetMode", "cell");
        } catch (error) {
          console.warn(
            "[DatasetStore] Failed to persist last dataset mode:",
            error,
          );
        }
      }

      console.log("Dataset added to store:", {
        id: dataset.id,
        name: dataset.name,
        type: dataset.type,
        summary: dataset.getSummary(),
      });

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
}));
