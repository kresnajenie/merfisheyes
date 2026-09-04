import { useStore } from "zustand";

import { createLabelledMoleculeVisualizationStoreInstance } from "./createLabelledMoleculeVisualizationStore";

import type { LabelledMoleculeVisualizationState } from "./createLabelledMoleculeVisualizationStore";

/**
 * Global singleton for the labelled single-molecule viewer, mirroring the
 * single-cell / single-molecule stores. The factory exists separately so a
 * split panel can own an isolated instance later.
 */
const store = createLabelledMoleculeVisualizationStoreInstance();

export function useLabelledMoleculeVisualizationStore(): LabelledMoleculeVisualizationState;
export function useLabelledMoleculeVisualizationStore<T>(
  selector: (state: LabelledMoleculeVisualizationState) => T,
): T;
export function useLabelledMoleculeVisualizationStore<T>(
  selector?: (state: LabelledMoleculeVisualizationState) => T,
) {
  return useStore(store, selector as (s: LabelledMoleculeVisualizationState) => T);
}

export const labelledMoleculeVisualizationStore = store;
