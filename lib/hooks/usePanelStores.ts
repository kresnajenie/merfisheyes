"use client";

import { useContext } from "react";
import { useStore } from "zustand";

import { PanelContext } from "../contexts/PanelContext";
import { useVisualizationStore } from "../stores/visualizationStore";
import { useDatasetStore } from "../stores/datasetStore";
import { useSingleMoleculeStore } from "../stores/singleMoleculeStore";
import { useSingleMoleculeVisualizationStore } from "../stores/singleMoleculeVisualizationStore";
import type { VisualizationState } from "../stores/createVisualizationStore";
import type { DatasetState } from "../stores/createDatasetStore";
import type { SingleMoleculeState } from "../stores/createSingleMoleculeStore";
import type { SingleMoleculeVisualizationState } from "../stores/createSingleMoleculeVisualizationStore";

// Overloads: with selector returns T, without selector returns full state
export function usePanelVisualizationStore(): VisualizationState;
export function usePanelVisualizationStore<T>(
  selector: (s: VisualizationState) => T,
): T;
export function usePanelVisualizationStore<T>(
  selector?: (s: VisualizationState) => T,
) {
  const ctx = useContext(PanelContext);
  const sel = selector ?? ((s: VisualizationState) => s as unknown as T);

  if (ctx) {
    return useStore(ctx.visualizationStore, sel);
  }

  return useVisualizationStore(sel);
}

export function usePanelDatasetStore(): DatasetState;
export function usePanelDatasetStore<T>(
  selector: (s: DatasetState) => T,
): T;
export function usePanelDatasetStore<T>(selector?: (s: DatasetState) => T) {
  const ctx = useContext(PanelContext);
  const sel = selector ?? ((s: DatasetState) => s as unknown as T);

  if (ctx) {
    return useStore(ctx.datasetStore, sel);
  }

  return useDatasetStore(sel);
}

export function usePanelSingleMoleculeStore(): SingleMoleculeState;
export function usePanelSingleMoleculeStore<T>(
  selector: (s: SingleMoleculeState) => T,
): T;
export function usePanelSingleMoleculeStore<T>(
  selector?: (s: SingleMoleculeState) => T,
) {
  const ctx = useContext(PanelContext);
  const sel = selector ?? ((s: SingleMoleculeState) => s as unknown as T);

  if (ctx) {
    return useStore(ctx.singleMoleculeStore, sel);
  }

  return useSingleMoleculeStore(sel);
}

export function usePanelSingleMoleculeVisualizationStore(): SingleMoleculeVisualizationState;
export function usePanelSingleMoleculeVisualizationStore<T>(
  selector: (s: SingleMoleculeVisualizationState) => T,
): T;
export function usePanelSingleMoleculeVisualizationStore<T>(
  selector?: (s: SingleMoleculeVisualizationState) => T,
) {
  const ctx = useContext(PanelContext);
  const sel =
    selector ??
    ((s: SingleMoleculeVisualizationState) => s as unknown as T);

  if (ctx) {
    return useStore(ctx.singleMoleculeVisualizationStore, sel);
  }

  return useSingleMoleculeVisualizationStore(sel);
}

export function usePanelId(): string | null {
  const ctx = useContext(PanelContext);

  return ctx?.panelId ?? null;
}
