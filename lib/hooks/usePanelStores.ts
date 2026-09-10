"use client";

import type { VisualizationState } from "../stores/createVisualizationStore";
import type { DatasetState } from "../stores/createDatasetStore";
import type { SingleMoleculeState } from "../stores/createSingleMoleculeStore";
import type { SingleMoleculeVisualizationState } from "../stores/createSingleMoleculeVisualizationStore";
import type {
  LabelledMoleculeVisualizationState,
  LabelledMoleculeVisualizationStore,
} from "../stores/createLabelledMoleculeVisualizationStore";

import { useContext } from "react";
import { useStore } from "zustand";

import { PanelContext } from "../contexts/PanelContext";
import { useVisualizationStore } from "../stores/visualizationStore";
import { useDatasetStore } from "../stores/datasetStore";
import { useSingleMoleculeStore } from "../stores/singleMoleculeStore";
import { useSingleMoleculeVisualizationStore } from "../stores/singleMoleculeVisualizationStore";
import {
  labelledMoleculeVisualizationStore,
  useLabelledMoleculeVisualizationStore,
} from "../stores/labelledMoleculeVisualizationStore";

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
export function usePanelDatasetStore<T>(selector: (s: DatasetState) => T): T;
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
    selector ?? ((s: SingleMoleculeVisualizationState) => s as unknown as T);

  if (ctx) {
    return useStore(ctx.singleMoleculeVisualizationStore, sel);
  }

  return useSingleMoleculeVisualizationStore(sel);
}

export function usePanelLabelledMoleculeVisualizationStore(): LabelledMoleculeVisualizationState;
export function usePanelLabelledMoleculeVisualizationStore<T>(
  selector: (s: LabelledMoleculeVisualizationState) => T,
): T;
export function usePanelLabelledMoleculeVisualizationStore<T>(
  selector?: (s: LabelledMoleculeVisualizationState) => T,
) {
  const ctx = useContext(PanelContext);
  const sel =
    selector ?? ((s: LabelledMoleculeVisualizationState) => s as unknown as T);

  if (ctx) {
    return useStore(ctx.labelledMoleculeVisualizationStore, sel);
  }

  return useLabelledMoleculeVisualizationStore(sel);
}

/**
 * The vanilla store behind the hook above, for getState/setState outside
 * React — event handlers and the WebGL frame loop, which must not re-render
 * to read a value. Same panel-or-global resolution.
 */
export function usePanelLabelledMoleculeApi(): LabelledMoleculeVisualizationStore {
  const ctx = useContext(PanelContext);

  return (
    ctx?.labelledMoleculeVisualizationStore ??
    labelledMoleculeVisualizationStore
  );
}

export function usePanelId(): string | null {
  const ctx = useContext(PanelContext);

  return ctx?.panelId ?? null;
}
