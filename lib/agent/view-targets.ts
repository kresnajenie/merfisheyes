import type { StandardizedDataset } from "@/lib/StandardizedDataset";
import type { VisualizationState } from "@/lib/stores/createVisualizationStore";

/**
 * The views an agent can drive (docs/plans/MCP.md pillar D3): each open
 * cells panel registers here — from `useCellVizUrlSync`, the one hook every
 * cells viewer calls with its dataset and store — and unregisters when it
 * closes. Getters, not values: the store's state and the dataset change.
 */
export interface CellsTarget {
  kind: "cells";
  panel: "left" | "right";
  store: () => VisualizationState;
  dataset: () => StandardizedDataset | null;
}

export type ViewTarget = CellsTarget;

const targets = new Map<string, ViewTarget>();

export function registerViewTarget(key: string, target: ViewTarget) {
  targets.set(key, target);

  return () => {
    if (targets.get(key) === target) targets.delete(key);
  };
}

/** The view a command acts on: the left panel's, else any. */
export function primaryTarget(): ViewTarget | null {
  return targets.get("cells:left") ?? targets.values().next().value ?? null;
}
