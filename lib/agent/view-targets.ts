import type { StandardizedDataset } from "@/lib/StandardizedDataset";
import type { SingleMoleculeDataset } from "@/lib/SingleMoleculeDataset";
import type { VisualizationState } from "@/lib/stores/createVisualizationStore";
import type { SingleMoleculeVisualizationState } from "@/lib/stores/createSingleMoleculeVisualizationStore";

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

/** A molecules panel: the molecule viewer's own (`molecules:<panel>`), or
 * molecules drawn over cells (`molecules-overlay:<panel>`). */
export interface MoleculesTarget {
  kind: "molecules";
  panel: "left" | "right";
  overlay: boolean;
  store: () => SingleMoleculeVisualizationState;
  dataset: () => SingleMoleculeDataset | null;
}

export type ViewTarget = CellsTarget | MoleculesTarget;

const targets = new Map<string, ViewTarget>();

export function registerViewTarget(key: string, target: ViewTarget) {
  targets.set(key, target);

  return () => {
    if (targets.get(key) === target) targets.delete(key);
  };
}

/** The cells view a command acts on: the left panel's, else any. */
export function primaryTarget(): CellsTarget | null {
  const left = targets.get("cells:left");

  if (left?.kind === "cells") return left;
  for (const t of targets.values()) if (t.kind === "cells") return t;

  return null;
}

/** The molecules a command acts on: the molecule viewer's, else those
 * drawn over the cells (left panel first). */
export function moleculesTarget(): MoleculesTarget | null {
  for (const key of ["molecules:left", "molecules-overlay:left", "molecules:right", "molecules-overlay:right"]) {
    const t = targets.get(key);

    if (t?.kind === "molecules") return t;
  }

  return null;
}

/** What an agent needs of a cells panel's scene to frame it: its camera,
 * controls and point positions (registered by components/three-scene.tsx). */
export interface SceneHandle {
  camera: () => {
    position: { x: number; y: number; z: number; copy: (v: any) => any; clone: () => any };
    fov?: number;
    aspect?: number;
  } | null;
  controls: () => {
    target: { x: number; y: number; z: number; copy: (v: any) => any; clone: () => any; set: (x: number, y: number, z: number) => any };
    update: () => void;
  } | null;
  /** x, y, z per cell, in the scene's coordinates, in cell order. */
  positions: () => Float32Array | null;
}

const scenes = new Map<string, SceneHandle>();

export function registerScene(key: string, handle: SceneHandle) {
  scenes.set(key, handle);

  return () => {
    if (scenes.get(key) === handle) scenes.delete(key);
  };
}

export function sceneFor(panel: "left" | "right"): SceneHandle | null {
  return scenes.get(`cells:${panel}`) ?? null;
}
