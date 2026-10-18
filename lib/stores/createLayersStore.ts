import type { ElementKind, Layer } from "@/lib/scene/contract";

import { createStore } from "zustand";

/** The layers a panel's scene is drawing right now, by element kind. */
export type Layers = Partial<Record<ElementKind, Layer>>;

export interface LayersState {
  layers: Layers;
  setLayer: (kind: ElementKind, layer: Layer | null) => void;
}

/**
 * One per panel (the split screen's right panel has its own, the left one
 * uses the global instance). Written by `useLayers` in the scene that owns
 * the layers; read by the controls rails and panels, which sit beside the
 * scene rather than inside it.
 */
export function createLayersStoreInstance() {
  return createStore<LayersState>((set) => ({
    layers: {},
    setLayer: (kind, layer) =>
      set((state) => {
        if ((state.layers[kind] ?? null) === layer) return state;
        const layers = { ...state.layers };

        if (layer) layers[kind] = layer;
        else delete layers[kind];

        return { layers };
      }),
  }));
}

export type LayersStore = ReturnType<typeof createLayersStoreInstance>;
