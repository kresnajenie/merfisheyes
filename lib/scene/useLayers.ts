"use client";

import type { ElementKind, Layer, SceneHandle } from "./contract";
import type { LayerContext, OpenedElements } from "./registry";
import type { Layers } from "@/lib/stores/createLayersStore";

import { useCallback, useEffect, useRef } from "react";
import { toast } from "react-toastify";
import { useStore } from "zustand";

import { ELEMENTS } from "./registry";

import { usePanelLayersApi } from "@/lib/hooks/usePanelStores";

/**
 * The layers for a set of opened elements, in the scene that calls this.
 *
 * A layer is created when its element appears and disposed when the element
 * changes or the scene component unmounts; in between it outlives scene
 * rebuilds (a view-mode switch): `onSceneReady` attaches every layer to the
 * new scene and its cleanup detaches them, so what a layer loaded and how
 * it was set stay. The layers are published to the panel's layers store,
 * where `LayerControls` finds them.
 */
export function useLayers(
  elements: OpenedElements,
  ctx: LayerContext = {},
): {
  /** Give this to the scene; it attaches the layers and returns the detach. */
  onSceneReady: (scene: SceneHandle) => () => void;
  layers: Layers;
} {
  const api = usePanelLayersApi();
  const layers = useStore(api, (s) => s.layers);
  const ctxRef = useRef(ctx);

  ctxRef.current = ctx;
  const sceneRef = useRef<SceneHandle | null>(null);
  const instances = useRef(
    new Map<ElementKind, { opened: unknown; layer: Layer }>(),
  );
  const detachers = useRef(new Map<ElementKind, () => void>());

  const attach = (kind: ElementKind, layer: Layer) => {
    const scene = sceneRef.current;
    const def = ELEMENTS.find((d) => d.kind === kind)!;

    if (!scene || detachers.current.has(kind)) return;
    if (def.canAttach && !def.canAttach(scene)) return;
    detachers.current.set(kind, scene.attach(layer));
  };
  const detach = (kind: ElementKind) => {
    detachers.current.get(kind)?.();
    detachers.current.delete(kind);
  };
  const drop = (kind: ElementKind) => {
    const current = instances.current.get(kind);

    if (!current) return;
    detach(kind);
    current.layer.dispose();
    instances.current.delete(kind);
    api.getState().setLayer(kind, null);
  };

  // One layer per element, for as long as that element is the one given
  useEffect(() => {
    for (const def of ELEMENTS) {
      const opened = elements[def.kind];
      const current = instances.current.get(def.kind);

      if (current?.opened === opened) continue;
      drop(def.kind);
      if (!opened) continue;
      const layer = def.createLayer(opened as never, ctxRef.current);

      instances.current.set(def.kind, { opened, layer });
      attach(def.kind, layer);
      api.getState().setLayer(def.kind, layer);
    }
  }, [elements.image, elements.shapes, elements.points, elements.molecules]);

  useEffect(
    () => () => {
      for (const kind of [...instances.current.keys()]) drop(kind);
    },
    [],
  );

  // What the store holds but cannot be drawn is said, not just left out
  useEffect(() => {
    for (const problem of elements.problems ?? []) {
      toast.warning(`Not drawn — ${problem}`, {
        toastId: `element-problem:${problem}`,
        autoClose: 12_000,
      });
    }
  }, [elements.problems]);

  const onSceneReady = useCallback((scene: SceneHandle) => {
    sceneRef.current = scene;
    for (const [kind, { layer }] of instances.current) attach(kind, layer);

    return () => {
      for (const kind of [...detachers.current.keys()]) detach(kind);
      if (sceneRef.current === scene) sceneRef.current = null;
    };
  }, []);

  return { onSceneReady, layers };
}
