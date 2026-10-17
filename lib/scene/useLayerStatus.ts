"use client";

import type { Layer } from "./contract";

import { useEffect, useReducer } from "react";

/**
 * Re-render a controls component whenever its layer reports a change. The
 * layer pushes (contract: `subscribe`); nothing is polled.
 */
export function useLayerStatus(layer: Layer | null): number {
  const [version, bump] = useReducer((n: number) => n + 1, 0);

  useEffect(() => {
    if (!layer) return;
    bump();

    return layer.subscribe(bump);
  }, [layer]);

  return version;
}
