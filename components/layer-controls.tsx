"use client";

import type { LayerContext } from "@/lib/scene/registry";

import { ELEMENTS } from "@/lib/scene/registry";
import { usePanelLayers } from "@/lib/hooks/usePanelStores";

/**
 * The controls of every layer the panel's scene is drawing, for one slot:
 * the buttons of the viewer's rail, or the floating panels. Mounted once
 * per slot by the viewer; a layer's controls appear with the layer.
 */
export function LayerControls({ slot }: { slot: "rail" | "panel" }) {
  const layers = usePanelLayers();

  return (
    <>
      {ELEMENTS.map((def) => {
        const layer = layers[def.kind];

        if (!layer || !def.Controls || (def.slot ?? "panel") !== slot) {
          return null;
        }

        return <def.Controls key={def.kind} layer={layer} />;
      })}
    </>
  );
}

/**
 * Keeps every layer in step with the panel's stores (a gene selection, a
 * visibility toggle). Rendered by the scene that owns the layers.
 */
export function LayerBindings({ ctx }: { ctx: LayerContext }) {
  const layers = usePanelLayers();

  return (
    <>
      {ELEMENTS.map((def) => {
        const layer = layers[def.kind];

        if (!layer || !def.Bind) return null;

        return <def.Bind key={def.kind} ctx={ctx} layer={layer} />;
      })}
    </>
  );
}
