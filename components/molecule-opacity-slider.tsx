"use client";

import { Slider } from "@heroui/react";
import { Tooltip } from "@heroui/tooltip";

import { glassButton } from "@/components/primitives";
import { usePanelSingleMoleculeVisualizationStore } from "@/lib/hooks/usePanelStores";

/**
 * The molecules' opacity, for this panel: every gene, assigned and
 * unassigned, in the per-gene layer and in a store's tiled points. Below 100%
 * molecules blend, so where they pile up reads as denser — which is how a
 * spot's localization is judged once zoomed out. Sits beside the dot-size
 * slider and looks like it.
 */
export function MoleculeOpacitySlider() {
  const opacity = usePanelSingleMoleculeVisualizationStore(
    (s) => s.moleculeOpacity,
  );
  const setOpacity = usePanelSingleMoleculeVisualizationStore(
    (s) => s.setMoleculeOpacity,
  );

  return (
    <Tooltip
      content={`Molecule opacity: ${Math.round(opacity * 100)}%`}
      placement="right"
    >
      <div
        className={`w-14 h-32 rounded-full border-2 border-default-200 p-2 flex flex-col items-center justify-center ${glassButton()}`}
        data-testid="molecule-opacity"
      >
        <Slider
          aria-label="Molecule opacity"
          className="h-full"
          maxValue={1}
          minValue={0.05}
          orientation="vertical"
          size="sm"
          step={0.05}
          value={opacity}
          onChange={(value) => setOpacity(value as number)}
        />
      </div>
    </Tooltip>
  );
}
