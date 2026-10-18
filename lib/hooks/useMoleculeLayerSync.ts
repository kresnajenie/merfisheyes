"use client";

import type { MoleculeLayer } from "@/lib/webgl/MoleculeLayer";

import { useEffect } from "react";

import { usePanelSingleMoleculeVisualizationStore } from "./usePanelStores";

/**
 * Keep a `MoleculeLayer` in step with this panel's single-molecule
 * visualization store. Both the standalone molecule scene and the cell
 * scene's overlay wire their layer through this, so what is drawn for a
 * selection is decided in one place.
 *
 * `onSynced` runs after each reconciliation settles (whole loads finished,
 * streams started) — for the test hooks and a fallback camera fit.
 */
export function useMoleculeLayerSync(
  layer: MoleculeLayer | null,
  onSynced?: (layer: MoleculeLayer) => void,
): void {
  const selectedGenes = usePanelSingleMoleculeVisualizationStore(
    (s) => s.selectedGenes,
  );
  const legendOrder = usePanelSingleMoleculeVisualizationStore(
    (s) => s.selectedGenesLegend,
  );
  const globalScale = usePanelSingleMoleculeVisualizationStore(
    (s) => s.globalScale,
  );
  const showAssigned = usePanelSingleMoleculeVisualizationStore(
    (s) => s.showAssigned,
  );
  const showUnassigned = usePanelSingleMoleculeVisualizationStore(
    (s) => s.showUnassigned,
  );
  const visible = usePanelSingleMoleculeVisualizationStore(
    (s) => s.smLayerVisible,
  );

  useEffect(() => {
    if (!layer) return;
    let stale = false;

    layer
      .sync({
        selectedGenes,
        legendOrder,
        globalScale,
        showAssigned,
        showUnassigned,
        visible,
      })
      .then(() => {
        if (!stale) onSynced?.(layer);
      });

    return () => {
      stale = true;
    };
    // onSynced is a fresh closure each render; the sync is what matters.
  }, [
    layer,
    selectedGenes,
    legendOrder,
    globalScale,
    showAssigned,
    showUnassigned,
    visible,
  ]);
}
