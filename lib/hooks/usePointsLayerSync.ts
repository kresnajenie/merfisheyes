"use client";

import type { GeneStyle, PointsLayer } from "@/lib/spatialdata/PointsLayer";

import { useEffect } from "react";

import { usePanelSingleMoleculeVisualizationStore } from "./usePanelStores";

import { VISUALIZATION_CONFIG } from "@/lib/config/visualization.config";
import { SM_DOT_SIZE_FACTOR } from "@/lib/webgl/point-cloud";

/**
 * Render order of a gene's molecules: later in the legend → drawn later
 * (on top). Everything stays below the cells' point cloud (render order 0).
 */
export function moleculeRenderOrder(legendIndex: number): number {
  return -1000 + legendIndex * 2 + 1;
}

/**
 * Keep a `PointsLayer` in step with this panel's single-molecule
 * visualization store: which genes to draw, in which colour, shape and size,
 * and whether the layer is shown at all.
 */
export function usePointsLayerSync(layer: PointsLayer | null): void {
  const selectedGenes = usePanelSingleMoleculeVisualizationStore(
    (s) => s.selectedGenes,
  );
  const legend = usePanelSingleMoleculeVisualizationStore(
    (s) => s.selectedGenesLegend,
  );
  const geneDataCache = usePanelSingleMoleculeVisualizationStore(
    (s) => s.geneDataCache,
  );
  const globalScale = usePanelSingleMoleculeVisualizationStore(
    (s) => s.globalScale,
  );
  const opacity = usePanelSingleMoleculeVisualizationStore(
    (s) => s.moleculeOpacity,
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
    const genes = new Map<string, GeneStyle>();

    // Every gene in the legend, hidden ones (the eye, or another gene's
    // cmd+click solo) included: they keep their share of the budget and what
    // they have loaded, so showing them again is instant
    let i = 0;
    const size = (localScale: number) =>
      VISUALIZATION_CONFIG.SINGLE_MOLECULE_POINT_BASE_SIZE *
      globalScale *
      localScale *
      SM_DOT_SIZE_FACTOR;

    for (const gene of legend) {
      const index = i++;
      const viz = selectedGenes.get(gene) ?? geneDataCache.get(gene);

      if (!viz) continue;
      // The same rules as MoleculeLayer.sync: each kind by the panel's
      // switch and the gene's own
      const assigned = showAssigned && viz.showAssigned;
      const unassigned =
        layer.marksUnassigned && showUnassigned && viz.showUnassigned;

      genes.set(gene, {
        visible: selectedGenes.has(gene) && (assigned || unassigned),
        assigned,
        unassigned: unassigned
          ? {
              color: viz.unassignedColor,
              shape: viz.unassignedShape,
              dotSize: size(viz.unassignedLocalScale),
            }
          : null,
        color: viz.color,
        shape: viz.assignedShape,
        dotSize: size(viz.localScale),
        opacity,
        renderOrder: moleculeRenderOrder(index),
      });
    }
    layer.visible = visible;
    layer.setGenes(genes);
  }, [
    layer,
    selectedGenes,
    legend,
    geneDataCache,
    globalScale,
    opacity,
    visible,
    showAssigned,
    showUnassigned,
  ]);
}
