"use client";

import React, { useEffect, useRef, useState } from "react";

import { usePanelId } from "@/lib/hooks/usePanelStores";
import {
  encodeCellVizState,
  decodeCellVizState,
  encodeSMVizState,
  decodeSMVizState,
  type CellVizUrlState,
  type SMVizUrlState,
} from "@/lib/utils/url-visualization-state";
import {
  scheduleUrlUpdate,
  readUrlVizState,
} from "@/lib/utils/url-state-writer";
import { selectBestClusterColumn } from "@/lib/utils/dataset-utils";
import type { StandardizedDataset } from "@/lib/StandardizedDataset";
import type { SingleMoleculeDataset } from "@/lib/SingleMoleculeDataset";
import type { VisualizationState } from "@/lib/stores/createVisualizationStore";
import type { SingleMoleculeVisualizationState } from "@/lib/stores/createSingleMoleculeVisualizationStore";

// --- Panel helpers ---

function getPanel(panelId: string | null): "left" | "right" {
  return panelId === "right" ? "right" : "left";
}

// --- Apply decoded state to stores (exported for imperative use) ---

export function applyCellVizState(
  decoded: CellVizUrlState,
  store: VisualizationState,
  dataset: StandardizedDataset,
) {
  const columnNames = dataset.clusters?.map((c) => c.column) ?? [];
  let validColumn: string | null = null;

  if (decoded.c && columnNames.includes(decoded.c)) {
    validColumn = decoded.c;
  }

  if (!validColumn) {
    validColumn = selectBestClusterColumn(dataset);
  }
  store.setSelectedColumn(validColumn);

  if (decoded.e) {
    const embeddingNames = Object.keys(dataset.embeddings ?? {});

    if (embeddingNames.includes(decoded.e)) {
      store.setSelectedEmbedding(decoded.e);
    }
  }

  if (decoded.g && dataset.genes?.includes(decoded.g)) {
    store.setSelectedGene(decoded.g);
  }

  if (decoded.ct && decoded.ct.length > 0 && validColumn) {
    const cluster = dataset.clusters?.find((c) => c.column === validColumn);

    if (cluster) {
      const validValues = new Set(cluster.values.map(String));

      decoded.ct.forEach((ct) => {
        if (validValues.has(ct)) {
          store.toggleCelltype(ct);
        }
      });
    }
  }

  if (decoded.m && Array.isArray(decoded.m)) {
    const validModes = decoded.m.filter(
      (m) => m === "gene" || m === "celltype",
    );

    if (validModes.includes("gene") && !decoded.g) {
      const filtered = validModes.filter((m) => m !== "gene");

      store.setMode(filtered.length > 0 ? filtered : ["celltype"]);
    } else if (validModes.length > 0) {
      store.setMode(validModes);
    }
  }

  if (decoded.gs) {
    store.setGeneScaleMin(decoded.gs[0]);
    store.setGeneScaleMax(decoded.gs[1]);
  }
  if (decoded.ns) {
    store.setNumericalScaleMin(decoded.ns[0]);
    store.setNumericalScaleMax(decoded.ns[1]);
  }

  if (decoded.sz !== undefined) {
    store.setSizeScale(decoded.sz);
  }
}

export function applySMVizState(
  decoded: SMVizUrlState,
  store: SingleMoleculeVisualizationState,
  dataset: SingleMoleculeDataset,
) {
  const validGenes = new Set(dataset.uniqueGenes);

  store.clearGenes();

  if (decoded.genes && decoded.genes.length > 0) {
    decoded.genes.forEach(([name, color, localScale, isVisible]) => {
      if (!validGenes.has(name)) return;
      store.addGene(name, color, localScale);

      if (!isVisible) {
        store.toggleGeneVisibility(name);
      }
    });
  }

  if (decoded.gs !== undefined) {
    store.setGlobalScale(decoded.gs);
  }

  if (decoded.vm) {
    store.setViewMode(decoded.vm);
  }
}

/**
 * Try to read and decode cell viz state from URL for a given panel.
 * Returns the decoded state if valid, null otherwise.
 * This is a synchronous, non-hook function for use in imperative code (load callbacks).
 */
export function tryReadCellVizFromUrl(
  panel: "left" | "right",
): CellVizUrlState | null {
  const urlState = readUrlVizState();
  const encoded = panel === "left" ? urlState.left : urlState.right;

  if (!encoded) return null;

  return decodeCellVizState(encoded);
}

/**
 * Try to read and decode SM viz state from URL for a given panel.
 */
export function tryReadSMVizFromUrl(
  panel: "left" | "right",
): SMVizUrlState | null {
  const urlState = readUrlVizState();
  const encoded = panel === "left" ? urlState.left : urlState.right;

  if (!encoded) return null;

  return decodeSMVizState(encoded);
}

// --- Hooks ---

/**
 * Syncs single cell visualization state to/from URL params.
 * Call `tryApplyUrlState(dataset)` once after dataset loads to apply URL state.
 * Returns { hasUrlState, tryApplyUrlState }.
 */
export function useCellVizUrlSync(
  datasetReady: boolean,
  dataset: StandardizedDataset | null,
  store: VisualizationState,
): { hasUrlState: boolean; hasUrlStateRef: React.RefObject<boolean> } {
  const panelId = usePanelId();
  const panel = getPanel(panelId);
  const [hasUrlState, setHasUrlState] = useState(false);
  const hasUrlStateRef = useRef(false);
  const appliedRef = useRef(false);

  // Reading: apply URL state once after dataset is ready
  useEffect(() => {
    if (!datasetReady || !dataset || appliedRef.current) return;
    appliedRef.current = true;

    const decoded = tryReadCellVizFromUrl(panel);

    if (!decoded) return;

    applyCellVizState(decoded, store, dataset);
    hasUrlStateRef.current = true;
    setHasUrlState(true);
  }, [datasetReady, dataset, store, panel]);

  // Writing: encode state changes to URL (debounced)
  const {
    selectedGene,
    selectedColumn,
    selectedCelltypes,
    mode,
    geneScaleMin,
    geneScaleMax,
    numericalScaleMin,
    numericalScaleMax,
    sizeScale,
    selectedEmbedding,
  } = store;

  useEffect(() => {
    if (!datasetReady) return;

    const encoded = encodeCellVizState({
      selectedGene,
      selectedColumn,
      selectedCelltypes,
      mode,
      geneScaleMin,
      geneScaleMax,
      numericalScaleMin,
      numericalScaleMax,
      sizeScale,
      selectedEmbedding,
    });

    scheduleUrlUpdate(panel, encoded);
  }, [
    datasetReady,
    panel,
    selectedGene,
    selectedColumn,
    selectedCelltypes,
    mode,
    geneScaleMin,
    geneScaleMax,
    numericalScaleMin,
    numericalScaleMax,
    sizeScale,
    selectedEmbedding,
  ]);

  return { hasUrlState, hasUrlStateRef };
}

/**
 * Syncs single molecule visualization state to/from URL params.
 */
export function useSMVizUrlSync(
  datasetReady: boolean,
  dataset: SingleMoleculeDataset | null,
  store: SingleMoleculeVisualizationState,
): { hasUrlState: boolean; hasUrlStateRef: React.RefObject<boolean> } {
  const panelId = usePanelId();
  const panel = getPanel(panelId);
  const [hasUrlState, setHasUrlState] = useState(false);
  const hasUrlStateRef = useRef(false);
  const appliedRef = useRef(false);

  // Reading: apply URL state once after dataset is ready
  useEffect(() => {
    if (!datasetReady || !dataset || appliedRef.current) return;
    appliedRef.current = true;

    const decoded = tryReadSMVizFromUrl(panel);

    if (!decoded) return;

    applySMVizState(decoded, store, dataset);
    hasUrlStateRef.current = true;
    setHasUrlState(true);
  }, [datasetReady, dataset, store, panel]);

  // Writing: encode state changes to URL (debounced)
  const { selectedGenes, geneDataCache, globalScale, viewMode } = store;

  useEffect(() => {
    if (!datasetReady) return;

    const encoded = encodeSMVizState({
      selectedGenes,
      geneDataCache,
      globalScale,
      viewMode,
    });

    scheduleUrlUpdate(panel, encoded);
  }, [
    datasetReady,
    panel,
    selectedGenes,
    geneDataCache,
    globalScale,
    viewMode,
  ]);

  return { hasUrlState, hasUrlStateRef };
}
