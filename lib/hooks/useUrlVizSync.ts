"use client";

import type { StandardizedDataset } from "@/lib/StandardizedDataset";
import type { SingleMoleculeDataset } from "@/lib/SingleMoleculeDataset";
import type { VisualizationState } from "@/lib/stores/createVisualizationStore";
import type { SingleMoleculeVisualizationState } from "@/lib/stores/createSingleMoleculeVisualizationStore";

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
import { pickDefaultGenes } from "@/lib/utils/auto-select-genes";

// --- Panel helpers ---

function getPanel(panelId: string | null): "left" | "right" {
  return panelId === "right" ? "right" : "left";
}

// --- Apply decoded state to stores (exported for imperative use) ---

export async function applyCellVizState(
  decoded: CellVizUrlState,
  store: VisualizationState,
  dataset: StandardizedDataset,
) {
  // Check loaded clusters first, then fall back to all known column names
  const loadedColumnNames = dataset.clusters?.map((c) => c.column) ?? [];
  const allColumnNames =
    dataset.allClusterColumnNames && dataset.allClusterColumnNames.length > 0
      ? dataset.allClusterColumnNames
      : loadedColumnNames;
  let validColumn: string | null = null;

  if (decoded.c && allColumnNames.includes(decoded.c)) {
    validColumn = decoded.c;
  }

  if (!validColumn) {
    validColumn = selectBestClusterColumn(dataset);
  }

  // Apply type overrides before resolving isNumerical
  if (decoded.to && typeof decoded.to === "object") {
    store.setColumnTypeOverrides(decoded.to);
  }

  const overrides = decoded.to ?? {};
  const isNumerical = overrides[validColumn ?? ""]
    ? overrides[validColumn ?? ""] === "numerical"
    : dataset.allClusterColumnTypes?.[validColumn ?? ""] === "numerical";

  store.setSelectedColumn(validColumn, isNumerical);

  // If the selected column is known but not loaded, fetch it on demand
  if (
    validColumn &&
    !dataset.clusters?.some((c) => c.column === validColumn) &&
    dataset.adapter
  ) {
    try {
      let newClusters: Array<{
        column: string;
        type: string;
        values: any[];
        palette: Record<string, string> | null;
      }> | null = null;

      if (dataset.adapter.mode === "local") {
        newClusters = await dataset.adapter.loadClusters([validColumn]);
      } else {
        const { getStandardizedDatasetWorker } = await import(
          "@/lib/workers/standardizedDatasetWorkerManager"
        );
        const worker = await getStandardizedDatasetWorker();

        newClusters = await worker.loadClusterFromS3(
          dataset.id,
          [validColumn],
          dataset.metadata?.customS3BaseUrl,
        );
      }

      if (newClusters && newClusters.length > 0) {
        dataset.addClusters(newClusters);
        store.incrementClusterVersion();
      }
    } catch (error) {
      console.warn(
        `Failed to load cluster column "${validColumn}" from URL:`,
        error,
      );
    }
  }

  if (decoded.e) {
    const embeddingNames =
      dataset.allEmbeddingNames && dataset.allEmbeddingNames.length > 0
        ? dataset.allEmbeddingNames
        : Object.keys(dataset.embeddings ?? {});

    if (embeddingNames.includes(decoded.e)) {
      store.setSelectedEmbedding(decoded.e);
    }
  }

  if (decoded.g && dataset.genes?.includes(decoded.g)) {
    store.setSelectedGene(decoded.g);
  }

  // Apply celltypes (cluster data should now be loaded from the on-demand fetch above)
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

  // If URL had no genes (or none were valid), fall back to auto-select defaults
  if (store.selectedGenes.size === 0) {
    pickDefaultGenes(dataset.uniqueGenes).forEach((gene) =>
      store.addGene(gene),
    );
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

    // Set flag synchronously so auto-select effects skip
    hasUrlStateRef.current = true;
    setHasUrlState(true);

    // Apply state (may trigger async cluster loading for unloaded columns)
    applyCellVizState(decoded, store, dataset);
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
    columnTypeOverrides,
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
      columnTypeOverrides,
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
    columnTypeOverrides,
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
