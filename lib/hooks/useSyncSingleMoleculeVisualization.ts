"use client";

import type { SingleMoleculeVisualizationStore } from "../stores/createSingleMoleculeVisualizationStore";
import type { SingleMoleculeVisualizationState } from "../stores/createSingleMoleculeVisualizationStore";
import type { SingleMoleculeStore } from "../stores/createSingleMoleculeStore";
import type { SingleMoleculeDataset } from "../SingleMoleculeDataset";

import { toast } from "react-toastify";
import { useEffect, useRef } from "react";

import { useSingleMoleculeVisualizationStore } from "../stores/singleMoleculeVisualizationStore";
import { useSingleMoleculeStore } from "../stores/singleMoleculeStore";
import { useSplitScreenStore } from "../stores/splitScreenStore";

type SMSyncFields = {
  selectedGenes: Map<
    string,
    { gene: string; color: string; localScale: number }
  >;
  selectedGenesLegend: Set<string>;
  geneDataCache: Map<
    string,
    { gene: string; color: string; localScale: number }
  >;
};

function pickSMSyncFields(
  state: SingleMoleculeVisualizationState,
): SMSyncFields {
  return {
    selectedGenes: state.selectedGenes,
    selectedGenesLegend: state.selectedGenesLegend,
    geneDataCache: state.geneDataCache,
  };
}

function getSMDataset(
  store: SingleMoleculeStore | typeof useSingleMoleculeStore,
): SingleMoleculeDataset | null {
  const state = store.getState();
  const dataset = state.currentDatasetId
    ? state.datasets.get(state.currentDatasetId)
    : null;

  return dataset ?? null;
}

let lastToastTime = 0;

function throttledToast(message: string) {
  const now = Date.now();

  if (now - lastToastTime < 2000) return;
  lastToastTime = now;
  toast.warning(message, { autoClose: 3000 });
}

/**
 * Bidirectional sync between left (global) and right (panel) SM visualization stores.
 * Syncs gene selection and colors. Only active when syncEnabled and rightPanelType === "sm".
 */
export function useSyncSingleMoleculeVisualization(
  rightSMVizStore: SingleMoleculeVisualizationStore,
  rightSMDatasetStore: SingleMoleculeStore,
) {
  const syncingRef = useRef(false);
  const prevLeftRef = useRef<SMSyncFields | null>(null);
  const prevRightRef = useRef<SMSyncFields | null>(null);

  const syncEnabled = useSplitScreenStore((s) => s.syncEnabled);
  const rightPanelType = useSplitScreenStore((s) => s.rightPanelType);
  const syncFromUrl = useSplitScreenStore((s) => s.syncFromUrl);
  const settlingRef = useRef(false);

  const isActive = syncEnabled && rightPanelType === "sm";

  useEffect(() => {
    if (!isActive) {
      prevLeftRef.current = null;
      prevRightRef.current = null;
      settlingRef.current = false;

      return;
    }

    const propagate = (
      sourceFields: SMSyncFields,
      prevFields: SMSyncFields | null,
      targetStore:
        | SingleMoleculeVisualizationStore
        | typeof useSingleMoleculeVisualizationStore,
      targetDatasetStore: SingleMoleculeStore | typeof useSingleMoleculeStore,
    ) => {
      if (syncingRef.current) return;
      syncingRef.current = true;
      try {
        const targetDataset = getSMDataset(targetDatasetStore);

        if (!targetDataset) return;

        const targetState = targetStore.getState();
        const prevLegend = prevFields?.selectedGenesLegend ?? new Set<string>();

        // Detect added genes
        for (const gene of sourceFields.selectedGenesLegend) {
          if (!prevLegend.has(gene)) {
            if (!targetDataset.uniqueGenes.includes(gene)) {
              throttledToast(`Gene "${gene}" not found in the other dataset`);
            } else if (!targetState.selectedGenesLegend.has(gene)) {
              // Add with same color from source
              const geneViz = sourceFields.geneDataCache.get(gene);

              targetStore
                .getState()
                .addGene(gene, geneViz?.color, geneViz?.localScale);
            }
          }
        }

        // Detect removed genes
        for (const gene of prevLegend) {
          if (!sourceFields.selectedGenesLegend.has(gene)) {
            if (targetState.selectedGenesLegend.has(gene)) {
              targetStore.getState().removeGene(gene);
            }
          }
        }

        // Detect visibility toggles (gene in legend but checked/unchecked)
        if (prevFields) {
          for (const gene of sourceFields.selectedGenesLegend) {
            if (!prevLegend.has(gene)) continue; // newly added, already handled above

            const wasVisible = prevFields.selectedGenes.has(gene);
            const isVisible = sourceFields.selectedGenes.has(gene);

            if (
              wasVisible !== isVisible &&
              targetState.selectedGenesLegend.has(gene)
            ) {
              const targetVisible = targetStore
                .getState()
                .selectedGenes.has(gene);

              if (targetVisible !== isVisible) {
                targetStore.getState().toggleGeneVisibility(gene);
              }
            }
          }
        }

        // Detect color changes for existing genes
        const prevCache = prevFields?.geneDataCache ?? new Map();

        for (const gene of sourceFields.selectedGenesLegend) {
          const prevViz = prevCache.get(gene);
          const currViz = sourceFields.geneDataCache.get(gene);

          if (
            prevViz &&
            currViz &&
            prevViz.color !== currViz.color &&
            targetState.selectedGenesLegend.has(gene)
          ) {
            targetStore.getState().setGeneColor(gene, currViz.color);
          }
        }
      } finally {
        syncingRef.current = false;
      }
    };

    // When sync was restored from URL, both panels need time to restore their
    // own URL state before sync subscriptions start propagating changes.
    const isFromUrl = syncFromUrl;

    if (isFromUrl) {
      settlingRef.current = true;
    }

    if (!isFromUrl) {
      // Manual toggle: push left panel state → right panel immediately
      const leftFields = pickSMSyncFields(
        useSingleMoleculeVisualizationStore.getState(),
      );

      propagate(leftFields, null, rightSMVizStore, rightSMDatasetStore);
    }

    // Initialize prev snapshots
    prevLeftRef.current = pickSMSyncFields(
      useSingleMoleculeVisualizationStore.getState(),
    );
    prevRightRef.current = pickSMSyncFields(rightSMVizStore.getState());

    const unsubLeft = useSingleMoleculeVisualizationStore.subscribe((state) => {
      const fields = pickSMSyncFields(state);

      if (settlingRef.current) {
        prevLeftRef.current = fields;

        return;
      }
      propagate(
        fields,
        prevLeftRef.current,
        rightSMVizStore,
        rightSMDatasetStore,
      );
      prevLeftRef.current = fields;
    });

    const unsubRight = rightSMVizStore.subscribe((state) => {
      const fields = pickSMSyncFields(state);

      if (settlingRef.current) {
        prevRightRef.current = fields;

        return;
      }
      propagate(
        fields,
        prevRightRef.current,
        useSingleMoleculeVisualizationStore,
        useSingleMoleculeStore,
      );
      prevRightRef.current = fields;
    });

    // After settling period, take fresh snapshots and enable propagation
    let settleTimer: ReturnType<typeof setTimeout> | undefined;

    if (isFromUrl) {
      settleTimer = setTimeout(() => {
        settlingRef.current = false;
        prevLeftRef.current = pickSMSyncFields(
          useSingleMoleculeVisualizationStore.getState(),
        );
        prevRightRef.current = pickSMSyncFields(rightSMVizStore.getState());
        useSplitScreenStore.getState().setSyncFromUrl(false);
      }, 3000);
    }

    return () => {
      unsubLeft();
      unsubRight();
      if (settleTimer) clearTimeout(settleTimer);
    };
  }, [isActive, rightSMVizStore, rightSMDatasetStore]);
}
