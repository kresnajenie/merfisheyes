"use client";

import type { VisualizationStore } from "../stores/createVisualizationStore";
import type { VisualizationState } from "../stores/createVisualizationStore";
import type { DatasetStore } from "../stores/createDatasetStore";
import type { StandardizedDataset } from "../StandardizedDataset";

import { toast } from "react-toastify";
import { useEffect, useRef } from "react";

import { useVisualizationStore } from "../stores/visualizationStore";
import { useDatasetStore } from "../stores/datasetStore";
import { useSplitScreenStore } from "../stores/splitScreenStore";

type SyncFields = Pick<
  VisualizationState,
  "selectedColumn" | "selectedCelltypes" | "selectedGene" | "colorPalette"
>;

function pickSyncFields(state: VisualizationState): SyncFields {
  return {
    selectedColumn: state.selectedColumn,
    selectedCelltypes: state.selectedCelltypes,
    selectedGene: state.selectedGene,
    colorPalette: state.colorPalette,
  };
}

function getStandardizedDataset(
  store: DatasetStore | typeof useDatasetStore,
): StandardizedDataset | null {
  const state = store.getState();
  const dataset = state.currentDatasetId
    ? state.datasets.get(state.currentDatasetId)
    : null;

  if (dataset && "clusters" in dataset) {
    return dataset as StandardizedDataset;
  }

  return null;
}

function datasetHasColumn(
  dataset: StandardizedDataset,
  column: string,
): boolean {
  return dataset.clusters?.some((c) => c.column === column) ?? false;
}

function datasetHasGene(dataset: StandardizedDataset, gene: string): boolean {
  return dataset.genes?.includes(gene) ?? false;
}

function getColumnValues(
  dataset: StandardizedDataset,
  column: string,
  overrides: Record<string, "categorical" | "numerical"> = {},
): Set<string> {
  // If override says numerical, treat as no categorical values
  if (overrides[column] === "numerical") return new Set();

  const cluster = dataset.clusters?.find((c) => c.column === column);

  if (!cluster) return new Set();

  // If override says categorical, ignore the stored type
  if (overrides[column] === "categorical") {
    return new Set(cluster.values?.map(String) ?? []);
  }

  if (cluster.type === "numerical") return new Set();

  return new Set(cluster.values?.map(String) ?? []);
}

function isNumericalColumn(
  dataset: StandardizedDataset,
  column: string,
  overrides: Record<string, "categorical" | "numerical"> = {},
): boolean {
  if (overrides[column]) return overrides[column] === "numerical";

  const cluster = dataset.clusters?.find((c) => c.column === column);

  return cluster?.type === "numerical";
}

let lastToastTime = 0;

function throttledToast(message: string) {
  const now = Date.now();

  if (now - lastToastTime < 2000) return;
  lastToastTime = now;
  toast.warning(message, { autoClose: 3000 });
}

/**
 * Bidirectional sync between left (global) and right (panel) visualization stores.
 * Only active when syncEnabled and rightPanelType === "cell".
 */
export function useSyncVisualization(
  rightVizStore: VisualizationStore,
  rightDatasetStore: DatasetStore,
) {
  const syncingRef = useRef(false);
  const prevLeftRef = useRef<SyncFields | null>(null);
  const prevRightRef = useRef<SyncFields | null>(null);

  const syncEnabled = useSplitScreenStore((s) => s.syncEnabled);
  const rightPanelType = useSplitScreenStore((s) => s.rightPanelType);
  const syncFromUrl = useSplitScreenStore((s) => s.syncFromUrl);
  const settlingRef = useRef(false);

  const isActive = syncEnabled && rightPanelType === "cell";

  useEffect(() => {
    if (!isActive) {
      prevLeftRef.current = null;
      prevRightRef.current = null;
      settlingRef.current = false;

      return;
    }

    const propagate = (
      sourceFields: SyncFields,
      prevFields: SyncFields | null,
      targetStore: VisualizationStore | typeof useVisualizationStore,
      targetDatasetStore: DatasetStore | typeof useDatasetStore,
    ) => {
      if (syncingRef.current) return;
      syncingRef.current = true;
      try {
        const targetDataset = getStandardizedDataset(targetDatasetStore);

        if (!targetDataset) return;

        const targetState = targetStore.getState();
        const overrides = targetState.columnTypeOverrides ?? {};

        // 1. selectedColumn changed
        if (sourceFields.selectedColumn !== prevFields?.selectedColumn) {
          const column = sourceFields.selectedColumn;

          if (column) {
            if (!datasetHasColumn(targetDataset, column)) {
              throttledToast(
                `Column "${column}" not found in the other dataset`,
              );
            } else {
              const numerical = isNumericalColumn(targetDataset, column, overrides);

              targetStore.getState().setSelectedColumn(column, numerical);

              // Re-apply celltypes after column change (setSelectedColumn clears them)
              if (!numerical) {
                const targetValues = getColumnValues(targetDataset, column, overrides);

                for (const ct of sourceFields.selectedCelltypes) {
                  if (targetValues.has(ct)) {
                    targetStore.getState().toggleCelltype(ct);
                  }
                }
              }
            }
          } else {
            targetStore.getState().setSelectedColumn(null);
          }
        }

        // 2. selectedCelltypes changed (but column didn't change)
        else if (
          !setsEqual(
            sourceFields.selectedCelltypes,
            prevFields?.selectedCelltypes,
          )
        ) {
          const column = sourceFields.selectedColumn;

          if (column && datasetHasColumn(targetDataset, column)) {
            const targetValues = getColumnValues(targetDataset, column, overrides);
            const prevCelltypes = prevFields?.selectedCelltypes ?? new Set();
            const currentTargetCelltypes = targetState.selectedCelltypes;

            // Find added celltypes
            for (const ct of sourceFields.selectedCelltypes) {
              if (!prevCelltypes.has(ct)) {
                if (!targetValues.has(ct)) {
                  throttledToast(
                    `Celltype "${ct}" not found in the other dataset`,
                  );
                } else if (!currentTargetCelltypes.has(ct)) {
                  targetStore.getState().toggleCelltype(ct);
                }
              }
            }

            // Find removed celltypes
            for (const ct of prevCelltypes) {
              if (!sourceFields.selectedCelltypes.has(ct)) {
                if (currentTargetCelltypes.has(ct)) {
                  targetStore.getState().toggleCelltype(ct);
                }
              }
            }
          }
        }

        // 3. selectedGene changed
        if (sourceFields.selectedGene !== prevFields?.selectedGene) {
          const gene = sourceFields.selectedGene;

          if (gene) {
            if (!datasetHasGene(targetDataset, gene)) {
              throttledToast(`Gene "${gene}" not found in the other dataset`);
            } else {
              targetStore.getState().setSelectedGene(gene);
            }
          } else {
            targetStore.getState().setSelectedGene(null);
          }
        }

        // 4. colorPalette changed
        if (
          !shallowEqual(sourceFields.colorPalette, prevFields?.colorPalette)
        ) {
          const column = sourceFields.selectedColumn;

          if (column && datasetHasColumn(targetDataset, column)) {
            const targetValues = getColumnValues(targetDataset, column, overrides);
            // Only sync colors for celltypes that exist in both datasets
            const currentTargetPalette = targetState.colorPalette;
            const newPalette = { ...currentTargetPalette };
            let changed = false;

            for (const [key, value] of Object.entries(
              sourceFields.colorPalette,
            )) {
              if (targetValues.has(key) && newPalette[key] !== value) {
                newPalette[key] = value;
                changed = true;
              }
            }
            if (changed) {
              targetStore.getState().setColorPalette(newPalette);

              // Also update in-memory palette on the target dataset
              const targetCluster = targetDataset.clusters?.find(
                (c) => c.column === column,
              );

              if (targetCluster) {
                targetCluster.palette = newPalette;
              }
            }
          }
        }
      } finally {
        syncingRef.current = false;
      }
    };

    // When sync was restored from URL, both panels need time to restore their
    // own URL state before sync subscriptions start propagating changes.
    // During settling, subscriptions only update prev-state snapshots.
    const isFromUrl = syncFromUrl;

    if (isFromUrl) {
      settlingRef.current = true;
    }

    if (!isFromUrl) {
      // Manual toggle: push left panel state → right panel immediately
      const leftFields = pickSyncFields(useVisualizationStore.getState());

      propagate(leftFields, null, rightVizStore, rightDatasetStore);
    }

    // Initialize prev snapshots
    prevLeftRef.current = pickSyncFields(useVisualizationStore.getState());
    prevRightRef.current = pickSyncFields(rightVizStore.getState());

    const unsubLeft = useVisualizationStore.subscribe((state) => {
      const fields = pickSyncFields(state);

      if (settlingRef.current) {
        // During settling, only track state — don't propagate
        prevLeftRef.current = fields;

        return;
      }
      propagate(fields, prevLeftRef.current, rightVizStore, rightDatasetStore);
      prevLeftRef.current = fields;
    });

    const unsubRight = rightVizStore.subscribe((state) => {
      const fields = pickSyncFields(state);

      if (settlingRef.current) {
        // During settling, only track state — don't propagate
        prevRightRef.current = fields;

        return;
      }
      propagate(
        fields,
        prevRightRef.current,
        useVisualizationStore,
        useDatasetStore,
      );
      prevRightRef.current = fields;
    });

    // After settling period, take fresh snapshots and enable propagation
    let settleTimer: ReturnType<typeof setTimeout> | undefined;

    if (isFromUrl) {
      settleTimer = setTimeout(() => {
        settlingRef.current = false;
        prevLeftRef.current = pickSyncFields(useVisualizationStore.getState());
        prevRightRef.current = pickSyncFields(rightVizStore.getState());
        useSplitScreenStore.getState().setSyncFromUrl(false);
      }, 3000);
    }

    return () => {
      unsubLeft();
      unsubRight();
      if (settleTimer) clearTimeout(settleTimer);
    };
  }, [isActive, rightVizStore, rightDatasetStore]);
}

function setsEqual(a: Set<string>, b: Set<string> | undefined): boolean {
  if (!b) return false;
  if (a.size !== b.size) return false;
  for (const item of a) {
    if (!b.has(item)) return false;
  }

  return true;
}

function shallowEqual(
  a: Record<string, string>,
  b: Record<string, string> | undefined,
): boolean {
  if (!b) return false;
  const aKeys = Object.keys(a);
  const bKeys = Object.keys(b);

  if (aKeys.length !== bKeys.length) return false;
  for (const key of aKeys) {
    if (a[key] !== b[key]) return false;
  }

  return true;
}
