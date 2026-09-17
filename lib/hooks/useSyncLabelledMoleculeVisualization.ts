"use client";

import type {
  LabelledMoleculeVisualizationState,
  LabelledMoleculeVisualizationStore,
  LmMenu,
} from "../stores/createLabelledMoleculeVisualizationStore";
import type { StandardizedDataset } from "../StandardizedDataset";

import { useEffect, useRef } from "react";

import { labelledMoleculeVisualizationStore } from "../stores/labelledMoleculeVisualizationStore";
import { getLmDataset } from "../stores/lmDatasetRegistry";
import {
  carrySelections,
  sameCarry,
  type LmCarryFields,
} from "../utils/lm-selection";
import { useSplitScreenStore } from "../stores/splitScreenStore";

/** What the two panels keep in step. Camera is deliberately not in here. */
type LmSyncFields = LmCarryFields & { colorBy: LmMenu };

function pick(state: LabelledMoleculeVisualizationState): LmSyncFields {
  return {
    selections: state.selections,
    hiddenValues: state.hiddenValues,
    geneColorSlots: state.geneColorSlots,
    colorOverrides: state.colorOverrides,
    colorBy: state.colorBy,
  };
}

/**
 * Bidirectional selection sync between the two labelled-molecule panels.
 *
 * Selections only — not the camera. Two embryos at different stages are
 * different sizes and orientations, so a shared pose is meaningless, and the
 * user asked for the selections instead.
 *
 * Incoming values are intersected with the target's own vocabulary, which is
 * what makes this work across the whole series without special-casing.
 * Measured over all 45 spiralia datasets, pairwise:
 *
 *   gene    96% overlap same stage, 94% across stages
 *   cell    80% same stage, 0% across
 *   domain  67% same stage, 0% across
 *
 * No cell or domain name is common to all 45 — a 2-cell embryo has AB/CD, a
 * 24-cell has 1a1 — so intersecting degrades exactly the way it should: genes
 * carry everywhere, cell and domain carry within a stage and simply drop
 * across one.
 */
export function useSyncLabelledMoleculeVisualization(
  rightStore: LabelledMoleculeVisualizationStore,
) {
  const syncingRef = useRef(false);
  const prevLeftRef = useRef<LmSyncFields | null>(null);
  const prevRightRef = useRef<LmSyncFields | null>(null);
  const settlingRef = useRef(false);

  const syncEnabled = useSplitScreenStore((s) => s.syncEnabled);
  const rightPanelType = useSplitScreenStore((s) => s.rightPanelType);
  const syncFromUrl = useSplitScreenStore((s) => s.syncFromUrl);

  const isActive = syncEnabled && rightPanelType === "lm";

  useEffect(() => {
    if (!isActive) {
      prevLeftRef.current = null;
      prevRightRef.current = null;
      settlingRef.current = false;

      return;
    }

    const propagate = (
      source: LmSyncFields,
      target: LabelledMoleculeVisualizationStore,
      targetDataset: StandardizedDataset | null,
    ) => {
      // A write to the target re-enters this subscriber; the guard is what
      // stops the two panels ping-ponging a selection back and forth.
      if (syncingRef.current) return;
      syncingRef.current = true;
      try {
        const state = target.getState();
        const carried = carrySelections(source, targetDataset);

        // Colours come across with the selection, not just which values are
        // selected: the same gene showing in two colours across a split is
        // worse than it not carrying at all.
        if (!sameCarry(state, carried)) {
          target.getState().applyUrlState(carried);
        }

        if (state.colorBy !== source.colorBy) {
          target.getState().setColorBy(source.colorBy);
        }
      } finally {
        syncingRef.current = false;
      }
    };

    // Restored from a URL: both panels are still applying their own `v=`/`rv=`
    // state, so hold off propagating until that has landed.
    const isFromUrl = syncFromUrl;

    if (isFromUrl) settlingRef.current = true;

    if (!isFromUrl) {
      // Turned on by hand: the left panel is the one the user was working in,
      // so push its selections rightward once, immediately.
      propagate(
        pick(labelledMoleculeVisualizationStore.getState()),
        rightStore,
        getLmDataset("right"),
      );
    }

    prevLeftRef.current = pick(labelledMoleculeVisualizationStore.getState());
    prevRightRef.current = pick(rightStore.getState());

    const unsubLeft = labelledMoleculeVisualizationStore.subscribe((state) => {
      const fields = pick(state);

      if (settlingRef.current) {
        prevLeftRef.current = fields;

        return;
      }
      if (
        !prevLeftRef.current ||
        !sameCarry(fields, prevLeftRef.current) ||
        fields.colorBy !== prevLeftRef.current.colorBy
      ) {
        propagate(fields, rightStore, getLmDataset("right"));
      }
      prevLeftRef.current = fields;
    });

    const unsubRight = rightStore.subscribe((state) => {
      const fields = pick(state);

      if (settlingRef.current) {
        prevRightRef.current = fields;

        return;
      }
      if (
        !prevRightRef.current ||
        !sameCarry(fields, prevRightRef.current) ||
        fields.colorBy !== prevRightRef.current.colorBy
      ) {
        propagate(
          fields,
          labelledMoleculeVisualizationStore,
          getLmDataset(null),
        );
      }
      prevRightRef.current = fields;
    });

    let settleTimer: ReturnType<typeof setTimeout> | undefined;

    if (isFromUrl) {
      settleTimer = setTimeout(() => {
        settlingRef.current = false;
        prevLeftRef.current = pick(
          labelledMoleculeVisualizationStore.getState(),
        );
        prevRightRef.current = pick(rightStore.getState());
        useSplitScreenStore.getState().setSyncFromUrl(false);
      }, 3000);
    }

    return () => {
      unsubLeft();
      unsubRight();
      if (settleTimer) clearTimeout(settleTimer);
    };
  }, [isActive, rightStore]);
}
