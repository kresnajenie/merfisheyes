"use client";

import type {
  LabelledMoleculeVisualizationState,
  LabelledMoleculeVisualizationStore,
  LmMenu,
} from "../stores/createLabelledMoleculeVisualizationStore";
import type { StandardizedDataset } from "../StandardizedDataset";

import { useEffect, useRef } from "react";

import { LM_MENUS } from "../stores/createLabelledMoleculeVisualizationStore";
import { labelledMoleculeVisualizationStore } from "../stores/labelledMoleculeVisualizationStore";
import { getLmDataset } from "../stores/lmDatasetRegistry";
import { useSplitScreenStore } from "../stores/splitScreenStore";

/** What the two panels keep in step. Camera is deliberately not in here. */
type LmSyncFields = {
  selections: Record<LmMenu, Set<string>>;
  colorBy: LmMenu;
};

function pick(state: LabelledMoleculeVisualizationState): LmSyncFields {
  return { selections: state.selections, colorBy: state.colorBy };
}

/**
 * Values of `menu` the dataset actually has, as a set for membership tests.
 *
 * A column is only present once it has been lazily loaded; until then the
 * dataset can answer "no values", which would wrongly drop every incoming
 * selection. Returning null for that case means "cannot tell yet, do nothing".
 */
function vocabulary(
  dataset: StandardizedDataset | null,
  menu: LmMenu,
): Set<string> | null {
  const col = dataset?.clusters?.find((c) => c.column === menu);

  if (!col?.uniqueValues?.length) return null;

  return new Set(col.uniqueValues);
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

        for (const menu of LM_MENUS) {
          const vocab = vocabulary(targetDataset, menu);

          // Column not loaded yet — leave this menu alone rather than clearing
          // it against an empty vocabulary.
          if (!vocab) continue;

          const next = new Set<string>();

          for (const v of source.selections[menu]) {
            if (vocab.has(v)) next.add(v);
          }

          const current = state.selections[menu];
          const same =
            current.size === next.size &&
            [...next].every((v) => current.has(v));

          if (!same) target.getState().setSelection(menu, next);
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
        fields.selections !== prevLeftRef.current?.selections ||
        fields.colorBy !== prevLeftRef.current?.colorBy
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
        fields.selections !== prevRightRef.current?.selections ||
        fields.colorBy !== prevRightRef.current?.colorBy
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
