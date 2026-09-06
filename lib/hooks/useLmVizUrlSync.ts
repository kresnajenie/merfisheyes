"use client";

import type { LabelledMoleculeVisualizationState } from "@/lib/stores/createLabelledMoleculeVisualizationStore";

import { useEffect, useRef } from "react";

import { labelledMoleculeVisualizationStore } from "@/lib/stores/labelledMoleculeVisualizationStore";
import {
  decodeLmVizState,
  encodeLmVizState,
} from "@/lib/utils/url-visualization-state";
import {
  readUrlVizState,
  scheduleUrlUpdate,
} from "@/lib/utils/url-state-writer";

/**
 * Two-way sync between the labelled-molecule viewer state and the `v=` URL
 * param, so the Share button copies a link that reproduces the whole view.
 *
 * Reads once after the dataset is ready (selections reference values that only
 * exist once the columns are loaded), then writes on every change, debounced
 * through the shared writer that the cell and molecule viewers already use.
 */
export function useLmVizUrlSync(
  datasetReady: boolean,
  store: LabelledMoleculeVisualizationState,
) {
  const appliedRef = useRef(false);

  // ── Read: restore a shared link, once.
  useEffect(() => {
    if (!datasetReady || appliedRef.current) return;
    appliedRef.current = true;

    const encoded = readUrlVizState().left;

    if (!encoded) return;

    const d = decodeLmVizState(encoded);

    if (!d) return;

    const patch: Partial<LabelledMoleculeVisualizationState> = {};

    if (d.cb) patch.colorBy = d.cb;
    if (d.gs !== undefined) patch.globalScale = d.gs;
    if (d.ss !== undefined) patch.selectedScale = d.ss;
    if (d.us !== undefined) patch.unselectedScale = d.us;
    if (d.cam?.length === 6) {
      patch.camera = {
        position: [d.cam[0], d.cam[1], d.cam[2]],
        target: [d.cam[3], d.cam[4], d.cam[5]],
      };
      patch.pendingCamera = patch.camera;
    }

    if (d.g || d.d || d.c) {
      patch.selections = {
        gene: new Set((d.g ?? []).map(([gene]) => gene)),
        domain: new Set(d.d ?? []),
        cell: new Set(d.c ?? []),
      };
    }
    // Restore the exact palette slots so shared links keep their colours.
    if (d.g) patch.geneColorSlots = new Map(d.g);

    labelledMoleculeVisualizationStore.getState().applyUrlState(patch);
  }, [datasetReady]);

  // ── Write: mirror state into the URL.
  const {
    colorBy,
    selections,
    geneColorSlots,
    globalScale,
    selectedScale,
    unselectedScale,
    camera,
  } = store;

  useEffect(() => {
    if (!datasetReady) return;

    scheduleUrlUpdate(
      "left",
      encodeLmVizState({
        colorBy,
        selections,
        geneColorSlots,
        globalScale,
        selectedScale,
        unselectedScale,
        camera,
      }),
    );
  }, [
    datasetReady,
    colorBy,
    selections,
    geneColorSlots,
    globalScale,
    selectedScale,
    unselectedScale,
    camera,
  ]);
}
