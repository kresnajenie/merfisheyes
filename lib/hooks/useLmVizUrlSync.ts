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
    if (d.dv) patch.domainVariant = d.dv;
    if (d.um) patch.unselectedMode = d.um;
    if (d.ua !== undefined) patch.unselectedAlpha = d.ua;
    if (d.gs !== undefined) patch.globalScale = d.gs;
    if (d.ga !== undefined) patch.globalAlpha = d.ga;
    if (d.vm) patch.viewMode = d.vm;

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
    domainVariant,
    selections,
    geneColorSlots,
    unselectedMode,
    unselectedAlpha,
    globalScale,
    globalAlpha,
    viewMode,
  } = store;

  useEffect(() => {
    if (!datasetReady) return;

    scheduleUrlUpdate(
      "left",
      encodeLmVizState({
        colorBy,
        domainVariant,
        selections,
        geneColorSlots,
        unselectedMode,
        unselectedAlpha,
        globalScale,
        globalAlpha,
        viewMode,
      }),
    );
  }, [
    datasetReady,
    colorBy,
    domainVariant,
    selections,
    geneColorSlots,
    unselectedMode,
    unselectedAlpha,
    globalScale,
    globalAlpha,
    viewMode,
  ]);
}
