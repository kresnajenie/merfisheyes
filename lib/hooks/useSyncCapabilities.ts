"use client";

import type { DatasetStore } from "../stores/createDatasetStore";
import type { LayersStore } from "../stores/createLayersStore";
import type { SingleMoleculeStore } from "../stores/createSingleMoleculeStore";
import type { PanelContent } from "../stores/splitScreenStore";

import { usePathname } from "next/navigation";
import { useStore } from "zustand";

import { useDatasetStore } from "../stores/datasetStore";
import { layersStore } from "../stores/layersStore";
import { useSingleMoleculeStore } from "../stores/singleMoleculeStore";
import {
  hasCells,
  hasImage,
  hasMolecules,
  panelTypeForViewer,
  syncCapabilities,
  useSplitScreenStore,
} from "../stores/splitScreenStore";

/**
 * What the split's two panels can sync, from what each has loaded: the left
 * panel's global stores and the right panel's own. Both the sync button and
 * the sync hooks read this, so they always agree.
 */
export function useSyncCapabilities(right: {
  datasetStore: DatasetStore;
  singleMoleculeStore: SingleMoleculeStore;
  layersStore: LayersStore;
}): PanelContent {
  const pathname = usePathname();
  const rightType = useSplitScreenStore((s) => s.rightPanelType);

  const left: PanelContent = {
    cells: useDatasetStore(hasCells),
    molecules: useSingleMoleculeStore(hasMolecules),
    // The labelled-molecule viewer keeps its own stores, so it is known by
    // the page (left) and by what the right panel was opened as
    labelled: panelTypeForViewer(pathname.split("/")[1] ?? "") === "lm",
    image: useStore(layersStore, hasImage),
  };
  const rightContent: PanelContent = {
    cells: useStore(right.datasetStore, hasCells),
    molecules: useStore(right.singleMoleculeStore, hasMolecules),
    labelled: rightType === "lm",
    image: useStore(right.layersStore, hasImage),
  };

  return syncCapabilities(left, rightContent);
}
