"use client";

import type { VisualizationStore } from "../stores/createVisualizationStore";
import type { DatasetStore } from "../stores/createDatasetStore";
import type { SingleMoleculeStore } from "../stores/createSingleMoleculeStore";
import type { SingleMoleculeVisualizationStore } from "../stores/createSingleMoleculeVisualizationStore";
import type { LabelledMoleculeVisualizationStore } from "../stores/createLabelledMoleculeVisualizationStore";
import type { LayersStore } from "../stores/createLayersStore";

import { createContext } from "react";

export interface PanelContextValue {
  panelId: string;
  visualizationStore: VisualizationStore;
  datasetStore: DatasetStore;
  singleMoleculeStore: SingleMoleculeStore;
  singleMoleculeVisualizationStore: SingleMoleculeVisualizationStore;
  labelledMoleculeVisualizationStore: LabelledMoleculeVisualizationStore;
  /** The layers the panel's scene draws, for its controls. */
  layersStore: LayersStore;
}

export const PanelContext = createContext<PanelContextValue | null>(null);
