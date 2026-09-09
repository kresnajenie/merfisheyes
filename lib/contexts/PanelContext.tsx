"use client";

import type { VisualizationStore } from "../stores/createVisualizationStore";
import type { DatasetStore } from "../stores/createDatasetStore";
import type { SingleMoleculeStore } from "../stores/createSingleMoleculeStore";
import type { SingleMoleculeVisualizationStore } from "../stores/createSingleMoleculeVisualizationStore";
import type { LabelledMoleculeVisualizationStore } from "../stores/createLabelledMoleculeVisualizationStore";

import { createContext } from "react";

export interface PanelContextValue {
  panelId: string;
  visualizationStore: VisualizationStore;
  datasetStore: DatasetStore;
  singleMoleculeStore: SingleMoleculeStore;
  singleMoleculeVisualizationStore: SingleMoleculeVisualizationStore;
  labelledMoleculeVisualizationStore: LabelledMoleculeVisualizationStore;
}

export const PanelContext = createContext<PanelContextValue | null>(null);
