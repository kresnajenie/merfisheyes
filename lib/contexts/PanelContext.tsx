"use client";

import { createContext } from "react";

import type { VisualizationStore } from "../stores/createVisualizationStore";
import type { DatasetStore } from "../stores/createDatasetStore";
import type { SingleMoleculeStore } from "../stores/createSingleMoleculeStore";
import type { SingleMoleculeVisualizationStore } from "../stores/createSingleMoleculeVisualizationStore";

export interface PanelContextValue {
  panelId: string;
  visualizationStore: VisualizationStore;
  datasetStore: DatasetStore;
  singleMoleculeStore: SingleMoleculeStore;
  singleMoleculeVisualizationStore: SingleMoleculeVisualizationStore;
}

export const PanelContext = createContext<PanelContextValue | null>(null);
