"use client";

import type { PanelType } from "@/lib/stores/splitScreenStore";

import { ThreeScene } from "./three-scene";
import { VisualizationControls } from "./visualization-controls";
import UMAPPanel from "./umap-panel";
import { SingleMoleculeThreeScene } from "./single-molecule-three-scene";
import { SingleMoleculeControls } from "./single-molecule-controls";
import { SingleMoleculeLegends } from "./single-molecule-legends";
import LabelledMoleculeViewerPage from "./labelled-molecule-viewer";
import {
  DatasetErrorScreen,
  DatasetLoadingScreen,
} from "./dataset-loading-screen";

import {
  useCellVizUrlSync,
  useSMOverlayUrlSync,
  useSMVizUrlSync,
} from "@/lib/hooks/useUrlVizSync";
import { useDatasetLoader } from "@/lib/hooks/useDatasetLoader";
import {
  usePanelVisualizationStore,
  usePanelSingleMoleculeStore,
  usePanelSingleMoleculeVisualizationStore,
} from "@/lib/hooks/usePanelStores";

interface SplitPanelViewerProps {
  datasetId: string | null;
  s3Url: string | null;
  type: PanelType;
}

export function SplitPanelViewer({
  datasetId,
  s3Url,
  type,
}: SplitPanelViewerProps) {
  if (type === "cell") {
    return <CellViewer datasetId={datasetId} s3Url={s3Url} />;
  }
  if (type === "lm") {
    // The LM viewer already owns its loader, and inside a PanelProvider its
    // hooks resolve to this panel's stores, so it needs no separate copy —
    // only the URL to load and the flag that strips the page-level chrome.
    return s3Url ? <LabelledMoleculeViewerPage embedded s3Url={s3Url} /> : null;
  }

  return <SingleMoleculeViewer datasetId={datasetId} s3Url={s3Url} />;
}

function CellViewer({
  datasetId,
  s3Url,
}: {
  datasetId: string | null;
  s3Url: string | null;
}) {
  const vizStore = usePanelVisualizationStore();
  // The same opener as the full page, into this panel's stores; a split
  // panel does not count as a view of the dataset
  const { loaded, loading, overlaysLoading, progress, message, error } =
    useDatasetLoader({ kind: "cell", datasetId, s3Url }, { countViews: false });
  const dataset = loaded?.cells ?? null;

  // URL sync hook (handles reading after the dataset is ready + writing)
  useCellVizUrlSync(!!dataset, dataset, vizStore);

  // SM overlay of this panel: three-scene loads the linked molecule dataset
  // into the panel's SM store; this restores / seeds its genes (rov=).
  const smOverlayDataset = usePanelSingleMoleculeStore((s) =>
    s.currentDatasetId ? (s.datasets.get(s.currentDatasetId) ?? null) : null,
  );
  const smOverlayVizStore = usePanelSingleMoleculeVisualizationStore();

  useSMOverlayUrlSync(!!smOverlayDataset, smOverlayDataset, smOverlayVizStore);

  if (loading) {
    return (
      <DatasetLoadingScreen
        message={message}
        progress={progress}
        variant="panel"
      />
    );
  }
  if (error) return <DatasetErrorScreen error={error} variant="panel" />;
  if (!dataset) return null;

  return (
    <>
      <VisualizationControls />
      <ThreeScene
        dataset={dataset}
        elements={loaded?.elements}
        linkConfig={loaded?.linkConfig}
        loadingOverlay={overlaysLoading}
      />
      <UMAPPanel />
    </>
  );
}

function SingleMoleculeViewer({
  datasetId,
  s3Url,
}: {
  datasetId: string | null;
  s3Url: string | null;
}) {
  const smVizStore = usePanelSingleMoleculeVisualizationStore();
  const { loaded, loading, progress, message, error } = useDatasetLoader(
    { kind: "sm", datasetId, s3Url },
    { countViews: false },
  );
  const smDataset = loaded?.molecules ?? null;

  // URL sync hook (handles reading after the dataset is ready + writing)
  useSMVizUrlSync(!!smDataset, smDataset, smVizStore);

  if (loading) {
    return (
      <DatasetLoadingScreen
        message={message}
        progress={progress}
        variant="panel"
      />
    );
  }
  if (error) return <DatasetErrorScreen error={error} variant="panel" />;
  if (!smDataset) return null;

  return (
    <>
      <SingleMoleculeControls />
      <SingleMoleculeLegends />
      <SingleMoleculeThreeScene />
    </>
  );
}
