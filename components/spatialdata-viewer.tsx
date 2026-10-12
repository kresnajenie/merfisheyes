"use client";

import type { SingleMoleculeDataset } from "@/lib/SingleMoleculeDataset";
import type { StandardizedDataset } from "@/lib/StandardizedDataset";

import { useCallback, useEffect, useState } from "react";
import { Progress, Spinner } from "@heroui/react";

import { ImageScene } from "@/components/image-viewer";
import { SingleMoleculeControls } from "@/components/single-molecule-controls";
import { SingleMoleculeLegends } from "@/components/single-molecule-legends";
import { SingleMoleculeThreeScene } from "@/components/single-molecule-three-scene";
import { SpatialDataImageControls } from "@/components/spatialdata-image-controls";
import { type SceneHandle, ThreeScene } from "@/components/three-scene";
import { VisualizationControls } from "@/components/visualization-controls";
import { SplitScreenContainer } from "@/components/split-screen-container";
import { subtitle, title } from "@/components/primitives";
import { ImageLayer } from "@/lib/spatialdata/ImageLayer";
import {
  openSpatialDataImage,
  type SpatialDataImage,
} from "@/lib/spatialdata/image";
import { loadSpatialDataCells } from "@/lib/spatialdata/loadSpatialDataCells";
import { loadSpatialDataMolecules } from "@/lib/spatialdata/loadSpatialDataMolecules";
import {
  elementNames,
  openSpatialDataStore,
  openSpatialDataStoreFromUrl,
} from "@/lib/spatialdata/store";
import { useSMOverlayUrlSync } from "@/lib/hooks/useUrlVizSync";
import { useDatasetStore } from "@/lib/stores/datasetStore";
import { useSingleMoleculeStore } from "@/lib/stores/singleMoleculeStore";
import { useSingleMoleculeVisualizationStore } from "@/lib/stores/singleMoleculeVisualizationStore";
import { useVisualizationStore } from "@/lib/stores/visualizationStore";
import { applyCellOpenState } from "@/lib/viewer/open-dataset";

/**
 * Viewer for a SpatialData zarr store, read lazily and as written — either a
 * dataset on our S3 (`datasetId`) or a store at a public URL (`url`).
 *
 * A store may hold cells (its table), molecules (a linked single molecule
 * dataset, converted from its points with scripts/spatialdata_points_to_sm.py
 * and named by the mapping.json beside it) and an image, in any combination:
 *   - with cells: the single cell scene, molecules overlaid, image underneath
 *   - molecules but no cells: the single molecule scene, image underneath
 *   - an image alone: just the image
 */
export function SpatialDataViewer({
  datasetId,
  url,
}: {
  datasetId?: string;
  url?: string;
}) {
  const vizStore = useVisualizationStore();
  const addDataset = useDatasetStore((s) => s.addDataset);
  const addMolecules = useSingleMoleculeStore((s) => s.addDataset);
  const [dataset, setDataset] = useState<StandardizedDataset | null>(null);
  // Molecules shown on their own, when the store has no cells
  const [molecules, setMolecules] = useState<SingleMoleculeDataset | null>(
    null,
  );
  const [image, setImage] = useState<SpatialDataImage | null>(null);
  const [loaded, setLoaded] = useState(false);
  const [imageLayer, setImageLayer] = useState<ImageLayer | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [progress, setProgress] = useState(0);
  const [message, setMessage] = useState("Initializing...");

  // The scene loads the linked molecule dataset into the global SM store;
  // this picks its default genes (or restores them from the URL).
  const smOverlayDataset = useSingleMoleculeStore((s) =>
    s.currentDatasetId ? (s.datasets.get(s.currentDatasetId) ?? null) : null,
  );
  const smVizStore = useSingleMoleculeVisualizationStore();

  useSMOverlayUrlSync(!!smOverlayDataset, smOverlayDataset, smVizStore);

  useEffect(() => {
    let cancelled = false;

    (async () => {
      try {
        const onProgress = (p: number, m: string) => {
          setProgress(p);
          setMessage(m);
        };
        let store;
        let id: string;
        let name: string;
        let folder: string | undefined;

        setMessage("Opening store...");
        if (url) {
          const base = url.replace(/\/+$/, "");

          folder = base.slice(0, base.lastIndexOf("/"));
          name = base.slice(folder.length + 1).replace(/\.zarr$/, "");
          id = `spatialdata_${name}`;
          store = await openSpatialDataStoreFromUrl(base);
        } else if (datasetId) {
          const metaRes = await fetch(`/api/datasets/${datasetId}`);
          const meta = await metaRes.json().catch(() => ({}));

          if (metaRes.status !== 200) {
            throw new Error(
              meta.message || `Dataset fetch failed: ${metaRes.status}`,
            );
          }
          id = datasetId;
          name = meta.title || datasetId;
          store = await openSpatialDataStore(datasetId);
        } else {
          throw new Error("No dataset given");
        }

        let cells: StandardizedDataset | null = null;
        let ownMolecules: SingleMoleculeDataset | null = null;

        if (elementNames(store.keys, "tables").length > 0) {
          cells = await loadSpatialDataCells(
            store,
            id,
            name,
            onProgress,
            folder,
          );
          if (cancelled) return;
          await applyCellOpenState({
            dataset: cells,
            config: null,
            store: vizStore,
          });
        } else {
          setMessage("Loading molecules...");
          ownMolecules = await loadSpatialDataMolecules(store);
        }
        // Only the image's metadata; a store without one shows the rest
        const img = await openSpatialDataImage(store).catch((e) => {
          console.error("Error opening image:", e);

          return null;
        });

        if (cancelled) return;
        if (!cells && !ownMolecules && !img) {
          throw new Error(
            "This SpatialData store has no cells, molecules or image to show",
          );
        }
        setImage(img);
        if (cells) {
          addDataset(cells);
          setDataset(cells);
        }
        if (ownMolecules) {
          addMolecules(ownMolecules);
          setMolecules(ownMolecules);
        }
        setLoaded(true);
      } catch (err) {
        console.error("Error loading SpatialData store:", err);
        if (!cancelled) {
          setError(err instanceof Error ? err.message : "Failed to load");
        }
      }
    })();

    return () => {
      cancelled = true;
    };
  }, [datasetId, url]);

  // The scene is rebuilt on view changes; the layer lives and dies with it.
  const addImageLayer = useCallback(
    ({ group, camera, renderer }: SceneHandle) => {
      if (!image) return;
      const layer = new ImageLayer(image, group, camera, renderer);

      setImageLayer(layer);

      return () => {
        layer.dispose();
        setImageLayer(null);
      };
    },
    [image],
  );

  if (error) {
    return (
      <div className="flex flex-col items-center justify-center h-full gap-4 p-8 text-center">
        <h2 className={title({ size: "md", color: "pink" })}>
          Failed to load dataset
        </h2>
        <p className={subtitle()}>{error}</p>
      </div>
    );
  }

  if (!loaded) {
    return (
      <div className="flex items-center justify-center h-full">
        <div className="flex flex-col items-center gap-4 w-full max-w-md px-4">
          <Spinner color="primary" size="lg" />
          <p className={subtitle({ class: "text-center" })}>
            Loading SpatialData store...
          </p>
          <Progress
            aria-label="Loading progress"
            className="w-full"
            color="primary"
            value={progress}
          />
          <p className="text-sm text-default-500 text-center">{message}</p>
        </div>
      </div>
    );
  }

  if (dataset) {
    return (
      <SplitScreenContainer>
        <VisualizationControls />
        <ThreeScene dataset={dataset} onSceneReady={addImageLayer} />
        {imageLayer && <SpatialDataImageControls layer={imageLayer} />}
      </SplitScreenContainer>
    );
  }
  if (molecules) {
    return (
      <SplitScreenContainer>
        <SingleMoleculeControls />
        <SingleMoleculeLegends />
        <SingleMoleculeThreeScene onSceneReady={addImageLayer} />
        {imageLayer && <SpatialDataImageControls layer={imageLayer} />}
      </SplitScreenContainer>
    );
  }

  return image && <ImageScene image={image} />;
}
