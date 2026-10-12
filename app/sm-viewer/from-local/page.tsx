"use client";

import { Suspense, useCallback, useEffect, useState } from "react";
import { useRouter } from "next/navigation";

import { SingleMoleculeThreeScene } from "@/components/single-molecule-three-scene";
import { SpatialDataImageControls } from "@/components/spatialdata-image-controls";
import { type SceneHandle } from "@/components/three-scene";
import { ImageLayer } from "@/lib/spatialdata/ImageLayer";
import { useLocalBundleStore } from "@/lib/stores/localBundleStore";
import { SingleMoleculeControls } from "@/components/single-molecule-controls";
import { SingleMoleculeLegends } from "@/components/single-molecule-legends";
import { SplitScreenContainer } from "@/components/split-screen-container";
import { useSingleMoleculeStore } from "@/lib/stores/singleMoleculeStore";
import { useSingleMoleculeVisualizationStore } from "@/lib/stores/singleMoleculeVisualizationStore";
import { useSMVizUrlSync } from "@/lib/hooks/useUrlVizSync";
import { pickDefaultGenes } from "@/lib/utils/auto-select-genes";

function ViewerFromLocalContent() {
  const router = useRouter();
  const dataset = useSingleMoleculeStore((state) => {
    const id = state.currentDatasetId;

    return id ? state.datasets.get(id) : null;
  });

  const smVizStore = useSingleMoleculeVisualizationStore();
  const { addGene, clearGenes } = smVizStore;

  // URL visualization state sync
  const { hasUrlStateRef } = useSMVizUrlSync(
    !!dataset,
    dataset ?? null,
    smVizStore,
  );

  // Auto-select default genes when dataset loads (skip if URL state was applied)
  useEffect(() => {
    if (!dataset || hasUrlStateRef.current) return;

    clearGenes();

    const genesToSelect = pickDefaultGenes(dataset.uniqueGenes);

    genesToSelect.forEach((gene) => {
      addGene(gene);
    });
  }, [dataset, addGene, clearGenes]);

  useEffect(() => {
    if (!dataset) {
      router.replace("/?mode=sm");
    }
  }, [dataset, router]);

  // An image dropped in the same folder as these molecules goes under them.
  // The scene is rebuilt on view changes; the layer lives and dies with it.
  const image = useLocalBundleStore((s) =>
    dataset && s.anchorId === dataset.id ? s.image : null,
  );
  const [imageLayer, setImageLayer] = useState<ImageLayer | null>(null);
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

  if (!dataset) {
    return null;
  }

  return (
    <SplitScreenContainer>
      <SingleMoleculeControls />
      <SingleMoleculeLegends />
      <SingleMoleculeThreeScene onSceneReady={addImageLayer} />
      {imageLayer && <SpatialDataImageControls layer={imageLayer} />}
    </SplitScreenContainer>
  );
}

export default function SingleMoleculeViewerFromLocalPage() {
  return (
    <Suspense fallback={<div>Loading...</div>}>
      <ViewerFromLocalContent />
    </Suspense>
  );
}
