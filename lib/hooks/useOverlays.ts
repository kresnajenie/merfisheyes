"use client";

import type { OpenedElements } from "@/lib/scene/registry";
import type { SingleMoleculeDataset } from "@/lib/SingleMoleculeDataset";
import type { StandardizedDataset } from "@/lib/StandardizedDataset";
import type { DatasetLinkConfig } from "@/lib/config/dataset-links";

import { useEffect, useState } from "react";
import { toast } from "react-toastify";

import {
  usePanelSingleMoleculeApi,
  usePanelSingleMoleculeVisualizationApi,
  usePanelVisualizationStore,
} from "./usePanelStores";

import {
  fetchMappingConfig,
  getDatasetLinkConfig,
} from "@/lib/config/dataset-links";
import { openImageRef } from "@/lib/spatialdata/image";
import { useLocalBundleStore } from "@/lib/stores/localBundleStore";
import { useSingleMoleculeStore } from "@/lib/stores/singleMoleculeStore";
import { loadClusterColumn } from "@/lib/utils/load-cluster-column";
import { type OverlayIO, resolveOverlays } from "@/lib/viewer/resolve-overlays";

const io: OverlayIO = {
  fetchMapping: fetchMappingConfig,
  // Revalidated: the route is cached for a minute, which would hide an
  // overlay that was just attached or removed until the cache expired
  fetchAppMapping: (id) =>
    fetch(`/api/datasets/${id}/mapping`, { cache: "no-cache" })
      .then((r) => (r.ok ? (r.json() as Promise<DatasetLinkConfig>) : null))
      .catch(() => null),
  openImage: openImageRef,
  loadMolecules: async (ref, source, linkFormat) => {
    const { openLinkedMolecules } = await import(
      "@/lib/spatialdata/loadSpatialDataMolecules"
    );

    return openLinkedMolecules(ref, source, linkFormat);
  },
  loadColumn: loadClusterColumn,
};

/**
 * The overlays of a cell dataset (lib/viewer/resolve-overlays.ts), kept in
 * step with this panel's stores: linked molecules go into the panel's
 * molecule store while they are shown and leave it, with their gene
 * selection, when the dataset changes.
 */
export function useOverlays(dataset: StandardizedDataset | null): {
  elements: OpenedElements;
  linkConfig: DatasetLinkConfig | null;
  loading: boolean;
} {
  const smStoreApi = usePanelSingleMoleculeApi();
  const smVizApi = usePanelSingleMoleculeVisualizationApi();
  const incrementClusterVersion = usePanelVisualizationStore(
    (s) => s.incrementClusterVersion,
  );
  // Molecules picked to check a local preview against: drawn once picked
  const check = useLocalBundleStore((s) =>
    dataset && s.anchorId === dataset.id ? s.check : null,
  );
  const [elements, setElements] = useState<OpenedElements>({});
  const [linkConfig, setLinkConfig] = useState<DatasetLinkConfig | null>(null);
  const [loading, setLoading] = useState(false);

  useEffect(() => {
    setElements({});
    setLinkConfig(dataset ? getDatasetLinkConfig(dataset) : null);
    setLoading(false);
    if (!dataset) return;
    let stale = false;
    let molecules: SingleMoleculeDataset | null = null;
    const local = useLocalBundleStore.getState();

    resolveOverlays(
      dataset,
      local.anchorId === dataset.id
        ? {
            image: local.image,
            shapes: local.shapes,
            molecules: local.molecules,
            check,
          }
        : null,
      getDatasetLinkConfig(dataset),
      io,
      {
        element: (patch) => {
          setElements((e) => ({ ...e, ...patch }));
          const ds = patch.points ?? patch.molecules;

          if (ds) {
            molecules = ds;
            smStoreApi.getState().addDataset(ds);
          }
        },
        link: setLinkConfig,
        loading: setLoading,
        // The saved-default-genes fallback keys off a public URL, read by
        // the main (left) viewer's page from the global store
        source: (url) => {
          if (smStoreApi === useSingleMoleculeStore) {
            useSingleMoleculeStore.getState().setOverlaySourceUrl(url);
          }
        },
        column: incrementClusterVersion,
        warn: (message) =>
          toast.warn(message, {
            autoClose: 12000,
            toastId: "local-image-normalized",
          }),
        stale: () => stale,
      },
    );

    return () => {
      stale = true;
      // The gene picker and legends are gated on the selection, not on the
      // dataset: both go with the overlay
      if (molecules) {
        smVizApi.getState().clearGenes();
        smStoreApi.getState().removeDataset(molecules.id);
      }
    };
  }, [dataset, check]);

  return { elements, linkConfig, loading };
}
