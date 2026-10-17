"use client";

import type { DatasetSource, LoadedDataset } from "@/lib/viewer/load-dataset";

import { useEffect, useMemo, useState } from "react";

import { useOverlays } from "./useOverlays";
import {
  usePanelDatasetStore,
  usePanelId,
  usePanelSingleMoleculeStore,
  usePanelSingleMoleculeVisualizationStore,
  usePanelVisualizationStore,
} from "./usePanelStores";

import { useDatasetStore } from "@/lib/stores/datasetStore";
import { useSingleMoleculeStore } from "@/lib/stores/singleMoleculeStore";
import { useViewerRegistrationStore } from "@/lib/stores/viewerRegistrationStore";
import { loadDataset } from "@/lib/viewer/load-dataset";

export interface DatasetLoader {
  loaded: LoadedDataset | null;
  loading: boolean;
  /** The linked molecules of a cell dataset are still being read. */
  overlaysLoading: boolean;
  progress: number;
  message: string;
  error: string | null;
  retry: () => void;
}

/**
 * Open a dataset into this panel's stores, with the progress and the
 * failure a loading screen shows. The main (left) panel also publishes the
 * dataset's registration, for the owner controls and the claim banner, and
 * counts one view per load for a registered dataset.
 */
export function useDatasetLoader(
  source: DatasetSource | null,
  options: { countViews?: boolean } = {},
): DatasetLoader {
  const panel = usePanelId() === "right" ? "right" : "left";
  const viz = usePanelVisualizationStore();
  const smViz = usePanelSingleMoleculeVisualizationStore();
  const addCells = usePanelDatasetStore((s) => s.addDataset);
  const addMolecules = usePanelSingleMoleculeStore((s) => s.addDataset);
  const setRegistration = useViewerRegistrationStore((s) => s.set);
  const resetRegistration = useViewerRegistrationStore((s) => s.reset);
  const [loaded, setLoaded] = useState<LoadedDataset | null>(null);
  const [loading, setLoading] = useState(true);
  const [progress, setProgress] = useState(0);
  const [message, setMessage] = useState("Initializing...");
  const [error, setError] = useState<string | null>(null);
  const [attempt, setAttempt] = useState(0);
  // A cell dataset's overlays (a mapping.json, a dropped folder) are
  // resolved once it is open and join its elements as they arrive
  const overlays = useOverlays(
    source?.kind === "cell" ? (loaded?.cells ?? null) : null,
  );
  const withOverlays = useMemo(
    () =>
      loaded && source?.kind === "cell"
        ? {
            ...loaded,
            elements: overlays.elements,
            linkConfig: overlays.linkConfig,
          }
        : loaded,
    [loaded, overlays.elements, overlays.linkConfig, source?.kind],
  );
  const key = source
    ? `${source.kind}|${source.datasetId ?? ""}|${source.s3Url ?? ""}`
    : null;

  useEffect(() => {
    if (!source || (!source.datasetId && !source.s3Url)) {
      setError("No dataset given");
      setLoading(false);

      return;
    }
    let cancelled = false;

    setLoading(true);
    setError(null);
    setProgress(0);
    setMessage("Initializing...");
    loadDataset(source, {
      panel,
      viz,
      smViz,
      // Already open in this tab (just uploaded, or shown in the other panel)
      findCells: (id) => {
        const ds = useDatasetStore.getState().datasets.get(id);

        return ds && "spatial" in ds ? ds : undefined;
      },
      findMolecules: (id) => {
        const ds = useSingleMoleculeStore.getState().datasets.get(id);

        return ds && "uniqueGenes" in ds ? ds : undefined;
      },
      onProgress: (p, m) => {
        if (cancelled) return;
        setProgress(p);
        setMessage(m);
      },
      isCancelled: () => cancelled,
    })
      .then((result) => {
        if (cancelled || !result) return;
        if (result.cells) addCells(result.cells);
        if (result.molecules) addMolecules(result.molecules);
        // A store's own molecules over its cells: the gene picker and the
        // legends read them from the panel's molecule store
        const over = result.elements.points ?? result.elements.molecules;

        if (over) addMolecules(over);
        if (panel === "left") setRegistration(result.registration);
        setLoaded(result);
        setLoading(false);
        // The server dedups per session and day
        if (options.countViews !== false && result.registration.dbId) {
          fetch(`/api/datasets/${result.registration.dbId}/view`, {
            method: "POST",
          }).catch(() => {});
        }
      })
      .catch((err) => {
        if (cancelled) return;
        console.error("Error loading dataset:", err);
        setError(err instanceof Error ? err.message : "Failed to load dataset");
        setLoading(false);
      });

    return () => {
      cancelled = true;
    };
  }, [key, attempt]);

  // Leaving the viewer clears what the owner controls read
  useEffect(
    () => () => {
      if (panel === "left") resetRegistration();
    },
    [],
  );

  return {
    loaded: withOverlays,
    loading,
    overlaysLoading: overlays.loading,
    progress,
    message,
    error,
    retry: () => setAttempt((n) => n + 1),
  };
}
