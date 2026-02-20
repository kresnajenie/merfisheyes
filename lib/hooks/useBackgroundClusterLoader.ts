"use client";

import { useEffect, useRef } from "react";

import type { StandardizedDataset } from "../StandardizedDataset";

/**
 * Background-loads remaining cluster columns after the initial priority column.
 * Non-fatal: errors are logged but don't break the UI.
 */
export function useBackgroundClusterLoader(
  dataset: StandardizedDataset | null,
  incrementClusterVersion: () => void,
) {
  const loadingRef = useRef(false);

  useEffect(() => {
    if (!dataset) return;
    if (dataset.clustersFullyLoaded) return;
    if (!dataset.adapter) return;
    if (dataset.allClusterColumnNames.length === 0) return;

    // Prevent concurrent loads for the same dataset
    if (loadingRef.current) return;

    let cancelled = false;

    const loadRemaining = async () => {
      loadingRef.current = true;

      try {
        const loadedNames = new Set(
          (dataset.clusters || []).map((c) => c.column),
        );
        const remaining = dataset.allClusterColumnNames.filter(
          (name) => !loadedNames.has(name),
        );

        if (remaining.length === 0) {
          dataset.clustersFullyLoaded = true;
          return;
        }

        console.log(
          `[BackgroundClusterLoader] Loading ${remaining.length} remaining cluster columns...`,
        );

        const newClusters = await dataset.adapter.loadClusters(remaining);

        if (cancelled) return;

        if (newClusters && newClusters.length > 0) {
          dataset.addClusters(newClusters);
          incrementClusterVersion();
        }

        dataset.clustersFullyLoaded = true;
        console.log(
          "[BackgroundClusterLoader] All cluster columns loaded.",
        );
      } catch (error) {
        if (!cancelled) {
          console.warn(
            "[BackgroundClusterLoader] Failed to load remaining clusters:",
            error,
          );
        }
      } finally {
        loadingRef.current = false;
      }
    };

    loadRemaining();

    return () => {
      cancelled = true;
    };
  }, [dataset, incrementClusterVersion]);
}
