"use client";

import type { StandardizedDataset } from "@/lib/StandardizedDataset";
import type { PanelType } from "@/lib/stores/splitScreenStore";

import { Suspense, useEffect, useRef, useState } from "react";
import { useRouter, useSearchParams } from "next/navigation";
import { Button, Progress, Spinner } from "@heroui/react";

import { ThreeScene } from "@/components/three-scene";
import { VisualizationControls } from "@/components/visualization-controls";
import UMAPPanel from "@/components/umap-panel";
import { SplitScreenContainer } from "@/components/split-screen-container";
import { ClaimDatasetBanner } from "@/components/claim-dataset-banner";
import { useVisualizationStore } from "@/lib/stores/visualizationStore";
import { useDatasetStore } from "@/lib/stores/datasetStore";
import { useSingleMoleculeStore } from "@/lib/stores/singleMoleculeStore";
import { useSingleMoleculeVisualizationStore } from "@/lib/stores/singleMoleculeVisualizationStore";
import { useSplitScreenStore } from "@/lib/stores/splitScreenStore";
import { useViewerRegistrationStore } from "@/lib/stores/viewerRegistrationStore";
import { selectBestClusterColumn } from "@/lib/utils/dataset-utils";
import {
  applyViewerConfig,
  type ViewerConfig,
} from "@/lib/utils/viewer-config";
import { loadClusterColumn } from "@/lib/utils/load-cluster-column";
import {
  useCellVizUrlSync,
  useSMOverlayUrlSync,
} from "@/lib/hooks/useUrlVizSync";
import LightRays from "@/components/react-bits/LightRays";
import { subtitle, title } from "@/components/primitives";

function ViewerFromS3Content() {
  const searchParams = useSearchParams();
  const router = useRouter();
  const vizStore = useVisualizationStore();
  const { addDataset } = useDatasetStore();
  const {
    isSplitMode,
    rightPanelDatasetId,
    rightPanelS3Url,
    rightPanelType,
    syncEnabled,
    enableSplit,
    setRightPanel,
    setRightPanelS3,
    setSyncEnabled,
    setSyncFromUrl,
  } = useSplitScreenStore();
  const [dataset, setDataset] = useState<StandardizedDataset | null>(null);
  const [isLoading, setIsLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [loadingProgress, setLoadingProgress] = useState(0);
  const [loadingMessage, setLoadingMessage] = useState("Initializing...");
  // DB registration for this S3 URL (null until the by-url lookup resolves).
  const [registration, setRegistration] = useState<{
    id: string;
    ownerId: string | null;
    adminOwned: boolean;
    viewerConfig: ViewerConfig | null;
  } | null>(null);
  const setViewerRegistration = useViewerRegistrationStore((s) => s.set);
  const resetViewerRegistration = useViewerRegistrationStore((s) => s.reset);
  const defaultsAppliedRef = useRef(false);
  const viewCountedRef = useRef(false);

  const s3Url = searchParams.get("url");

  // URL visualization state sync (SC viz writes to v=)
  const { hasUrlStateRef } = useCellVizUrlSync(!!dataset, dataset, vizStore);

  // SM overlay URL state sync — uses a dedicated ov= slot so it doesn't
  // collide with the SC viz. The SM dataset is loaded asynchronously by
  // three-scene.tsx when mapping.json has linkColumn "__all__" and pushed
  // into the global SM store, so we read it from there.
  const smDataset = useSingleMoleculeStore((s) => {
    const id = s.currentDatasetId;

    return id ? (s.datasets.get(id) ?? null) : null;
  });
  const smVizStore = useSingleMoleculeVisualizationStore();

  // Resolve the overlay's default genes (Model 3): this SC dataset's saved
  // overlay defaults first, else the SM overlay dataset's own saved defaults,
  // else null (the hook falls back to the pickDefaultGenes heuristic).
  const resolveOverlayDefaultGenes = async (smDs: {
    uniqueGenes: string[];
  }): Promise<string[] | null> => {
    const sc = (registration?.viewerConfig?.defaultGenes ?? []).filter((g) =>
      smDs.uniqueGenes.includes(g),
    );

    if (sc.length > 0) return sc;

    const url = useSingleMoleculeStore.getState().overlaySourceUrl;

    if (url) {
      try {
        const res = await fetch(
          `/api/datasets/by-url?url=${encodeURIComponent(url)}`,
        );

        if (res.ok) {
          const j = await res.json();
          const sm = (
            (j?.viewerConfig?.defaultGenes as string[] | undefined) ?? []
          ).filter((g) => smDs.uniqueGenes.includes(g));

          if (sm.length > 0) return sm;
        }
      } catch {
        // ignore — fall through to the heuristic
      }
    }

    return null;
  };

  useSMOverlayUrlSync(
    !!smDataset,
    smDataset,
    smVizStore,
    resolveOverlayDefaultGenes,
  );

  // Read split params from URL on mount
  useEffect(() => {
    const splitId = searchParams.get("split");
    const splitS3Url = searchParams.get("splitS3Url");
    const splitType = searchParams.get("splitType") as PanelType | null;

    if (splitS3Url && splitType) {
      enableSplit();
      setRightPanelS3(decodeURIComponent(splitS3Url), splitType);
    } else if (splitId && splitType) {
      enableSplit();
      setRightPanel(splitId, splitType);
    }

    if (searchParams.get("sync") === "1") {
      setSyncEnabled(true);
      setSyncFromUrl(true);
    }
  }, []);

  // Write split params to URL when split state changes
  useEffect(() => {
    // Use window.location.search as base to avoid stale Next.js searchParams
    const newParams = new URLSearchParams(window.location.search);

    if (isSplitMode && rightPanelType) {
      if (rightPanelS3Url) {
        newParams.set("splitS3Url", encodeURIComponent(rightPanelS3Url));
        newParams.delete("split");
      } else if (rightPanelDatasetId) {
        newParams.set("split", rightPanelDatasetId);
        newParams.delete("splitS3Url");
      }
      newParams.set("splitType", rightPanelType);
      if (syncEnabled) {
        newParams.set("sync", "1");
      } else {
        newParams.delete("sync");
      }
      router.replace(`?${newParams.toString()}`, { scroll: false });
    } else if (!isSplitMode) {
      newParams.delete("split");
      newParams.delete("splitS3Url");
      newParams.delete("splitType");
      newParams.delete("sync");
      const paramStr = newParams.toString();

      router.replace(paramStr ? `?${paramStr}` : window.location.pathname, {
        scroll: false,
      });
    }
  }, [
    isSplitMode,
    rightPanelDatasetId,
    rightPanelS3Url,
    rightPanelType,
    syncEnabled,
  ]);

  useEffect(() => {
    if (!s3Url) {
      setError("No S3 URL provided");
      setIsLoading(false);

      return;
    }

    loadDatasetFromCustomS3(decodeURIComponent(s3Url));
  }, [s3Url]);

  const loadDatasetFromCustomS3 = async (baseUrl: string) => {
    try {
      setIsLoading(true);
      setError(null);
      setLoadingProgress(0);
      setLoadingMessage("Initializing...");
      defaultsAppliedRef.current = false;
      viewCountedRef.current = false;

      console.log("Loading dataset from custom S3:", baseUrl);

      // Import StandardizedDataset and URL state reader
      const { StandardizedDataset } = await import("@/lib/StandardizedDataset");
      const { tryReadCellVizFromUrl } = await import(
        "@/lib/hooks/useUrlVizSync"
      );

      // Look up the DB registration for this URL (owner-saved defaults,
      // ownership, and the id to count views against). 404 = unregistered.
      let reg: {
        id: string;
        ownerId: string | null;
        adminOwned: boolean;
        viewerConfig: ViewerConfig | null;
      } | null = null;

      try {
        const r = await fetch(
          `/api/datasets/by-url?url=${encodeURIComponent(baseUrl)}`,
        );

        if (r.ok) {
          const j = await r.json();

          if (j.registered) {
            reg = {
              id: j.id,
              ownerId: j.ownerId ?? null,
              adminOwned: !!j.adminOwned,
              viewerConfig: (j.viewerConfig as ViewerConfig | null) ?? null,
            };
          }
        }
      } catch {
        // by-url lookup is best-effort; viewer still loads unregistered.
      }
      setRegistration(reg);
      setViewerRegistration({
        dbId: reg?.id ?? null,
        ownerId: reg?.ownerId ?? null,
        adminOwned: reg?.adminOwned ?? false,
        registered: !!reg,
        viewerConfig: reg?.viewerConfig ?? null,
        s3Url: baseUrl,
      });

      // Read URL state to get priority column hint (bypass auto-detection).
      // Precedence: URL state > owner-saved priority column > heuristic.
      const urlState = tryReadCellVizFromUrl("left");
      const priorityColumnHint =
        urlState?.c || reg?.viewerConfig?.priorityColumn || undefined;

      // Load dataset using fromCustomS3 method
      const standardizedDataset = await StandardizedDataset.fromCustomS3(
        baseUrl,
        (progress, message) => {
          console.log(`${progress}%: ${message}`);
          setLoadingProgress(progress);
          setLoadingMessage(message);
        },
        priorityColumnHint,
      );

      console.log(
        "StandardizedDataset created from custom S3:",
        standardizedDataset,
      );

      // Store dataset in both local state and global store
      setDataset(standardizedDataset);
      addDataset(standardizedDataset);
      console.log("Dataset added to datasetStore");

      setIsLoading(false);
    } catch (err) {
      console.error("Error loading dataset from custom S3:", err);
      setError(
        err instanceof Error
          ? err.message
          : "Failed to load dataset from custom S3",
      );
      setIsLoading(false);
    }
  };

  // Apply defaults when the dataset changes, unless a shared-link URL state
  // already set the view. Precedence: URL state > owner-saved config > heuristic.
  useEffect(() => {
    if (!dataset || hasUrlStateRef.current || defaultsAppliedRef.current)
      return;
    defaultsAppliedRef.current = true;

    if (registration?.viewerConfig) {
      applyViewerConfig(registration.viewerConfig, vizStore, dataset);

      // Per-sample align keys off the transform column's per-cell values, which
      // are lazy-loaded on S3 datasets. applyViewerConfig only sets the column
      // name, so fetch its data (if absent) and bump clusterVersion — the
      // scene's transform effect then re-applies the saved sample transforms.
      const tc = registration.viewerConfig.transformColumn;

      if (tc) {
        loadClusterColumn(dataset, tc)
          .then((fetched) => {
            if (fetched) vizStore.incrementClusterVersion();
          })
          .catch(() => {});
      }
    } else {
      vizStore.setSelectedColumn(selectBestClusterColumn(dataset));
    }
  }, [dataset, registration]);

  // Count one view per load for registered datasets (server dedups per
  // session/day; the ref guards React strict-mode double-invoke).
  useEffect(() => {
    if (!dataset || !registration?.id || viewCountedRef.current) return;
    viewCountedRef.current = true;
    fetch(`/api/datasets/${registration.id}/view`, { method: "POST" }).catch(
      () => {},
    );
  }, [dataset, registration]);

  // Clear the shared registration when leaving the viewer.
  useEffect(() => () => resetViewerRegistration(), []);

  // Loading state
  if (isLoading) {
    return (
      <>
        <div className="fixed inset-0 w-full h-full z-0">
          <LightRays
            lightSpread={1.0}
            mouseInfluence={0.1}
            pulsating={false}
            rayLength={10}
            raysColor="#667eea"
            raysOrigin="top-left"
            raysSpeed={0.8}
          />
        </div>
        <div className="fixed inset-0 w-full h-full z-0">
          <LightRays
            lightSpread={1.0}
            mouseInfluence={0.1}
            pulsating={false}
            rayLength={10}
            raysColor="#764ba2"
            raysOrigin="top-right"
            raysSpeed={0.8}
          />
        </div>
        <div className="relative z-10 flex items-center justify-center h-full">
          <div className="flex flex-col items-center gap-4 w-full max-w-md px-4">
            <Spinner color="secondary" size="lg" />
            <p className={subtitle()}>Loading dataset from S3...</p>
            <Progress
              aria-label="Loading progress"
              className="w-full"
              color="secondary"
              size="md"
              value={loadingProgress}
            />
            <p className="text-sm text-default-500">{loadingMessage}</p>
          </div>
        </div>
      </>
    );
  }

  // Error state
  if (error) {
    return (
      <>
        <div className="fixed inset-0 w-full h-full z-0">
          <LightRays
            lightSpread={1.0}
            mouseInfluence={0.1}
            pulsating={false}
            rayLength={10}
            raysColor="#FF72E1"
            raysOrigin="top-center"
            raysSpeed={0.8}
          />
        </div>
        <div className="relative z-10 flex items-center justify-center h-full p-8">
          <div className="flex flex-col items-center gap-6 max-w-2xl w-full">
            <div className="text-center">
              <h2 className={title({ size: "md", color: "pink" })}>
                Failed to load dataset from S3
              </h2>
              <p className={subtitle({ class: "mt-4" })}>{error}</p>
            </div>
            <div className="flex gap-4">
              <Button color="secondary" onPress={() => router.push("/viewer")}>
                Go to Home
              </Button>
              <Button
                color="default"
                variant="bordered"
                onPress={() =>
                  s3Url && loadDatasetFromCustomS3(decodeURIComponent(s3Url))
                }
              >
                Retry
              </Button>
            </div>
          </div>
        </div>
      </>
    );
  }

  // Dataset loaded
  if (!dataset) {
    return null;
  }

  return (
    <SplitScreenContainer>
      <ClaimDatasetBanner />
      <VisualizationControls />
      <ThreeScene dataset={dataset} />
      <UMAPPanel />
    </SplitScreenContainer>
  );
}

export default function ViewerFromS3Page() {
  return (
    <Suspense
      fallback={
        <div className="flex items-center justify-center h-full">
          <Spinner size="lg" />
        </div>
      }
    >
      <ViewerFromS3Content />
    </Suspense>
  );
}
