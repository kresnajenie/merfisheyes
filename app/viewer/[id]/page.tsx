"use client";

import type { StandardizedDataset } from "@/lib/StandardizedDataset";
import type { PanelType } from "@/lib/stores/splitScreenStore";

import { Suspense, useEffect, useRef, useState } from "react";
import { useParams, useRouter, useSearchParams } from "next/navigation";
import { Button, Progress, Spinner } from "@heroui/react";

import { ThreeScene } from "@/components/three-scene";
import { VisualizationControls } from "@/components/visualization-controls";
import UMAPPanel from "@/components/umap-panel";
import { SplitScreenContainer } from "@/components/split-screen-container";
import { useVisualizationStore } from "@/lib/stores/visualizationStore";
import { useDatasetStore } from "@/lib/stores/datasetStore";
import { useSplitScreenStore } from "@/lib/stores/splitScreenStore";
import { useViewerRegistrationStore } from "@/lib/stores/viewerRegistrationStore";
import { type ViewerConfig } from "@/lib/utils/viewer-config";
import { applyCellOpenState } from "@/lib/viewer/open-dataset";
import {
  useCellVizUrlSync,
  tryReadCellVizFromUrl,
  useSMOverlayUrlSync,
} from "@/lib/hooks/useUrlVizSync";
import { useSingleMoleculeStore } from "@/lib/stores/singleMoleculeStore";
import { useSingleMoleculeVisualizationStore } from "@/lib/stores/singleMoleculeVisualizationStore";
import LightRays from "@/components/react-bits/LightRays";
import { subtitle, title } from "@/components/primitives";

function ViewerByIdContent() {
  const params = useParams();
  const router = useRouter();
  const searchParams = useSearchParams();
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
  const [viewerConfig, setViewerConfig] = useState<ViewerConfig | null>(null);
  const setViewerRegistration = useViewerRegistrationStore((s) => s.set);
  const resetViewerRegistration = useViewerRegistrationStore((s) => s.reset);
  const viewCountedRef = useRef(false);

  const datasetId = params.id as string;

  // URL visualization state sync
  useCellVizUrlSync(!!dataset, dataset, vizStore);

  // SM overlay URL state sync — uses a dedicated ov= slot so the overlay's
  // selected molecule genes persist in the link. The SM dataset is loaded
  // asynchronously by three-scene.tsx when mapping.json has linkColumn
  // "__all__" and pushed into the global SM store, so we read it from there.
  const smOverlayDataset = useSingleMoleculeStore((s) => {
    const id = s.currentDatasetId;

    return id ? (s.datasets.get(id) ?? null) : null;
  });
  const smVizStore = useSingleMoleculeVisualizationStore();

  // Resolve the overlay's default genes: this SC dataset's saved overlay
  // defaults first, else the SM overlay dataset's own saved defaults, else
  // null (the hook falls back to the pickDefaultGenes heuristic).
  const resolveOverlayDefaultGenes = async (smDs: {
    uniqueGenes: string[];
  }): Promise<string[] | null> => {
    const sc = (viewerConfig?.defaultGenes ?? []).filter((g) =>
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
    !!smOverlayDataset,
    smOverlayDataset,
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
    if (!datasetId) {
      setError("No dataset ID provided");
      setIsLoading(false);

      return;
    }

    resolveDataset(datasetId);
  }, [datasetId]);

  const resolveDataset = async (id: string) => {
    // Step 1: Check Zustand store first
    const storeDataset = useDatasetStore.getState().datasets.get(id);

    if (storeDataset && "spatial" in storeDataset) {
      console.log("Dataset found in store:", id);
      setDataset(storeDataset as StandardizedDataset);
      setIsLoading(false);

      return;
    }

    // Step 2: Load from S3
    loadDatasetFromS3(id);
  };

  const loadDatasetFromS3 = async (id: string) => {
    try {
      setIsLoading(true);
      setError(null);
      setLoadingProgress(0);
      setLoadingMessage("Initializing...");

      console.log("Loading dataset from server:", id);

      const { StandardizedDataset } = await import("@/lib/StandardizedDataset");

      viewCountedRef.current = false;

      const urlState = tryReadCellVizFromUrl("left");

      // One pre-flight call so we can route based on storage format.
      // 202 = "still uploading/processing" — surface clearly.
      const metaRes = await fetch(`/api/datasets/${id}`);

      if (metaRes.status === 202) {
        const body = await metaRes.json().catch(() => ({}));

        throw new Error(
          body.message || "Dataset is not ready yet. Try again shortly.",
        );
      }
      if (!metaRes.ok) {
        const body = await metaRes.json().catch(() => ({}));

        throw new Error(
          body.message || `Dataset fetch failed: ${metaRes.status}`,
        );
      }
      const meta = await metaRes.json();
      const config = (meta.viewerConfig as ViewerConfig | null) ?? null;

      setViewerConfig(config);
      setViewerRegistration({
        dbId: id,
        ownerId: meta.ownerId ?? null,
        adminOwned: !!meta.adminOwned,
        registered: true,
        viewerConfig: config,
        s3Url: null,
      });

      // Priority column precedence: URL state > owner-saved config > heuristic.
      const priorityColumn = urlState?.c || config?.priorityColumn || undefined;
      const onProg = (progress: number, message: string) => {
        console.log(`${progress}%: ${message}`);
        setLoadingProgress(progress);
        setLoadingMessage(message);
      };

      const standardizedDataset =
        meta.formatVersion === "zarr"
          ? await StandardizedDataset.fromS3Zarr(id, onProg, priorityColumn)
          : await StandardizedDataset.fromS3(id, onProg, priorityColumn);

      console.log("StandardizedDataset created:", standardizedDataset);

      // Shared open pipeline (reset → owner config incl. per-sample align →
      // heuristic). Runs BEFORE setDataset so the URL sync hook overlays any
      // v= state on top afterwards — same ordering as the split panels.
      // Awaited so the loading screen only clears once the first frame can
      // paint already-aligned coordinates.
      setLoadingMessage("Applying alignment...");
      await applyCellOpenState({
        dataset: standardizedDataset,
        config,
        store: vizStore,
      });

      setDataset(standardizedDataset);
      addDataset(standardizedDataset);
      console.log("Dataset added to datasetStore");

      setIsLoading(false);
    } catch (err) {
      console.error("Error loading dataset:", err);
      setError(err instanceof Error ? err.message : "Failed to load dataset");
      setIsLoading(false);
    }
  };

  // Count one view per load (server dedups per session/day).
  useEffect(() => {
    if (!dataset || viewCountedRef.current) return;
    viewCountedRef.current = true;
    fetch(`/api/datasets/${datasetId}/view`, { method: "POST" }).catch(
      () => {},
    );
  }, [dataset]);

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
            raysColor="#FFD700"
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
            raysColor="#FFD700"
            raysOrigin="top-right"
            raysSpeed={0.8}
          />
        </div>
        <div className="relative z-10 flex items-center justify-center h-full">
          <div className="flex flex-col items-center gap-4 w-full max-w-md px-4">
            <Spinner color="primary" size="lg" />
            <p className={subtitle()}>Loading dataset...</p>
            <Progress
              aria-label="Loading progress"
              className="w-full"
              color="primary"
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
                Failed to load dataset
              </h2>
              <p className={subtitle({ class: "mt-4" })}>{error}</p>
            </div>
            <div className="flex gap-4">
              <Button color="primary" onPress={() => router.push("/viewer")}>
                Go to Home
              </Button>
              <Button
                color="default"
                variant="bordered"
                onPress={() => resolveDataset(datasetId)}
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
      <VisualizationControls />
      <ThreeScene dataset={dataset} />
      <UMAPPanel />
    </SplitScreenContainer>
  );
}

export default function ViewerByIdPage() {
  return (
    <Suspense
      fallback={
        <div className="flex items-center justify-center h-full">
          <Spinner size="lg" />
        </div>
      }
    >
      <ViewerByIdContent />
    </Suspense>
  );
}
