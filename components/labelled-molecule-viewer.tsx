"use client";

import type { ViewerConfig } from "@/lib/utils/viewer-config";

import { Progress, Spinner } from "@heroui/react";
import { Suspense, useCallback, useEffect, useRef, useState } from "react";
import { useSearchParams } from "next/navigation";

import LabelledMoleculeControls from "@/components/labelled-molecule-controls";
import LabelledMoleculeLegends from "@/components/labelled-molecule-legends";
import LabelledMoleculeTopControls from "@/components/labelled-molecule-top-controls";
import StageRail, { type ProjectDatasetSummary } from "@/components/stage-rail";
import { ClaimDatasetBanner } from "@/components/claim-dataset-banner";
import { subtitle } from "@/components/primitives";
import LabelledMoleculeThreeScene from "@/components/labelled-molecule-three-scene";
import { StandardizedDataset } from "@/lib/StandardizedDataset";
import { useLmVizUrlSync } from "@/lib/hooks/useLmVizUrlSync";
import { labelledMoleculeVisualizationStore } from "@/lib/stores/labelledMoleculeVisualizationStore";
import { useLabelledMoleculeVisualizationStore } from "@/lib/stores/labelledMoleculeVisualizationStore";
import { useViewerRegistrationStore } from "@/lib/stores/viewerRegistrationStore";
import { loadClusterColumn } from "@/lib/utils/load-cluster-column";

/** Columns the three menus need before the scene can draw. */
const REQUIRED_COLUMNS = ["gene", "domain", "cell"];

function LabelledMoleculeViewer() {
  const searchParams = useSearchParams();
  const baseUrl = searchParams.get("url");

  // The dataset actually being shown. Seeded from `?url=`, then swapped in
  // place by the stage rail — `useSearchParams` deliberately does not drive
  // this, because the viz-state writer rewrites the query string constantly
  // and a swap must not be triggered by one of those writes.
  const [activeUrl, setActiveUrl] = useState<string | null>(baseUrl);
  const [dataset, setDataset] = useState<StandardizedDataset | null>(null);
  /** True while a second dataset loads over an already-drawn one. */
  const [swapping, setSwapping] = useState(false);
  const [clusterVersion, setClusterVersion] = useState(0);
  const [progress, setProgress] = useState(0);
  const [message, setMessage] = useState("Loading…");
  const [error, setError] = useState<string | null>(null);
  // The project this dataset belongs to, if any — drives the stage rail.
  const [projectId, setProjectId] = useState<string | null>(null);
  const loadedFor = useRef<string | null>(null);

  // Restores a shared link once the columns exist, then mirrors state into `v=`.
  const vizStore = useLabelledMoleculeVisualizationStore();

  useLmVizUrlSync(!!dataset, vizStore);

  // Back/forward across swapped datasets.
  useEffect(() => {
    const onPop = () => {
      const u = new URLSearchParams(window.location.search).get("url");

      if (u) setActiveUrl(u);
    };

    window.addEventListener("popstate", onPop);

    return () => window.removeEventListener("popstate", onPop);
  }, []);

  const load = useCallback(async (url: string) => {
    try {
      setError(null);
      // Colour defaults to cell, so ask for that column up front rather than
      // letting the generic heuristic pick `gene`.
      const ds = await StandardizedDataset.fromCustomS3(
        url,
        (p, m) => {
          setProgress(p);
          setMessage(m);
        },
        "cell",
      );

      setMessage("Loading label columns…");
      // Concurrently: each is its own worker round trip, and serialising them
      // added two avoidable trips to every load.
      await Promise.all(REQUIRED_COLUMNS.map((c) => loadClusterColumn(ds, c)));

      setDataset(ds);
      setClusterVersion((v) => v + 1);
      setSwapping(false);

      // Ownership: look the dataset up by its S3 URL so an owner sees their
      // saved defaults and everyone else gets the claim banner. Best effort —
      // an unregistered dataset still opens.
      let config: ViewerConfig | null = null;

      try {
        const res = await fetch(
          `/api/datasets/by-url?url=${encodeURIComponent(url)}`,
        );

        if (res.ok) {
          const j = await res.json();

          config = (j?.viewerConfig as ViewerConfig | null) ?? null;
          setProjectId((j?.projectId as string | null) ?? null);
          useViewerRegistrationStore.getState().set({
            dbId: j?.id ?? null,
            ownerId: j?.ownerId ?? null,
            adminOwned: !!j?.adminOwned,
            registered: !!j?.id,
            viewerConfig: config,
            s3Url: url,
          });
        } else {
          useViewerRegistrationStore
            .getState()
            .set({ registered: false, s3Url: url });
        }
      } catch {
        useViewerRegistrationStore
          .getState()
          .set({ registered: false, s3Url: url });
      }

      // Apply the owner's saved camera unless the link already carries one —
      // an explicit shared view beats the dataset default.
      const st = labelledMoleculeVisualizationStore.getState();

      if (config?.camera && !st.pendingCamera && !st.camera) {
        st.applyCamera(config.camera);
      }
    } catch (e) {
      setError(e instanceof Error ? e.message : String(e));
      setSwapping(false);
    }
  }, []);

  useEffect(() => {
    if (!activeUrl || loadedFor.current === activeUrl) return;
    // A second load keeps the outgoing scene on screen behind a progress
    // overlay; only the first shows a bare spinner.
    if (loadedFor.current !== null) setSwapping(true);
    loadedFor.current = activeUrl;
    setProgress(0);
    setMessage("Loading…");
    load(activeUrl);
  }, [activeUrl, load]);

  /**
   * Open another embryo without a navigation.
   *
   * Nothing carries over. Cell and domain names are disjoint across stages —
   * no value of either column appears in all 45 datasets — so a kept selection
   * would mostly name values the incoming embryo has not got.
   */
  const switchDataset = useCallback(
    (d: ProjectDatasetSummary) => {
      if (d.s3BaseUrl === activeUrl) return;

      labelledMoleculeVisualizationStore.getState().reset();

      const next = new URL(window.location.href);

      next.searchParams.set("url", d.s3BaseUrl);
      // `v=` encodes selections of the outgoing embryo.
      next.searchParams.delete("v");
      window.history.pushState(null, "", next.toString());

      setActiveUrl(d.s3BaseUrl);
    },
    [activeUrl],
  );

  if (!baseUrl) {
    return (
      <div className="flex h-screen items-center justify-center text-sm text-default-400">
        Pass a dataset folder as <code className="mx-1">?url=</code>
      </div>
    );
  }

  if (error) {
    return (
      <div className="flex h-screen flex-col items-center justify-center gap-2">
        <p className="text-white">Failed to load dataset</p>
        <p className="max-w-lg text-center font-mono text-xs text-default-500">
          {error}
        </p>
      </div>
    );
  }

  if (!dataset) {
    return (
      <div className="flex h-screen items-center justify-center">
        <div className="flex w-full max-w-md flex-col items-center gap-4 px-4">
          <Spinner color="secondary" size="lg" />
          <p className={subtitle()}>Loading labelled molecules…</p>
          <Progress
            aria-label="Loading progress"
            className="w-full"
            color="secondary"
            // The coordinate download dominates and now reports real bytes, so
            // the bar moves continuously rather than sitting at 30%.
            size="md"
            value={progress}
          />
          <div className="flex w-full items-baseline justify-between text-sm">
            <span className="text-default-500">{message}</span>
            <span className="tabular-nums text-default-400">{progress}%</span>
          </div>
        </div>
      </div>
    );
  }

  return (
    <div className="absolute inset-0 overflow-hidden bg-black">
      {/* Keyed on the dataset: a swap tears the WebGL scene and every panel's
          local state down rather than relying on each effect to notice. */}
      <LabelledMoleculeThreeScene
        key={activeUrl}
        clusterVersion={clusterVersion}
        dataset={dataset}
        hasStageRail={!!projectId}
      />
      <LabelledMoleculeControls
        key={`c-${activeUrl}`}
        clusterVersion={clusterVersion}
        dataset={dataset}
        hasStageRail={!!projectId}
      />
      <LabelledMoleculeLegends
        key={`l-${activeUrl}`}
        clusterVersion={clusterVersion}
        dataset={dataset}
      />
      <LabelledMoleculeTopControls />
      <ClaimDatasetBanner />
      {swapping && (
        <div className="pointer-events-none absolute inset-x-0 top-0 z-[var(--z-chrome)] flex justify-center p-3">
          <div className="flex w-72 flex-col gap-1 rounded-xl bg-black/70 px-4 py-3 backdrop-blur">
            <Progress
              aria-label="Loading dataset"
              color="secondary"
              size="sm"
              value={progress}
            />
            <div className="flex items-baseline justify-between text-[11px]">
              <span className="text-default-400">{message}</span>
              <span className="tabular-nums text-default-500">{progress}%</span>
            </div>
          </div>
        </div>
      )}
      {projectId && (
        <StageRail
          currentUrl={activeUrl}
          projectId={projectId}
          onOpen={switchDataset}
          onSplit={(d: ProjectDatasetSummary) => {
            const url = new URL(window.location.href);

            url.searchParams.set("splitS3Url", d.s3BaseUrl);
            url.searchParams.set("splitType", "lm");
            window.location.href = url.toString();
          }}
        />
      )}
    </div>
  );
}

export default function LabelledMoleculeViewerPage() {
  return (
    <Suspense
      fallback={
        <div className="flex h-screen items-center justify-center bg-black" />
      }
    >
      <LabelledMoleculeViewer />
    </Suspense>
  );
}
