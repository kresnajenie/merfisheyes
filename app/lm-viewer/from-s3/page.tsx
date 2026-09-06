"use client";

import { Progress, Spinner } from "@heroui/react";
import { Suspense, useCallback, useEffect, useRef, useState } from "react";
import { useSearchParams } from "next/navigation";

import LabelledMoleculeControls from "@/components/labelled-molecule-controls";
import LabelledMoleculeLegends from "@/components/labelled-molecule-legends";
import LabelledMoleculeTopControls from "@/components/labelled-molecule-top-controls";
import { subtitle } from "@/components/primitives";
import LabelledMoleculeThreeScene from "@/components/labelled-molecule-three-scene";
import { StandardizedDataset } from "@/lib/StandardizedDataset";
import { useLmVizUrlSync } from "@/lib/hooks/useLmVizUrlSync";
import { useLabelledMoleculeVisualizationStore } from "@/lib/stores/labelledMoleculeVisualizationStore";
import { loadClusterColumn } from "@/lib/utils/load-cluster-column";

/** Columns the three menus need before the scene can draw. */
const REQUIRED_COLUMNS = ["gene", "domain", "cell"];

function LabelledMoleculeViewer() {
  const searchParams = useSearchParams();
  const baseUrl = searchParams.get("url");

  const [dataset, setDataset] = useState<StandardizedDataset | null>(null);
  const [clusterVersion, setClusterVersion] = useState(0);
  const [progress, setProgress] = useState(0);
  const [message, setMessage] = useState("Loading…");
  const [error, setError] = useState<string | null>(null);
  const loadedFor = useRef<string | null>(null);

  // Restores a shared link once the columns exist, then mirrors state into `v=`.
  const vizStore = useLabelledMoleculeVisualizationStore();

  useLmVizUrlSync(!!dataset, vizStore);

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
    } catch (e) {
      setError(e instanceof Error ? e.message : String(e));
    }
  }, []);

  useEffect(() => {
    if (!baseUrl || loadedFor.current === baseUrl) return;
    loadedFor.current = baseUrl;
    load(baseUrl);
  }, [baseUrl, load]);

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
            size="md"
            value={progress}
          />
          <p className="text-sm text-default-500">{message}</p>
        </div>
      </div>
    );
  }

  return (
    <div className="absolute inset-0 overflow-hidden bg-black">
      <LabelledMoleculeThreeScene
        clusterVersion={clusterVersion}
        dataset={dataset}
      />
      <LabelledMoleculeControls
        clusterVersion={clusterVersion}
        dataset={dataset}
      />
      <LabelledMoleculeLegends
        clusterVersion={clusterVersion}
        dataset={dataset}
      />
      <LabelledMoleculeTopControls />
    </div>
  );
}

export default function Page() {
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
