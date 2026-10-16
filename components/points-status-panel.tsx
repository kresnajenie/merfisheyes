"use client";

import type { PointsLayer, PointsStatus } from "@/lib/spatialdata/PointsLayer";

import { useEffect, useState } from "react";

import { usePanelSingleMoleculeVisualizationStore } from "@/lib/hooks/usePanelStores";

const count = (n: number) =>
  n >= 1e9
    ? `${(n / 1e9).toFixed(2)} B`
    : n >= 1e6
      ? `${(n / 1e6).toFixed(1)} M`
      : n >= 1e3
        ? `${(n / 1e3).toFixed(1)} k`
        : `${Math.round(n)}`;

const bytes = (n: number) =>
  n >= 1e9 ? `${(n / 1e9).toFixed(2)} GB` : `${(n / 1e6).toFixed(1)} MB`;

/**
 * What the molecule layer is showing and doing: per gene, the level drawn
 * (every molecule, or a uniform sample of 1 in N) and how many molecules are
 * on screen; above, whether it is fetching, waiting to swap in a new view, or
 * up to date; below, what has been read from the server and from the cache.
 */
export function PointsStatusPanel({ layer }: { layer: PointsLayer }) {
  const [status, setStatus] = useState<PointsStatus>(() => layer.status());
  const selected = usePanelSingleMoleculeVisualizationStore(
    (s) => s.selectedGenes,
  );

  useEffect(() => {
    const id = setInterval(() => setStatus(layer.status()), 250);

    return () => clearInterval(id);
  }, [layer]);

  if (status.genes.length === 0) return null;
  const { network } = status;
  const fetching = status.loading > 0 || network.inFlight > 0;
  const waiting = status.genes.some((g) => !g.ready);
  const state = fetching
    ? {
        label: `Fetching · ${network.inFlight} read${network.inFlight === 1 ? "" : "s"} in flight`,
        dot: "bg-amber-400 animate-pulse",
      }
    : waiting
      ? { label: "Drawing the new view", dot: "bg-sky-400 animate-pulse" }
      : { label: "Up to date", dot: "bg-emerald-400" };
  const drawn = status.genes.reduce((a, g) => a + g.drawn, 0);

  return (
    <div
      data-ui-overlay
      className="glass-panel absolute bottom-16 right-6 z-[var(--z-legends)] w-80 rounded-xl p-3 text-xs text-white"
      data-testid="points-status"
    >
      <div className="mb-2 flex items-center gap-2">
        <span className={`h-2 w-2 shrink-0 rounded-full ${state.dot}`} />
        <span className="font-medium">{state.label}</span>
        <span className="ml-auto text-white/60">{count(drawn)} drawn</span>
      </div>
      <ul className="flex flex-col gap-1.5">
        {status.genes.map((g) => (
          <li key={g.gene} className="flex items-center gap-2">
            <span
              className="h-2.5 w-2.5 shrink-0 rounded-full"
              style={{ background: selected.get(g.gene)?.color ?? "#fff" }}
            />
            <span className="w-16 truncate font-medium">{g.gene}</span>
            <span
              className={`rounded px-1.5 py-0.5 ${
                g.fraction >= 1
                  ? "bg-emerald-500/25 text-emerald-200"
                  : "bg-white/10 text-white/80"
              }`}
              title={`Level ${g.level} of ${g.levels - 1}`}
            >
              {g.fraction >= 1
                ? "every molecule"
                : `1 in ${Math.round(1 / g.fraction).toLocaleString()}`}
            </span>
            <span className="ml-auto tabular-nums text-white/70">
              {count(g.drawn)} / {count(g.total)}
            </span>
            {g.wantedLevel !== g.level && (
              <span
                className="text-sky-300"
                title={`Fetching level ${g.wantedLevel}`}
              >
                → L{g.wantedLevel}
              </span>
            )}
          </li>
        ))}
      </ul>
      <div className="mt-2 border-t border-white/10 pt-2 text-white/60">
        Fetched {bytes(network.bytes)} in {network.requests} reads ·{" "}
        {network.cacheHits} from cache
      </div>
    </div>
  );
}
