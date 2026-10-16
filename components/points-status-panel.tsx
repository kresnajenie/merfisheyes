"use client";

import { useEffect, useState } from "react";

import {
  DETAIL_BUDGET,
  type PointsDetail,
  type PointsLayer,
  type PointsStatus,
} from "@/lib/spatialdata/PointsLayer";
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

const DETAILS: { key: PointsDetail; label: string; hint: string }[] = [
  { key: "auto", label: "Auto", hint: "Up to 4 M molecules on screen" },
  { key: "more", label: "More", hint: "Up to 16 M molecules on screen" },
  {
    key: "all",
    label: "All",
    hint: "Every molecule in view, up to 64 M — needs a strong GPU",
  },
];
const DETAIL_KEY = "merfisheyes.points.detail";

function Spinner() {
  return (
    <span className="h-3 w-3 shrink-0 animate-spin rounded-full border-2 border-white/25 border-t-white" />
  );
}

/**
 * What the molecule layer is showing and doing: per gene, the level drawn
 * (every molecule, or a uniform sample of 1 in N), molecules on screen, and
 * whether its view is still loading; a progress bar for the view as a whole;
 * how much to draw (Auto, More, All); and what has been read so far.
 */
export function PointsStatusPanel({ layer }: { layer: PointsLayer }) {
  const [status, setStatus] = useState<PointsStatus>(() => layer.status());
  const selected = usePanelSingleMoleculeVisualizationStore(
    (s) => s.selectedGenes,
  );

  useEffect(() => {
    try {
      const saved = localStorage.getItem(DETAIL_KEY) as PointsDetail | null;

      if (saved && saved in { auto: 1, more: 1, all: 1 }) {
        layer.setDetail(saved);
      }
    } catch {
      // No storage: the default stays
    }
    const id = setInterval(() => setStatus(layer.status()), 200);

    return () => clearInterval(id);
  }, [layer]);

  const choose = (detail: PointsDetail) => {
    layer.setDetail(detail);
    setStatus(layer.status());
    try {
      localStorage.setItem(DETAIL_KEY, detail);
    } catch {
      // Not remembered, still applied
    }
  };

  if (status.genes.length === 0) return null;
  const wanted = status.genes.reduce((a, g) => a + g.wanted, 0);
  const loaded = status.genes.reduce((a, g) => a + g.loaded, 0);
  const drawn = status.genes.reduce((a, g) => a + g.drawn, 0);
  const allInView = status.genes.reduce((a, g) => a + g.allInView, 0);
  const loading = loaded < wanted;
  const pct = wanted ? Math.floor((loaded / wanted) * 100) : 100;

  return (
    <div
      data-ui-overlay
      className="glass-panel absolute bottom-16 right-6 z-[var(--z-legends)] w-96 rounded-xl p-3 text-xs text-white"
      data-testid="points-status"
    >
      <div className="flex items-center gap-2">
        {loading ? (
          <Spinner />
        ) : (
          <span className="h-2 w-2 shrink-0 rounded-full bg-emerald-400" />
        )}
        <span className="font-medium">
          {loading ? `Loading this view · ${pct}%` : "All loaded for this view"}
        </span>
        <span className="ml-auto tabular-nums text-white/60">
          {count(drawn)} drawn
        </span>
      </div>
      <div className="mt-2 h-1 overflow-hidden rounded-full bg-white/10">
        <div
          className={`h-full rounded-full transition-[width] duration-200 ${
            loading ? "bg-sky-400" : "bg-emerald-400"
          }`}
          style={{ width: `${pct}%` }}
        />
      </div>
      <div className="mt-1 flex text-[11px] tabular-nums text-white/50">
        <span>
          {count(loaded)} of {count(wanted)} molecules
        </span>
        <span className="ml-auto">{bytes(status.network.bytes)} fetched</span>
      </div>

      <ul className="mt-3 flex flex-col gap-1.5">
        {status.genes.map((g) => {
          const genePct = g.wanted
            ? Math.floor((g.loaded / g.wanted) * 100)
            : 100;

          return (
            <li key={g.gene} className="flex items-center gap-2">
              <span
                className="h-2.5 w-2.5 shrink-0 rounded-full"
                style={{ background: selected.get(g.gene)?.color ?? "#fff" }}
              />
              <span className="w-16 truncate font-medium">{g.gene}</span>
              <span
                className={`whitespace-nowrap rounded px-1.5 py-0.5 ${
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
              <span className="ml-auto whitespace-nowrap tabular-nums text-white/70">
                {count(g.drawn)} / {count(g.total)}
              </span>
              {g.loaded < g.wanted ? (
                <span
                  className="flex w-12 items-center justify-end gap-1 text-sky-300"
                  title="Loading"
                >
                  <Spinner />
                  {genePct}%
                </span>
              ) : (
                <span
                  className="w-12 text-right text-emerald-300"
                  title="Loaded"
                >
                  ✓
                </span>
              )}
            </li>
          );
        })}
      </ul>

      <div className="mt-3 border-t border-white/10 pt-2">
        <div className="flex items-center gap-2">
          <span className="text-white/60">Detail</span>
          <div className="ml-auto flex overflow-hidden rounded-md border border-white/15">
            {DETAILS.map((d) => (
              <button
                key={d.key}
                className={`px-2 py-0.5 ${
                  status.detail === d.key
                    ? "bg-white/20 text-white"
                    : "text-white/60 hover:bg-white/10"
                }`}
                title={d.hint}
                type="button"
                onClick={() => choose(d.key)}
              >
                {d.label}
              </button>
            ))}
          </div>
        </div>
        <div className="mt-1 text-[11px] text-white/50">
          Every molecule in this view: {count(allInView)}
          {allInView > DETAIL_BUDGET.all
            ? " — more than All can draw; zoom in to see every one"
            : status.detail !== "all" && allInView > DETAIL_BUDGET.more
              ? " — a lot for a GPU"
              : ""}
        </div>
      </div>
    </div>
  );
}
