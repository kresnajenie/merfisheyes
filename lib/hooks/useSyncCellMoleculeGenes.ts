"use client";

import type { DatasetStore } from "../stores/createDatasetStore";
import type { SingleMoleculeStore } from "../stores/createSingleMoleculeStore";
import type { SingleMoleculeVisualizationStore } from "../stores/createSingleMoleculeVisualizationStore";
import type { VisualizationStore } from "../stores/createVisualizationStore";

import { toast } from "react-toastify";
import { useEffect } from "react";

import { useDatasetStore } from "../stores/datasetStore";
import { useSingleMoleculeStore } from "../stores/singleMoleculeStore";
import { useSingleMoleculeVisualizationStore } from "../stores/singleMoleculeVisualizationStore";
import { hasCells, useSplitScreenStore } from "../stores/splitScreenStore";
import { useVisualizationStore } from "../stores/visualizationStore";

/**
 * Gene sync between a panel with cells and one with molecules only (the
 * split's "genes" kind, splitScreenStore.syncCapabilities). The cells show
 * one gene's expression, the molecules any number of genes, so:
 *
 * - the cells' gene is one of the molecules' genes: choosing another swaps
 *   it, clearing it takes it away, and the genes the molecules had besides
 *   stay;
 * - a gene added to the molecules becomes the cells' gene, and taking away
 *   the one the cells show clears it there.
 *
 * Genes are matched by name, ignoring case when that is the only match; one
 * the other side hasn't got is left out with a warning, as the other syncs do.
 */

/** A gene's name in `names`: the same, or the same but for case. */
export function matchGene(
  gene: string,
  names: readonly string[],
): string | null {
  if (names.includes(gene)) return gene;
  const lower = gene.toLowerCase();

  return names.find((n) => n.toLowerCase() === lower) ?? null;
}

/** What the molecules do when the cells' gene goes from `prev` to `gene`. */
export function cellGeneToMolecules(
  prev: string | null,
  gene: string | null,
  molecules: { legend: ReadonlySet<string>; genes: readonly string[] },
): { add: string | null; remove: string | null; missing: string | null } {
  const out = { add: null, remove: null, missing: null } as {
    add: string | null;
    remove: string | null;
    missing: string | null;
  };

  if (prev === gene) return out;
  const before = prev ? matchGene(prev, molecules.genes) : null;
  const after = gene ? matchGene(gene, molecules.genes) : null;

  if (gene && !after) out.missing = gene;
  if (before && before !== after && molecules.legend.has(before)) {
    out.remove = before;
  }
  if (after && !molecules.legend.has(after)) out.add = after;

  return out;
}

/**
 * The cells' gene once the molecules' genes go from `prev` to `legend`:
 * undefined to leave it, null to clear it. The last gene added wins.
 */
export function moleculeGenesToCell(
  prev: ReadonlySet<string>,
  legend: ReadonlySet<string>,
  cells: { gene: string | null; genes: readonly string[] },
): { gene: string | null | undefined; missing: string | null } {
  const added = [...legend].filter((g) => !prev.has(g));

  if (added.length > 0) {
    const last = added[added.length - 1];
    const match = matchGene(last, cells.genes);

    if (!match) return { gene: undefined, missing: last };

    return { gene: match === cells.gene ? undefined : match, missing: null };
  }
  const shown = cells.gene;
  const removedShown = [...prev].some(
    (g) => !legend.has(g) && shown && matchGene(g, [shown]) === shown,
  );

  return { gene: removedShown ? null : undefined, missing: null };
}

/** The part of a zustand store the binder needs: state and its changes. */
interface Bound<S> {
  getState(): S;
  subscribe(listener: (state: S) => void): () => void;
}

type CellViz = Bound<
  Pick<
    ReturnType<VisualizationStore["getState"]>,
    "selectedGene" | "setSelectedGene"
  >
>;
type MoleculeViz = Bound<
  Pick<
    ReturnType<SingleMoleculeVisualizationStore["getState"]>,
    "selectedGenesLegend" | "geneDataCache" | "addGene" | "removeGene"
  >
>;
type Datasets = Bound<{
  currentDatasetId: string | null;
  datasets: Map<string, unknown>;
}>;

export interface GeneSyncSides {
  cellViz: CellViz;
  cells: Datasets;
  moleculeViz: MoleculeViz;
  molecules: Datasets;
}

const current = (store: Datasets) => {
  const s = store.getState();

  return s.currentDatasetId ? s.datasets.get(s.currentDatasetId) : undefined;
};
const cellGenes = (store: Datasets): string[] =>
  (current(store) as { genes?: string[] } | undefined)?.genes ?? [];
const moleculeGenes = (store: Datasets): string[] =>
  (current(store) as { uniqueGenes?: string[] } | undefined)?.uniqueGenes ?? [];

/**
 * Keep the two sides' genes together until the returned function is called.
 * `push` first carries `from`'s current gene(s) over, as the other syncs push
 * the left panel on a manual toggle; `settling` holds propagation off (the
 * snapshots still follow) while both panels restore a shared link.
 */
export function bindCellMoleculeGenes(
  sides: GeneSyncSides,
  options: {
    push?: "cells" | "molecules" | null;
    settling?: () => boolean;
    warn?: (gene: string) => void;
  } = {},
): () => void {
  const { cellViz, cells, moleculeViz, molecules } = sides;
  const settling = options.settling ?? (() => false);
  const warn = options.warn ?? (() => {});
  let syncing = false;
  let prevGene = cellViz.getState().selectedGene;
  // The store replaces the set on every change, so a reference is a snapshot
  let prevLegend: ReadonlySet<string> =
    moleculeViz.getState().selectedGenesLegend;

  const toMolecules = (prev: string | null, gene: string | null) => {
    const m = moleculeViz.getState();
    const step = cellGeneToMolecules(prev, gene, {
      legend: m.selectedGenesLegend,
      genes: moleculeGenes(molecules),
    });

    if (step.missing) warn(step.missing);
    if (step.remove) m.removeGene(step.remove);
    if (step.add) moleculeViz.getState().addGene(step.add);
  };

  const toCells = (prev: ReadonlySet<string>, legend: ReadonlySet<string>) => {
    const c = cellViz.getState();
    const step = moleculeGenesToCell(prev, legend, {
      gene: c.selectedGene,
      genes: cellGenes(cells),
    });

    if (step.missing) warn(step.missing);
    if (step.gene !== undefined) c.setSelectedGene(step.gene);
  };

  const guarded = (fn: () => void) => {
    if (syncing) return;
    syncing = true;
    try {
      fn();
    } finally {
      syncing = false;
    }
  };

  if (options.push === "cells") {
    guarded(() => toMolecules(null, cellViz.getState().selectedGene));
  } else if (options.push === "molecules") {
    guarded(() =>
      toCells(new Set(), moleculeViz.getState().selectedGenesLegend),
    );
  }
  prevGene = cellViz.getState().selectedGene;
  prevLegend = moleculeViz.getState().selectedGenesLegend;

  const unsubCells = cellViz.subscribe((s) => {
    const gene = s.selectedGene;

    if (gene === prevGene) return;
    const prev = prevGene;

    prevGene = gene;
    if (!settling()) guarded(() => toMolecules(prev, gene));
  });
  const unsubMolecules = moleculeViz.subscribe((s) => {
    const legend = s.selectedGenesLegend;

    if (legend === prevLegend) return;
    const prev = prevLegend;

    prevLegend = legend;
    if (!settling()) guarded(() => toCells(prev, legend));
  });

  return () => {
    unsubCells();
    unsubMolecules();
  };
}

let lastToastTime = 0;

function throttledWarning(gene: string) {
  const now = Date.now();

  if (now - lastToastTime < 2000) return;
  lastToastTime = now;
  toast.warning(`Gene "${gene}" not found in the other dataset`, {
    autoClose: 3000,
  });
}

/**
 * The split's gene sync between cells and molecules, while sync is on and
 * the two panels are cells and molecules (`available`). Whichever side has
 * the cells, the left panel is pushed first, as the other syncs do.
 */
export function useSyncCellMoleculeGenes(
  right: {
    visualizationStore: VisualizationStore;
    datasetStore: DatasetStore;
    singleMoleculeVisualizationStore: SingleMoleculeVisualizationStore;
    singleMoleculeStore: SingleMoleculeStore;
  },
  available: boolean,
) {
  const syncEnabled = useSplitScreenStore((s) => s.syncEnabled);
  const isActive = syncEnabled && available;

  useEffect(() => {
    if (!isActive) return;
    const fromUrl = useSplitScreenStore.getState().syncFromUrl;
    let settling = fromUrl;
    const cellsOnRight = hasCells(right.datasetStore.getState());
    const left = {
      viz: useVisualizationStore as unknown as CellViz,
      cells: useDatasetStore as unknown as Datasets,
      moleculeViz:
        useSingleMoleculeVisualizationStore as unknown as MoleculeViz,
      molecules: useSingleMoleculeStore as unknown as Datasets,
    };
    const sides: GeneSyncSides = cellsOnRight
      ? {
          cellViz: right.visualizationStore as unknown as CellViz,
          cells: right.datasetStore as unknown as Datasets,
          moleculeViz: left.moleculeViz,
          molecules: left.molecules,
        }
      : {
          cellViz: left.viz,
          cells: left.cells,
          moleculeViz:
            right.singleMoleculeVisualizationStore as unknown as MoleculeViz,
          molecules: right.singleMoleculeStore as unknown as Datasets,
        };
    const unbind = bindCellMoleculeGenes(sides, {
      // The left panel is pushed on a manual toggle, whichever kind it is
      push: fromUrl ? null : cellsOnRight ? "molecules" : "cells",
      settling: () => settling,
      warn: throttledWarning,
    });
    // As the other syncs: a shared link gets 3 s to restore both panels
    const settleTimer = fromUrl
      ? setTimeout(() => {
          settling = false;
          useSplitScreenStore.getState().setSyncFromUrl(false);
        }, 3000)
      : undefined;

    return () => {
      unbind();
      if (settleTimer) clearTimeout(settleTimer);
    };
  }, [
    isActive,
    right.visualizationStore,
    right.datasetStore,
    right.singleMoleculeVisualizationStore,
    right.singleMoleculeStore,
  ]);
}
