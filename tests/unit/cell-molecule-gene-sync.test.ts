/**
 * Split-screen gene sync between a panel with cells and one with molecules
 * only (lib/hooks/useSyncCellMoleculeGenes.ts): the cells' one gene and the
 * molecules' genes follow each other by name. The decisions are pure; the
 * binder is driven on real store instances, as the split's panels hold them.
 */
import { describe, expect, it, vi } from "vitest";

import {
  bindCellMoleculeGenes,
  cellGeneToMolecules,
  matchGene,
  moleculeGenesToCell,
} from "@/lib/hooks/useSyncCellMoleculeGenes";
import { createVisualizationStoreInstance } from "@/lib/stores/createVisualizationStore";
import { createDatasetStoreInstance } from "@/lib/stores/createDatasetStore";
import { createSingleMoleculeStoreInstance } from "@/lib/stores/createSingleMoleculeStore";
import { createSingleMoleculeVisualizationStoreInstance } from "@/lib/stores/createSingleMoleculeVisualizationStore";

describe("matchGene", () => {
  it("matches by name, ignoring case only when nothing matches exactly", () => {
    expect(matchGene("Gad1", ["Gad1", "GAD1"])).toBe("Gad1");
    expect(matchGene("GAD1", ["Gad1", "Th"])).toBe("Gad1");
    expect(matchGene("Drd1", ["Gad1", "Th"])).toBeNull();
  });
});

describe("the cells' gene, over to the molecules", () => {
  const molecules = (legend: string[]) => ({
    legend: new Set(legend),
    genes: ["Gad1", "Th", "Drd1"],
  });

  it("adds the gene, keeping the molecules' others", () => {
    expect(cellGeneToMolecules(null, "Th", molecules(["Gad1"]))).toEqual({
      add: "Th",
      remove: null,
      missing: null,
    });
  });

  it("swaps the one it showed for the new one", () => {
    expect(cellGeneToMolecules("Th", "Drd1", molecules(["Gad1", "Th"]))).toEqual({
      add: "Drd1",
      remove: "Th",
      missing: null,
    });
  });

  it("takes it away when the cells' gene is cleared", () => {
    expect(cellGeneToMolecules("Th", null, molecules(["Th"]))).toEqual({
      add: null,
      remove: "Th",
      missing: null,
    });
  });

  it("says which gene the molecules haven't got, and still drops the old one", () => {
    expect(cellGeneToMolecules("Th", "Snap25", molecules(["Th"]))).toEqual({
      add: null,
      remove: "Th",
      missing: "Snap25",
    });
  });

  it("does nothing when the gene is already there", () => {
    expect(cellGeneToMolecules(null, "Gad1", molecules(["Gad1"]))).toEqual({
      add: null,
      remove: null,
      missing: null,
    });
  });
});

describe("the molecules' genes, over to the cells", () => {
  const cells = (gene: string | null) => ({
    gene,
    genes: ["Gad1", "Th", "Drd1"],
  });

  it("shows the gene added last", () => {
    expect(
      moleculeGenesToCell(new Set(["Gad1"]), new Set(["Gad1", "Th", "Drd1"]), cells(null)),
    ).toEqual({ gene: "Drd1", missing: null });
  });

  it("clears the cells' gene when the molecules drop it", () => {
    expect(
      moleculeGenesToCell(new Set(["Gad1", "Th"]), new Set(["Gad1"]), cells("Th")),
    ).toEqual({ gene: null, missing: null });
  });

  it("leaves the cells' gene when another gene is dropped", () => {
    expect(
      moleculeGenesToCell(new Set(["Gad1", "Th"]), new Set(["Th"]), cells("Th")),
    ).toEqual({ gene: undefined, missing: null });
  });

  it("says which gene the cells haven't got, and leaves theirs", () => {
    expect(
      moleculeGenesToCell(new Set(), new Set(["Snap25"]), cells("Th")),
    ).toEqual({ gene: undefined, missing: "Snap25" });
  });
});

/** A cell panel and a molecule panel, each with its dataset loaded. */
function panels() {
  const cellViz = createVisualizationStoreInstance();
  const cells = createDatasetStoreInstance();
  const moleculeViz = createSingleMoleculeVisualizationStoreInstance();
  const molecules = createSingleMoleculeStoreInstance();

  cells.getState().addDataset({
    id: "c",
    name: "cells",
    type: "h5ad",
    clusters: [],
    genes: ["Gad1", "Th", "Drd1", "Snap25"],
    getSummary: () => ({}),
  } as never);
  molecules.getState().addDataset({
    id: "m",
    name: "molecules",
    type: "spatialdata_points",
    uniqueGenes: ["Gad1", "Th", "Drd1"],
    getSummary: () => ({}),
  } as never);

  return { cellViz, cells, moleculeViz, molecules };
}

const legend = (p: ReturnType<typeof panels>) => [
  ...p.moleculeViz.getState().selectedGenesLegend,
];

describe("bindCellMoleculeGenes", () => {
  it("carries a gene chosen on either side to the other, both ways", () => {
    const p = panels();
    const unbind = bindCellMoleculeGenes(p);

    p.cellViz.getState().setSelectedGene("Th");
    expect(legend(p)).toEqual(["Th"]);

    p.moleculeViz.getState().addGene("Drd1");
    expect(p.cellViz.getState().selectedGene).toBe("Drd1");
    // Adding a gene to the molecules never takes one of theirs away
    expect(legend(p)).toEqual(["Th", "Drd1"]);

    // A change made on the cells' side swaps the gene they showed
    p.cellViz.getState().setSelectedGene("Gad1");
    expect(legend(p)).toEqual(["Th", "Gad1"]);

    p.moleculeViz.getState().removeGene("Gad1");
    expect(p.cellViz.getState().selectedGene).toBeNull();
    unbind();
  });

  it("keeps the molecules' own genes when the cells' gene changes", () => {
    const p = panels();

    p.moleculeViz.getState().addGene("Gad1");
    const unbind = bindCellMoleculeGenes(p);

    p.cellViz.getState().setSelectedGene("Th");
    p.cellViz.getState().setSelectedGene("Drd1");
    expect(legend(p)).toEqual(["Gad1", "Drd1"]);
    unbind();
  });

  it("warns of a gene the other side hasn't got and changes nothing there", () => {
    const p = panels();
    const warn = vi.fn();
    const unbind = bindCellMoleculeGenes(p, { warn });

    p.cellViz.getState().setSelectedGene("Snap25");
    expect(warn).toHaveBeenCalledWith("Snap25");
    expect(legend(p)).toEqual([]);
    unbind();
  });

  it("pushes the left panel's genes on a manual toggle, either kind", () => {
    const fromCells = panels();

    fromCells.cellViz.getState().setSelectedGene("Gad1");
    bindCellMoleculeGenes(fromCells, { push: "cells" })();
    expect(legend(fromCells)).toEqual(["Gad1"]);

    const fromMolecules = panels();

    fromMolecules.moleculeViz.getState().addGene("Th");
    bindCellMoleculeGenes(fromMolecules, { push: "molecules" })();
    expect(fromMolecules.cellViz.getState().selectedGene).toBe("Th");
  });

  it("holds off while a shared link restores both panels", () => {
    const p = panels();
    let settling = true;
    const unbind = bindCellMoleculeGenes(p, { settling: () => settling });

    // Each panel restores its own link state: nothing crosses over
    p.cellViz.getState().setSelectedGene("Th");
    p.moleculeViz.getState().addGene("Gad1");
    expect(legend(p)).toEqual(["Gad1"]);
    expect(p.cellViz.getState().selectedGene).toBe("Th");

    settling = false;
    p.cellViz.getState().setSelectedGene("Drd1");
    // Th was never carried over, so only Drd1 is added
    expect(legend(p)).toEqual(["Gad1", "Drd1"]);
    unbind();
  });

  it("stops when unbound", () => {
    const p = panels();

    bindCellMoleculeGenes(p)();
    p.cellViz.getState().setSelectedGene("Th");
    expect(legend(p)).toEqual([]);
  });
});
