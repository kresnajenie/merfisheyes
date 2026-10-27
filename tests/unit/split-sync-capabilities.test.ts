import { describe, expect, it } from "vitest";

import type { PanelContent } from "@/lib/stores/splitScreenStore";

import {
  canSync,
  hasCells,
  hasMolecules,
  syncCapabilities,
} from "@/lib/stores/splitScreenStore";
import { createDatasetStoreInstance } from "@/lib/stores/createDatasetStore";
import { createSingleMoleculeStoreInstance } from "@/lib/stores/createSingleMoleculeStore";

/**
 * Split-screen sync is offered for what the two panels have open, never for
 * which viewer each is. The button shows exactly when some sync hook will run
 * (the container hands each hook its kind and the button canSync of all), so
 * the two can no longer disagree, as they did when the button compared the
 * viewers: a cell folder beside a SpatialData store showed no button, and two
 * stores showed one that synced nothing.
 */

const none = { cells: false, molecules: false, labelled: false };

/** Each kind of panel, by what its stores hold once it has loaded. */
const PANELS: Record<string, PanelContent> = {
  "cell folder (/viewer)": { ...none, cells: true },
  "cell folder with linked molecules": { ...none, cells: true, molecules: true },
  "store with cells": { ...none, cells: true },
  "store with cells and points": { ...none, cells: true, molecules: true },
  "store with points only": { ...none, molecules: true },
  "molecule folder (/sm-viewer)": { ...none, molecules: true },
  "image (/image-viewer, or a store with only one)": none,
  "labelled molecules (/lm-viewer)": { ...none, labelled: true },
};

const names = Object.keys(PANELS);
const pairs = names.flatMap((l) => names.map((r) => [l, r] as const));

describe("split-screen sync follows what each panel has open", () => {
  it.each(pairs)("%s | %s", (left, right) => {
    const l = PANELS[left];
    const r = PANELS[right];
    const caps = syncCapabilities(l, r);

    // Each hook runs exactly when both panels have its kind
    expect(caps).toEqual({
      cells: l.cells && r.cells,
      molecules: l.molecules && r.molecules,
      labelled: l.labelled && r.labelled,
    });
    // The button shows exactly when some hook runs
    expect(canSync(caps)).toBe(caps.cells || caps.molecules || caps.labelled);
    // Order doesn't matter
    expect(syncCapabilities(r, l)).toEqual(caps);
  });

  it.each([
    // Hidden before: the viewers differed though both have cells
    ["cell folder (/viewer)", "store with cells", "cells"],
    ["store with cells", "cell folder (/viewer)", "cells"],
    // Shown before but synced nothing: neither panel was "cell" or "sm"
    ["store with cells", "store with cells", "cells"],
    ["store with points only", "store with points only", "molecules"],
    // Hidden before: molecules in a store beside a molecule folder
    ["store with points only", "molecule folder (/sm-viewer)", "molecules"],
    ["cell folder with linked molecules", "store with cells and points", "molecules"],
  ] as const)("%s beside %s syncs %s", (left, right, kind) => {
    const caps = syncCapabilities(PANELS[left], PANELS[right]);

    expect(caps[kind]).toBe(true);
    expect(canSync(caps)).toBe(true);
  });

  it("offers nothing between panels with nothing in common", () => {
    const image = PANELS["image (/image-viewer, or a store with only one)"];

    for (const other of names) {
      expect(canSync(syncCapabilities(image, PANELS[other]))).toBe(false);
    }
    expect(
      canSync(
        syncCapabilities(
          PANELS["store with cells"],
          PANELS["molecule folder (/sm-viewer)"],
        ),
      ),
    ).toBe(false);
  });
});

describe("a panel's content is read from its stores as the sync hooks read them", () => {
  // Shaped as the hooks tell them apart: cells carry `clusters`
  const cells = {
    id: "c",
    name: "c",
    type: "h5ad",
    clusters: [],
    getSummary: () => ({}),
  };
  const molecules = {
    id: "m",
    name: "m",
    type: "spatialdata_points",
    uniqueGenes: ["A"],
    getSummary: () => ({}),
  };

  it("has nothing until a dataset is added", () => {
    expect(hasCells(createDatasetStoreInstance().getState())).toBe(false);
    expect(hasMolecules(createSingleMoleculeStoreInstance().getState())).toBe(
      false,
    );
  });

  it("has cells once a cell dataset is current", () => {
    const store = createDatasetStoreInstance();

    store.getState().addDataset(cells as any);
    expect(hasCells(store.getState())).toBe(true);
  });

  it("does not count a molecule dataset in the cell store as cells", () => {
    const store = createDatasetStoreInstance();

    store.getState().addDataset(molecules as any);
    expect(hasCells(store.getState())).toBe(false);
  });

  it("has molecules once a molecule dataset is current", () => {
    const store = createSingleMoleculeStoreInstance();

    store.getState().addDataset(molecules as any);
    expect(hasMolecules(store.getState())).toBe(true);
  });
});
