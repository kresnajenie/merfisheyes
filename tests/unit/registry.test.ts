/**
 * The element registry over the tiny stores: what `openElements` finds in a
 * store, the rules it applies (a store's own points stand in for a linked
 * molecule folder; one element failing leaves the rest), and how a molecule
 * dataset is classified. The layers themselves are covered by
 * layer-contract.test.ts.
 */
import type { SpatialDataStore } from "@/lib/spatialdata/store";
import type { StandardizedDataset } from "@/lib/StandardizedDataset";

import { readdirSync, statSync } from "node:fs";
import { readFile } from "node:fs/promises";
import path from "node:path";

import { afterEach, describe, expect, it, vi } from "vitest";

import { SingleMoleculeDataset } from "@/lib/SingleMoleculeDataset";
import {
  ELEMENTS,
  elementDefinition,
  moleculeElements,
  openElements,
} from "@/lib/scene/registry";
import { pointsElementOf } from "@/lib/spatialdata/pointsMolecules";

const ZARR = path.resolve(__dirname, "../data/tiny/zarr");

function listFiles(dir: string, prefix = ""): string[] {
  return readdirSync(dir).flatMap((name) => {
    const full = path.join(dir, name);
    const key = prefix ? `${prefix}/${name}` : name;

    return statSync(full).isDirectory() ? listFiles(full, key) : [key];
  });
}

/**
 * A store composed the way the e2e fixture server composes one: the table
 * and image of spatialdata_v3.zarr, the points element of points_v1.zarr
 * and the polygons of shapes_v1, as `parts` says.
 */
function fixtureStore(
  parts: { points?: boolean; shapes?: boolean; image?: boolean },
  mapping?: object,
): SpatialDataStore {
  const main = path.join(ZARR, "spatialdata_v3.zarr");
  const points = path.join(ZARR, "points_v1.zarr");
  const shapes = path.join(ZARR, "shapes_v1");
  const keys = [
    ...listFiles(main).filter((k) => parts.image || !k.startsWith("images/")),
    ...(parts.points ? listFiles(points) : []),
    ...(parts.shapes ? listFiles(shapes).map((k) => `shapes/${k}`) : []),
  ];
  const file = (key: string) => {
    const k = key.replace(/^\/+/, "");

    if (k.startsWith("points/")) return path.join(points, k);
    if (k.startsWith("shapes/")) return path.join(shapes, k.slice(7));

    return path.join(main, k);
  };
  const read = async (key: string) =>
    readFile(file(key)).then(
      (b) => new Uint8Array(b),
      () => undefined,
    );

  return {
    keys,
    store: { get: read },
    sidecar: {
      get: async (key) =>
        key === "/mapping.json" && mapping
          ? new TextEncoder().encode(JSON.stringify(mapping))
          : undefined,
    },
  };
}

const cells = {
  adapter: { readInstanceIds: async () => ["0", "1"] },
} as unknown as StandardizedDataset;

afterEach(() => {
  vi.restoreAllMocks();
});

describe("ELEMENTS", () => {
  it("names every kind once, in draw order, each with an opener and a layer", () => {
    expect(ELEMENTS.map((d) => d.kind)).toEqual([
      "image",
      "shapes",
      "points",
      "molecules",
    ]);
    for (const def of ELEMENTS) {
      expect(typeof def.detect).toBe("function");
      expect(typeof def.open).toBe("function");
      expect(typeof def.createLayer).toBe("function");
      expect(elementDefinition(def.kind)).toBe(def);
    }
  });
});

describe("openElements", () => {
  it("opens the image, the polygons and the points of a store with cells", async () => {
    const store = fixtureStore({ points: true, shapes: true, image: true });
    const opened = await openElements(store, { id: "ds", name: "tiny", cells });

    expect(opened.image?.levels.length).toBeGreaterThan(0);
    expect(typeof opened.shapes).toBe("function");
    expect(opened.points).toBeInstanceOf(SingleMoleculeDataset);
    expect(opened.points!.id).toBe("ds_transcripts");
    expect(pointsElementOf(opened.points!)).not.toBeNull();
    expect(opened.molecules).toBeUndefined();
  });

  it("leaves out what the store does not hold", async () => {
    const opened = await openElements(fixtureStore({}), {
      id: "ds",
      name: "tiny",
      cells,
    });

    expect(opened).toEqual({});
  });

  it("opens no polygons without the cells they belong to", async () => {
    const opened = await openElements(fixtureStore({ shapes: true }), {
      id: "ds",
      name: "tiny",
    });

    expect(opened.shapes).toBeUndefined();
  });

  it("lets a store's own points stand in for a linked molecule folder", async () => {
    const fetch = vi
      .spyOn(globalThis, "fetch")
      .mockRejectedValue(new Error("no network in this test"));
    const store = fixtureStore(
      { points: true },
      { linkColumn: "__all__", links: { __all__: "https://x/molecules" } },
    );
    const opened = await openElements(store, { id: "ds", name: "tiny", cells });

    expect(opened.points).toBeDefined();
    expect(opened.molecules).toBeUndefined();
    expect(fetch).not.toHaveBeenCalled();
  });

  it("logs an element that fails to open and keeps the others", async () => {
    vi.spyOn(globalThis, "fetch").mockRejectedValue(new Error("offline"));
    const error = vi.spyOn(console, "error").mockImplementation(() => {});
    const store = fixtureStore(
      { image: true },
      { linkColumn: "__all__", links: { __all__: "https://x/molecules" } },
    );
    const opened = await openElements(store, { id: "ds", name: "tiny", cells });

    expect(opened.image).toBeDefined();
    expect(opened.molecules).toBeUndefined();
    expect(error).toHaveBeenCalledWith(
      "Error opening molecules:",
      expect.any(Error),
    );
    expect(opened.problems).toEqual([expect.stringMatching(/^molecules: /)]);
  });
});

/** The Visium-shaped store spatialdata wrote (make-tiny-visium.py), minus `drop`. */
function visiumStore(drop: (key: string) => boolean = () => false): SpatialDataStore {
  const dir = path.join(ZARR, "visium_v3.zarr");
  const keys = listFiles(dir).filter((k) => !drop(k));
  const read = async (key: string) => {
    const k = key.replace(/^\/+/, "");

    return drop(k)
      ? undefined
      : readFile(path.join(dir, k)).then(
          (b) => new Uint8Array(b),
          () => undefined,
        );
  };

  return { keys, store: { get: read }, sidecar: { get: async () => undefined } } as never;
}

describe("openElements on a Visium store (no global system)", () => {
  it("draws the full image at one pixel per spot unit, not the downscaled ones", async () => {
    const opened = await openElements(visiumStore(), { id: "v", name: "tiny", cells });

    expect(opened.image?.name).toBe("tiny_full_image");
    expect(opened.image?.levels.map((l) => [l.pixelWidth, l.pixelHeight])).toEqual([
      [1, 1],
      [2, 2],
    ]);
    expect(opened.problems).toBeUndefined();
  });

  it("draws the hires image, scaled into the spots' system, when there is no full image", async () => {
    const opened = await openElements(
      visiumStore((k) => k.startsWith("images/tiny_full_image")),
      { id: "v", name: "tiny", cells },
    );

    expect(opened.image?.name).toBe("tiny_hires_image");
    expect(opened.image?.levels[0].pixelWidth).toBe(4);
  });
});

describe("moleculeElements", () => {
  it("tells a points element from per-gene files, and nothing from null", async () => {
    const store = fixtureStore({ points: true });
    const { points } = await openElements(store, {
      id: "ds",
      name: "tiny",
      cells,
    });
    const legacy = new SingleMoleculeDataset({
      id: "legacy",
      name: "legacy",
      type: "custom",
      uniqueGenes: ["Gad1"],
      geneIndex: new Map(),
      dimensions: 2,
      scalingFactor: 1,
    });

    expect(moleculeElements(points)).toEqual({ points });
    expect(moleculeElements(legacy)).toEqual({ molecules: legacy });
    expect(moleculeElements(null)).toEqual({});
  });
});
