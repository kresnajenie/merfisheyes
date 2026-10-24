/**
 * lib/spatialdata/coordinates.ts on the roots of real SpatialData-DB stores
 * (tests/data/sddb/roots, one per technology): the target system each store
 * gets, and where its image lands in the cells' coordinates.
 *
 * Before this module the viewer multiplied every entry of the shapes'
 * transform list, so the Visium full image was drawn 148× too large (198× on
 * Visium HD): docs/plans/SDDB.md §1.
 */
import type { SpatialDataStore } from "@/lib/spatialdata/store";

import { readFileSync } from "node:fs";
import path from "node:path";

import { describe, expect, it } from "vitest";

import {
  applyMap,
  chain,
  chooseImage,
  coordinateSystems,
  elementTransforms,
  invert,
  linearScale,
  readCoordinateSystems,
  compose,
  type ElementTransforms,
  type Map2,
} from "@/lib/spatialdata/coordinates";
import { defaultTableName, elementNames } from "@/lib/spatialdata/store";

const ROOTS = path.join(__dirname, "..", "data", "sddb", "roots");

/** A store known by its root alone, as one opened by URL is. */
function rootStore(uid: string): SpatialDataStore {
  const root = JSON.parse(readFileSync(path.join(ROOTS, `${uid}.json`), "utf8"));
  const nodes = Object.keys(root.consolidated_metadata.metadata);
  const none = { get: async () => undefined } as unknown;

  return {
    root,
    keys: ["zarr.json", ...nodes.map((p) => `${p}/zarr.json`)],
    store: none as SpatialDataStore["store"],
    sidecar: none as SpatialDataStore["sidecar"],
  };
}

/** Level 0 of an image: its pixel → scene map. */
async function level0(sd: SpatialDataStore, name: string): Promise<Map2> {
  const cs = await readCoordinateSystems(sd);
  const placed = cs.place(`images/${name}`);

  if (!placed.ok) throw new Error(placed.reason);
  const attrs = sd.root.consolidated_metadata.metadata[`images/${name}`].attributes;
  const ms = (attrs.ome ?? attrs).multiscales[0];
  const axes = ms.axes.map((a: { name: string }) => a.name);

  return compose(chain(ms.datasets[0].coordinateTransformations, axes), placed.map);
}

async function chosenImage(sd: SpatialDataStore): Promise<string | null> {
  const chosen = chooseImage(
    await readCoordinateSystems(sd),
    elementNames(sd.keys, "images"),
  );

  return chosen.name;
}

describe("SpatialData coordinate systems (SpatialData-DB roots)", () => {
  it("Visium: targets the spots' system and draws the full image at 1 px per unit", async () => {
    const sd = rootStore("fJSOO03VzCTYtD1o");
    const sample = "Targeted_Visium_Human_Cerebellum_Neuroscience";
    const cs = await readCoordinateSystems(sd);

    expect(cs.target).toBe(sample);
    expect(cs.rule).toBe("region identity");
    expect(cs.region).toBe(`shapes/${sample}`);
    expect(await chosenImage(sd)).toBe(`${sample}_full_image`);
    const m = await level0(sd, `${sample}_full_image`);

    // Was 148.1 when the transform list was multiplied
    expect(m.scale[0]).toBeCloseTo(1, 9);
    expect(m.scale[1]).toBeCloseTo(1, 9);
    expect(m.offset).toEqual([0, 0]);
    // The hires image reaches the target through its scale
    const hires = await level0(sd, `${sample}_hires_image`);

    expect(hires.scale[0]).toBeCloseTo(6.666, 3);
  });

  it("Visium HD: the full image at 1 px per unit (was 197.8)", async () => {
    const sd = rootStore("haOzFbJ7Bc9x7NW6");
    const sample = "Visium_HD_3prime_Human_Ovarian_Cancer_FF";
    const cs = await readCoordinateSystems(sd);

    expect(cs.target).toBe(sample);
    // The curated default table's region, not the first table by name
    expect(cs.region).toBe(`shapes/${sample}_square_008um`);
    expect(await chosenImage(sd)).toBe(`${sample}_full_image`);
    expect((await level0(sd, `${sample}_full_image`)).scale[0]).toBeCloseTo(1, 9);
  });

  it("Atera: global, the image placed by undoing the cells' scale", async () => {
    const sd = rootStore("jCLxdY1tWKmbfWOB");
    const cs = await readCoordinateSystems(sd);

    expect([cs.target, cs.rule]).toEqual(["global", "global"]);
    expect(await chosenImage(sd)).toBe("morphology_focus");
    // As before this module: 1 / 4.70588
    expect((await level0(sd, "morphology_focus")).scale[0]).toBeCloseTo(0.2125, 6);
  });

  it("Xenium and spatial transcriptomics: global", async () => {
    for (const uid of ["0W4MBvLWKZXauWhZ", "gZwZDQTAM8vvm8sE"]) {
      const cs = await readCoordinateSystems(rootStore(uid));

      expect([cs.target, cs.rule], uid).toEqual(["global", "global"]);
    }
  });

  it("the points and the other shapes sit in the cells' frame", async () => {
    for (const uid of ["0W4MBvLWKZXauWhZ", "jCLxdY1tWKmbfWOB", "haOzFbJ7Bc9x7NW6"]) {
      const sd = rootStore(uid);
      const cs = await readCoordinateSystems(sd);

      for (const kind of ["shapes", "points"]) {
        for (const name of elementNames(sd.keys, kind)) {
          const p = cs.place(`${kind}/${name}`);

          expect(p.ok && p.map, `${uid} ${kind}/${name}`).toEqual({
            scale: [1, 1],
            offset: [0, 0],
          });
        }
      }
    }
  });
});

describe("coordinateSystems rules", () => {
  const el = (outputs: [string, Record<string, unknown>][]): ElementTransforms => ({
    axes: ["x", "y"],
    transforms: outputs.map(([output, t]) => ({ ...(t as object), output: { name: output } }) as never),
  });
  const id = { type: "identity" };

  it("picks the full image; the hires one only when there is no full image", () => {
    const spots = el([["s", id], ["s_hires", { type: "scale", scale: [0.25, 0.25] }]]);
    const image = (out: [string, Record<string, unknown>][]) => ({
      ...el(out),
      axes: ["c", "y", "x"],
    });
    const full = image([["s", id]]);
    const hires = image([["s_hires", id], ["s", { type: "scale", scale: [1, 4, 4] }]]);
    const withFull = coordinateSystems(
      new Map([["shapes/s", spots], ["images/full", full], ["images/hires", hires]]),
      "shapes/s",
    );

    expect(chooseImage(withFull, ["full", "hires"])).toEqual({ name: "full" });
    const withoutFull = coordinateSystems(
      new Map([["shapes/s", spots], ["images/hires", hires]]),
      "shapes/s",
    );

    expect(chooseImage(withoutFull, ["hires"])).toEqual({ name: "hires" });
    const p = withoutFull.place("images/hires");

    expect(p.ok && p.map.scale).toEqual([4, 4]);
  });

  it("places an affine image, rotation and all; rejects an element with no way into the target", () => {
    // A MALDI optical image as SpatialData-DB writes it: one row per output
    // axis (c, y, x), one column per input axis plus the translation; the
    // x, y part scales by 2 and turns by 90°
    const affine = {
      type: "affine",
      affine: [
        [1, 0, 0, 0],
        [0, 0, 2, 5], // y' = 2 x + 5
        [0, -2, 0, 7], // x' = -2 y + 7
      ],
      input: { name: "cyx", axes: [{ name: "c" }, { name: "y" }, { name: "x" }] },
      output: { name: "s", axes: [{ name: "c" }, { name: "y" }, { name: "x" }] },
    };
    const cs = coordinateSystems(
      new Map([
        ["shapes/s", el([["s", id]])],
        ["images/rotated", { ...el([["s", affine]]), axes: ["c", "y", "x"] }],
        ["images/elsewhere", { ...el([["other", id]]), axes: ["c", "y", "x"] }],
      ]),
      "shapes/s",
    );

    const rotated = cs.place("images/rotated");
    const elsewhere = cs.place("images/elsewhere");

    expect(rotated.ok).toBe(true);
    if (rotated.ok) {
      expect(rotated.map).toEqual({ scale: [0, 0], offset: [7, 5], skew: [-2, 2] });
      expect(applyMap(rotated.map, 1, 0)).toEqual([7, 7]);
      expect(applyMap(rotated.map, 0, 1)).toEqual([5, 5]);
      expect(linearScale(rotated.map)).toBeCloseTo(2);
      // Back again
      const back = invert(rotated.map);

      expect(applyMap(back, 7, 7).map((v) => +v.toFixed(9))).toEqual([1, 0]);
    }
    expect(!elsewhere.ok && elsewhere.reason).toMatch(/no transform into "s"/);
    expect(chooseImage(cs, ["elsewhere", "rotated"]).name).toBe("rotated");
  });

  it("MALDI (human kidney): the optical image placed by an affine with rotation", async () => {
    const sd = rootStore("fF2M7WUbkFnMMhpC");
    const cs = await readCoordinateSystems(sd);
    const p = cs.place("images/optical_image");

    expect(cs.target).toBe("global");
    expect(cs.region).toBe("shapes/maldi_spots");
    expect(p.ok).toBe(true);
    if (p.ok) {
      expect(linearScale(p.map)).toBeCloseTo(15.4258, 3);
      expect(p.map.skew).toBeDefined();
    }
  });

  it("CosMx: the cells' frame is the polygons', which carry a translation into global", async () => {
    const sd = rootStore("571yTEm4hqlkXRPf");
    const cs = await readCoordinateSystems(sd);

    expect(cs.region).toBe("shapes/cells_polygons");
    expect(cs.place("shapes/cells_polygons")).toEqual({ ok: true, map: { scale: [1, 1], offset: [0, 0] } });
    expect(cs.place("shapes/cells_polygons_for_vitessce")).toEqual({ ok: true, map: { scale: [1, 1], offset: [0, 0] } });
  });

  it("reads a plain (unnamed) transform list as a chain", () => {
    const image = elementTransforms("images", {
      ome: {
        multiscales: [
          {
            axes: [{ name: "c" }, { name: "y" }, { name: "x" }],
            coordinateTransformations: [
              { type: "scale", scale: [1, 2, 2] },
              { type: "translation", translation: [0, 10, 20] },
            ],
          },
        ],
      },
    });
    const cs = coordinateSystems(new Map([["images/i", image]]), null);
    const p = cs.place("images/i");

    expect(p.ok && p.map).toEqual({ scale: [2, 2], offset: [20, 10] });
  });
});

describe("defaultTableName", () => {
  it("prefers the curated default, then `table`, then the first by name", () => {
    const curated = { metadata_general: { "Default Table": "square_008um" } };

    expect(defaultTableName(curated, ["square_002um", "square_008um"])).toBe("square_008um");
    expect(defaultTableName(curated, ["square_002um", "table"])).toBe("table");
    expect(defaultTableName({}, ["b", "a"])).toBe("a");
    expect(defaultTableName({}, [])).toBeNull();
  });
});
