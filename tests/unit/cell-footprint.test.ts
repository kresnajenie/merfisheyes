/**
 * The cells' real shape and size (lib/spatialdata/footprint.ts): Visium
 * spots are circles of one radius, read from the shapes file's footer
 * alone; Visium HD bins are squares whose side is the grid's pitch; other
 * shapes (polygons, circles of varying radius) leave the cells as dots.
 */
import type { SpatialDataStore } from "@/lib/spatialdata/store";

import { readdirSync, readFileSync, statSync } from "node:fs";
import path from "node:path";

import { describe, expect, it } from "vitest";

import { gridPitch, readCellFootprint } from "@/lib/spatialdata/footprint";

const ZARR = path.resolve(__dirname, "../data/tiny/zarr");

function listFiles(dir: string, prefix = ""): string[] {
  return readdirSync(dir).flatMap((name) => {
    const full = path.join(dir, name);
    const key = prefix ? `${prefix}/${name}` : name;

    return statSync(full).isDirectory() ? listFiles(full, key) : [key];
  });
}

/** A store spatialdata wrote (make-tiny-visium.py), with the reads counted. */
function fixtureStore(name: string) {
  const VISIUM = path.join(ZARR, name);
  const reads: string[] = [];
  const file = (key: string) => {
    try {
      return new Uint8Array(readFileSync(path.join(VISIUM, key.replace(/^\/+/, ""))));
    } catch {
      return undefined;
    }
  };
  const sd = {
    root: undefined,
    keys: listFiles(VISIUM),
    store: {
      get: async (key: string) => {
        reads.push(`${key} (whole)`);

        return file(key);
      },
      getRange: async (key: string, range: { suffixLength?: number; offset?: number; length?: number }) => {
        const bytes = file(key);

        if (!bytes) return undefined;
        if (range.suffixLength !== undefined) {
          reads.push(`${key} (last ${range.suffixLength} bytes)`);

          return bytes.slice(Math.max(0, bytes.length - range.suffixLength));
        }

        return bytes.slice(range.offset!, range.offset! + range.length!);
      },
    },
    sidecar: { get: async () => undefined },
  } as unknown as SpatialDataStore;

  return { sd, reads };
}

describe("readCellFootprint", () => {
  it("draws Visium spots as circles of their radius, from the footer alone", async () => {
    const { sd, reads } = fixtureStore("visium_v3.zarr");
    // make-tiny-visium.py: a 6 × 6 grid of spots of radius 12
    const footprint = await readCellFootprint(sd, new Float32Array([28, 28, 68, 28]), 2);

    expect(footprint).toEqual({ shape: "circle", size: 24 });
    expect(reads.filter((r) => r.includes("shapes.parquet"))).toEqual([
      "/shapes/tiny/shapes.parquet (last 65536 bytes)",
    ]);
  });
});

describe("readCellFootprint on Visium HD bins", () => {
  it("draws them as squares whose side is the grid's pitch", async () => {
    const { sd, reads } = fixtureStore("visium_hd_v3.zarr");
    const root = JSON.parse(readFileSync(path.join(ZARR, "visium_hd_v3.zarr", "zarr.json"), "utf8"));

    (sd as { root: unknown }).root = root;
    // The bins' centres, every 16 units (as the table's obsm/spatial holds them)
    const coords = new Float32Array([8, 8, 24, 8, 56, 8, 8, 24]);

    expect(await readCellFootprint(sd, coords, 2)).toEqual({ shape: "square", size: 16 });
    expect(reads.filter((r) => r.includes("shapes.parquet"))).toEqual([
      "/shapes/tiny_square_016um/shapes.parquet (last 65536 bytes)",
    ]);
  });
});

describe("gridPitch", () => {
  const grid = (side: number, angle: number, skip: (c: number, r: number) => boolean) => {
    const coords: number[] = [];
    const [c, s] = [Math.cos(angle), Math.sin(angle)];

    for (let col = 0; col < 40; col++) {
      for (let row = 0; row < 40; row++) {
        if (skip(col, row)) continue;
        const [x, y] = [col * side, row * side];

        coords.push(19000 + x * c - y * s, 9000 + x * s + y * c);
      }
    }

    return new Float32Array(coords);
  };

  it("is a bin's side on a grid with bins missing", () => {
    const coords = grid(8, 0, (c, r) => (c * 7 + r * 3) % 5 === 0);

    expect(gridPitch(coords, 2)).toBeCloseTo(8, 3);
  });

  it("holds on a grid rotated in the image's pixels (as Visium HD's is)", () => {
    // 0.08° as measured on the ovary store; no two columns share an x
    const coords = grid(29.2225, (0.08 * Math.PI) / 180, () => false);

    expect(gridPitch(coords, 2)).toBeCloseTo(29.2225, 2);
  });

  it("is nothing for a single cell", () => {
    expect(gridPitch(new Float32Array([5, 1]), 2)).toBeNull();
  });
});
