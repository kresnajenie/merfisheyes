import type { SpatialDataStore } from "@/lib/spatialdata/store";

import { readFileSync, readdirSync, statSync } from "node:fs";
import { open, readFile } from "node:fs/promises";
import path from "node:path";

import { describe, expect, it } from "vitest";

import {
  findPointsElements,
  openPointsElement,
} from "@/lib/spatialdata/points";

// Fixture: python scripts/testdata/make-tiny-points.py
const DIR = path.resolve(__dirname, "../data/tiny/zarr");
const ROOT = path.join(DIR, "points_v1.zarr");
const expected = JSON.parse(
  readFileSync(path.join(DIR, "points_expected.json"), "utf8"),
);

type Range = { offset: number; length: number } | { suffixLength: number };

function listFiles(dir: string, prefix = ""): string[] {
  return readdirSync(dir).flatMap((name) => {
    const full = path.join(dir, name);
    const key = prefix ? `${prefix}/${name}` : name;

    return statSync(full).isDirectory() ? listFiles(full, key) : [key];
  });
}

/** The fixture as the viewer sees a store: keys, whole and ranged reads. */
function fixtureStore(): SpatialDataStore {
  const file = (key: string) => path.join(ROOT, key.replace(/^\/+/, ""));

  return {
    keys: listFiles(ROOT),
    sidecar: { get: async () => undefined },
    store: {
      get: async (key: string) =>
        readFile(file(key)).then(
          (b) => new Uint8Array(b),
          () => undefined,
        ),
      getRange: async (key: string, range: Range) => {
        const handle = await open(file(key)).catch(() => null);

        if (!handle) return undefined;
        try {
          const { size } = await handle.stat();
          const offset =
            "suffixLength" in range ? size - range.suffixLength : range.offset;
          const length =
            "suffixLength" in range ? range.suffixLength : range.length;
          const { buffer } = await handle.read(
            Buffer.alloc(length),
            0,
            length,
            offset,
          );

          return new Uint8Array(buffer);
        } finally {
          await handle.close();
        }
      },
    },
  };
}

describe("points in Zarr", () => {
  it("finds the element and reads its genes, levels and counts", async () => {
    const store = fixtureStore();

    expect(await findPointsElements(store)).toEqual(["transcripts"]);
    const el = await openPointsElement(store, "transcripts");

    expect(el.features).toEqual(expected.genes);
    expect(el.dims).toBe(3);
    expect(el.levels.map((l) => l.nPoints)).toEqual(expected.levels);
    expect(Array.from(await el.featureCounts())).toEqual(expected.counts);
  });

  it("reads a gene's every point at level 0, across chunks and shards", async () => {
    const el = await openPointsElement(fixtureStore(), "transcripts");
    const ft = await el.featureTiles(0, 1);
    const total = ft.counts.reduce((a, b) => a + b, 0);

    expect(total).toBe(expected.counts[1]);
    // Tiles ascend and their runs are back to back
    for (let i = 1; i < ft.tiles.length; i++) {
      expect(ft.tiles[i]).toBeGreaterThan(ft.tiles[i - 1]);
      expect(ft.starts[i]).toBe(ft.starts[i - 1] + ft.counts[i - 1]);
    }
    const xyz = await el.readXyz(0, ft.starts[0], ft.starts[0] + total);
    const got: number[][] = [];

    for (let i = 0; i < total; i++) {
      got.push([0, 1, 2].map((k) => Math.round(xyz[i * 3 + k] * 1000) / 1000));
    }
    got.sort((a, b) => a[0] - b[0] || a[1] - b[1] || a[2] - b[2]);
    expect(got).toEqual(expected.beta_xyz_sorted);
  });

  it("finds a gene's points in a region from the tiles alone", async () => {
    const el = await openPointsElement(fixtureStore(), "transcripts");
    const { origin, tileSize, nx } = el.levels[0].grid;
    const [x0, y0, x1, y1] = expected.alpha_bbox;
    const ft = await el.featureTiles(0, 0);
    let inside = 0;

    for (let i = 0; i < ft.tiles.length; i++) {
      const tx = ft.tiles[i] % nx;
      const ty = Math.floor(ft.tiles[i] / nx);

      // Only tiles overlapping the region are read
      if (
        (tx + 1) * tileSize + origin[0] <= x0 ||
        tx * tileSize + origin[0] >= x1 ||
        (ty + 1) * tileSize + origin[1] <= y0 ||
        ty * tileSize + origin[1] >= y1
      ) {
        continue;
      }
      const xyz = await el.readXyz(0, ft.starts[i], ft.starts[i] + ft.counts[i]);

      for (let j = 0; j < ft.counts[i]; j++) {
        const x = xyz[j * 3];
        const y = xyz[j * 3 + 1];

        if (x >= x0 && x < x1 && y >= y0 && y < y1) inside++;
      }
    }
    expect(inside).toBe(expected.alpha_in_bbox);
  });

  it("levels shrink, each a sample of the one below", async () => {
    const el = await openPointsElement(fixtureStore(), "transcripts");

    for (let l = 1; l < el.levels.length; l++) {
      expect(el.levels[l].nPoints).toBeLessThan(el.levels[l - 1].nPoints);
      expect(el.levels[l].grid.tileSize).toBe(el.levels[l - 1].grid.tileSize * 2);
    }
  });
});
