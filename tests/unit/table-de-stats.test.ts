import { readFileSync, readdirSync } from "node:fs";
import { readFile } from "node:fs/promises";
import path from "node:path";

import { describe, expect, it } from "vitest";

import { SpatialDataTableAdapter } from "@/lib/spatialdata/SpatialDataTableAdapter";

// Fixtures: python scripts/testdata/make-tiny-de.py
const DIR = path.resolve(__dirname, "../data/tiny/zarr");
const expected = JSON.parse(
  readFileSync(path.join(DIR, "de_expected.json"), "utf8"),
);

async function openAdapter(store: string) {
  const root = path.join(DIR, store);
  const keys = (readdirSync(root, { recursive: true }) as string[]).map((k) =>
    k.split(path.sep).join("/"),
  );
  const files = {
    get: (key: string) =>
      readFile(path.join(root, key)).then(
        (b) => new Uint8Array(b),
        () => undefined,
      ),
  };
  const adapter = new SpatialDataTableAdapter(files as any, keys);

  await adapter.initialize();

  return adapter;
}

describe("DE stats kept in the table", () => {
  it("are offered for the columns that have them", async () => {
    const adapter = await openAdapter("spatialdata_v3_de.zarr");

    expect(adapter.getAvailableDeStatsColumns()).toEqual(["celltype"]);
  });

  it("read back as computed by hand", async () => {
    const adapter = await openAdapter("spatialdata_v3_de.zarr");
    const genes = await adapter.loadGenes();
    const stats = await adapter.loadDeStats("celltype", genes);

    expect(stats).not.toBeNull();
    expect(stats!.celltypes).toEqual(expected.celltypes);
    expect(stats!.cellCounts).toEqual(expected.cellCounts);
    expect(stats!.genes).toEqual(genes);
    for (let i = 0; i < expected.means.length; i++) {
      expect(stats!.means[i]).toBeCloseTo(expected.means[i], 5);
      expect(stats!.pctExpressing[i]).toBeCloseTo(expected.pctExpressing[i], 6);
    }
  });

  it("are not offered for a table without them", async () => {
    const adapter = await openAdapter("spatialdata_v3.zarr");

    expect(adapter.getAvailableDeStatsColumns()).toEqual([]);
  });
});
