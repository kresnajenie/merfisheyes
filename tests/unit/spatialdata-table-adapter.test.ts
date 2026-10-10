import { readFileSync, readdirSync } from "node:fs";
import { readFile } from "node:fs/promises";
import path from "node:path";

import { describe, expect, it } from "vitest";

import { SpatialDataTableAdapter } from "@/lib/spatialdata/SpatialDataTableAdapter";

// Fixture: python scripts/testdata/make-tiny-zarr.py
const DIR = path.resolve(__dirname, "../data/tiny/zarr");
const ROOT = path.join(DIR, "spatialdata_v3.zarr");
const expected = JSON.parse(
  readFileSync(path.join(DIR, "expected.json"), "utf8"),
);

async function openAdapter() {
  const keys = (readdirSync(ROOT, { recursive: true }) as string[]).map((k) =>
    k.split(path.sep).join("/"),
  );
  const store = {
    get: (key: string) =>
      readFile(path.join(ROOT, key)).then(
        (b) => new Uint8Array(b),
        () => undefined,
      ),
  };
  const adapter = new SpatialDataTableAdapter(store as any, keys);

  await adapter.initialize();

  return adapter;
}

describe("SpatialDataTableAdapter", () => {
  it("finds the table and reads genes from the CSC layer of a CSR X", async () => {
    const adapter = await openAdapter();

    expect(adapter.tablePrefix).toBe("tables/table/");
    expect(adapter.getDatasetInfo()).toMatchObject({
      numCells: 30,
      numGenes: 8,
      xFormat: "csc",
    });
  });

  it("reads genes, coordinates and expression", async () => {
    const adapter = await openAdapter();

    expect(await adapter.loadGenes()).toEqual(expected.genes);

    const spatial = await adapter.loadSpatialCoordinates();

    expect(spatial.dimensions).toBe(2);
    expect(Array.from(spatial.coordinates)).toEqual(expected.spatial);

    for (let g = 0; g < expected.genes.length; g++) {
      expect(await adapter.fetchGeneExpression(expected.genes[g])).toEqual(
        expected.expression[g],
      );
    }
    expect(await adapter.fetchGeneExpression("not-a-gene")).toBeNull();
  });

  it("reads categorical and numerical obs columns", async () => {
    const adapter = await openAdapter();
    const [celltype, area] = await adapter.loadClusters(["celltype", "area"]);

    expect(celltype.type).toBe("categorical");
    expect(celltype.uniqueValues).toEqual(["astro", "microglia", "neuron"]);
    expect(
      Array.from(celltype.valueIndices!, (i) => celltype.uniqueValues![i]),
    ).toEqual(expected.celltype);
    expect(area.type).toBe("numerical");
  });

  it("hides the SpatialData join keys from the cluster columns", async () => {
    const adapter = await openAdapter();

    expect(adapter.getClusterColumnInfo().names).toEqual(["area", "celltype"]);
  });
});
