import { readFileSync, readdirSync } from "node:fs";
import { open, readFile } from "node:fs/promises";
import path from "node:path";
import { gzipSync } from "node:zlib";

import { describe, expect, it } from "vitest";

import { SpatialDataTableAdapter } from "@/lib/spatialdata/SpatialDataTableAdapter";
import { encodeDeStatsBuffer } from "@/lib/utils/de-stats";

// Fixtures: python scripts/testdata/make-tiny-zarr.py
const DIR = path.resolve(__dirname, "../data/tiny/zarr");
const expected = JSON.parse(
  readFileSync(path.join(DIR, "expected.json"), "utf8"),
);

type Range = { offset: number; length: number } | { suffixLength: number };

/** Part of a file, as a store reads a shard. */
async function readRange(file: string, range: Range) {
  const handle = await open(file).catch(() => null);

  if (!handle) return undefined;
  try {
    const { size } = await handle.stat();
    const offset =
      "suffixLength" in range ? size - range.suffixLength : range.offset;
    const length = "suffixLength" in range ? range.suffixLength : range.length;
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
}

// The same table written two ways; every read must give the same answers.
const STORES = [
  { name: "plain arrays and strings", store: "spatialdata_v3.zarr", extra: [] },
  {
    name: "sharded arrays and nullable strings",
    store: "spatialdata_v3_sharded.zarr",
    extra: ["note"],
  },
  // python scripts/testdata/make-tiny-csc-x.py: X itself CSC, no layer
  { name: "CSC X", store: "spatialdata_v3_csc_x.zarr", extra: [] },
];

describe.each(STORES)(
  "SpatialDataTableAdapter ($name)",
  ({ store: dir, extra }) => {
    const ROOT = path.join(DIR, dir);
    const columns = ["area", "celltype", ...extra];

    async function openAdapter(sidecarFiles?: Record<string, Uint8Array>) {
      const keys = (readdirSync(ROOT, { recursive: true }) as string[]).map(
        (k) => k.split(path.sep).join("/"),
      );
      const store = {
        get: (key: string) =>
          readFile(path.join(ROOT, key)).then(
            (b) => new Uint8Array(b),
            () => undefined,
          ),
        getRange: (key: string, range: Range) =>
          readRange(path.join(ROOT, key), range),
      };
      const sidecar = sidecarFiles && {
        get: async (key: string) => sidecarFiles[key],
      };
      const adapter = new SpatialDataTableAdapter(store as any, keys, sidecar);

      await adapter.initialize();

      return adapter;
    }

    it("finds the table and reads genes from CSC (X, or a layer of a CSR X)", async () => {
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

      expect(adapter.getClusterColumnInfo().names).toEqual(columns);
    });

    it.runIf(extra.includes("note"))(
      "reads a nullable string column, missing values as empty",
      async () => {
        const adapter = await openAdapter();
        const [note] = await adapter.loadClusters(["note"]);

        expect(
          Array.from(note.valueIndices!, (i) => note.uniqueValues![i]),
        ).toEqual(expected.note);
      },
    );

    it("adds columns and DE stats from annotations/ beside the store", async () => {
      const categories = ["0", "1"];
      const codes = Uint16Array.from({ length: 30 }, (_, i) =>
        i < 10 ? 1 : 0,
      );
      const stats = {
        column: "leiden",
        celltypes: categories,
        cellCounts: [20, 10],
        genes: expected.genes,
        means: Float32Array.from({ length: 16 }, (_, i) => i / 2),
        pctExpressing: Float32Array.from({ length: 16 }, (_, i) => i / 16),
      };
      const adapter = await openAdapter({
        "/annotations/annotations.json": new TextEncoder().encode(
          JSON.stringify({
            columns: { leiden: { categories } },
            de_stats: ["leiden"],
          }),
        ),
        "/annotations/leiden.bin": new Uint8Array(codes.buffer),
        "/annotations/de/leiden.bin.gz": gzipSync(encodeDeStatsBuffer(stats)),
      });

      expect(adapter.getClusterColumnInfo().names).toContain("leiden");

      const [leiden] = await adapter.loadClusters(["leiden"]);

      expect(leiden.uniqueValues).toEqual(categories);
      expect(Array.from(leiden.valueIndices!)).toEqual(Array.from(codes));

      expect(adapter.getAvailableDeStatsColumns()).toEqual(["leiden"]);
      expect(await adapter.loadDeStats("leiden", expected.genes)).toEqual(
        stats,
      );
      expect(await adapter.loadDeStats("celltype", expected.genes)).toBeNull();
    });

    it("has no annotations without the folder", async () => {
      const adapter = await openAdapter({});

      expect(adapter.getClusterColumnInfo().names).toEqual(columns);
      expect(adapter.getAvailableDeStatsColumns()).toEqual([]);
    });
  },
);
