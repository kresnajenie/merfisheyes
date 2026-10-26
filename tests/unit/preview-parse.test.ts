import fs from "node:fs";
import path from "node:path";

import { beforeAll, describe, expect, it } from "vitest";

/**
 * The browser preview's parse path (the worker's parseH5ad, then
 * StandardizedDataset) on the files it used to misread
 * (docs/plans/PREVIEW-SAVE.md, step 0; fixtures from
 * scripts/testdata/make-tiny-h5ad-edge.py):
 * - anndata 0.13's nullable encodings: 0 genes, columns silently skipped
 * - a var index named in var's `_index` attribute: 0 genes
 * - cells missing a position: dropped, shifting later positions onto the
 *   wrong cells
 */

// The worker reads Files with FileReader, which Node has not
beforeAll(() => {
  if (typeof (globalThis as any).FileReader !== "undefined") return;
  (globalThis as any).FileReader = class {
    result: unknown = null;
    onload: ((e: unknown) => void) | null = null;
    onerror: ((e: unknown) => void) | null = null;
    readAsArrayBuffer(blob: Blob) {
      blob.arrayBuffer().then(
        (b) => {
          this.result = b;
          this.onload?.({ target: this });
        },
        (e) => this.onerror?.({ target: { error: e } }),
      );
    }
    // Papa Parse reads CSVs as text
    readAsText(blob: Blob) {
      blob.text().then(
        (t) => {
          this.result = t;
          this.onload?.({ target: this });
        },
        (e) => this.onerror?.({ target: { error: e } }),
      );
    }
  };
});

const H5AD = path.join(__dirname, "..", "data", "tiny", "h5ad");

async function preview(name: string) {
  const { workerApi } = await import(
    "@/lib/workers/standardized-dataset.worker"
  );
  const { StandardizedDataset } = await import("@/lib/StandardizedDataset");
  const file = new File([fs.readFileSync(path.join(H5AD, name))], name);

  return StandardizedDataset.fromSerializedData(
    await workerApi.parseH5ad(file),
  );
}

/** Cell i's position, whichever shape the coordinates come in. */
function position(ds: any, i: number): number[] {
  const { coordinates: c, dimensions: d } = ds.spatial;

  return Array.isArray(c) ? c[i] : Array.from(c.slice(i * d, i * d + d));
}

function positionCount(ds: any): number {
  const { coordinates: c, dimensions: d } = ds.spatial;

  return Array.isArray(c) ? c.length : c.length / d;
}

describe("the preview reads anndata 0.13's nullable encodings", () => {
  it("lists the genes of a nullable-string var index", async () => {
    const ds = await preview("nullable.h5ad");

    expect(ds.genes).toEqual(["Gad1", "Sst", "Pvalb", "Vip"]);
  });

  it("offers nullable integer and boolean columns", async () => {
    const ds = await preview("nullable.h5ad");
    const columns = (ds.clusters ?? []).map((c: any) => c.column);

    expect(columns).toEqual(
      expect.arrayContaining(["celltype", "count", "flag"]),
    );
    const flag = ds.clusters!.find((c: any) => c.column === "flag") as any;

    // pandas' labels, and "" for the missing one, as a missing category has
    expect(flag.uniqueValues).toEqual(
      expect.arrayContaining(["True", "False", ""]),
    );
  });
});

describe("the preview takes the genes from the index var names", () => {
  it("reads var/gene_ids when var's _index attribute says so", async () => {
    const ds = await preview("named_index.h5ad");

    // The same list the store saved from this file reads
    // (SpatialDataTableAdapter.loadGenes)
    expect(ds.genes).toEqual(
      [0, 1, 2, 3].map((i) => `ENSMUSG${String(i).padStart(11, "0")}`),
    );
  });
});

describe("the preview keeps cells that have no position", () => {
  it("keeps one position per cell, in order", async () => {
    const ds = await preview("missing_xy.h5ad");

    expect(positionCount(ds)).toBe(24);
    expect(position(ds, 3).every(Number.isNaN)).toBe(true);
    expect(position(ds, 7).every(Number.isNaN)).toBe(true);
    // Cell 8 is still at its own position, not cell 10's
    expect(position(ds, 8)).toEqual([80, 40]);
    expect(position(ds, 23)).toEqual([230, 115]);
  });

  it("does the same for a MERSCOPE region", async () => {
    const { workerApi } = await import(
      "@/lib/workers/standardized-dataset.worker"
    );
    const { StandardizedDataset } = await import("@/lib/StandardizedDataset");
    const dir = path.join(__dirname, "..", "data", "tiny", "merscope", "tiny");
    const read = (name: string) =>
      fs.readFileSync(path.join(dir, name), "utf8");
    // Cells 3 and 7 lose their x (rows 4 and 8, after the header)
    const metadata = read("cell_metadata.csv").split("\n");
    const original = metadata.map((line) => line.split(","));

    for (const row of [4, 8]) {
      const cols = metadata[row].split(",");

      cols[1] = "";
      metadata[row] = cols.join(",");
    }
    const files = (
      [
        ["cell_metadata.csv", metadata.join("\n")],
        ["cell_by_gene.csv", read("cell_by_gene.csv")],
        ["cell_categories.csv", read("cell_categories.csv")],
      ] as const
    ).map(([name, text]) => {
      const f = new File([text], name);

      Object.defineProperty(f, "webkitRelativePath", { value: `tiny/${name}` });

      return f;
    });
    const ds = StandardizedDataset.fromSerializedData(
      await workerApi.parseMerscope(files),
    );

    expect(positionCount(ds)).toBe(300);
    expect(position(ds, 3).every(Number.isNaN)).toBe(true);
    // Cell 8 is at row 9's position, its own
    const [, x8, y8] = original[9].map(Number);

    expect(position(ds, 8)[0]).toBeCloseTo(x8, 2);
    expect(position(ds, 8)[1]).toBeCloseTo(y8, 2);
  });
});
