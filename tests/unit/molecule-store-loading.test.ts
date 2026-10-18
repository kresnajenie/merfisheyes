/**
 * A legacy molecule folder read through the one store: the manifest, a
 * gene whole with byte progress, a gene streamed in batches, a count from
 * the gzip trailer, a missing unassigned file, and a cancelled read.
 */
import type { Backend } from "@/lib/storage/ElementStore";

import { gzip } from "pako";
import { describe, expect, it, vi } from "vitest";

import { folderGeneKey, sanitizeGeneName } from "@/lib/molecules/storeLoading";
import { SingleMoleculeDataset } from "@/lib/SingleMoleculeDataset";
import { ElementStore } from "@/lib/storage/ElementStore";

const coords = (n: number, seed = 0) =>
  Float32Array.from({ length: n * 3 }, (_, i) => seed + i / 2);

const manifest = (genes: Record<string, number>, raw = true) =>
  gzip(
    JSON.stringify({
      name: "folder",
      type: "xenium",
      genes: {
        unique_gene_names: Object.keys(genes),
        molecule_counts: Object.fromEntries(
          Object.entries(genes).map(([g, n]) => [g, { assigned: n }]),
        ),
      },
      statistics: {
        spatial_dimensions: 3,
        total_molecules: Object.values(genes).reduce((a, b) => a + b, 0),
        unique_genes: Object.keys(genes).length,
      },
      processing: {
        scaling_factor: 2,
        ...(raw && { coordinate_range: "raw_rounded_2dp" }),
      },
      has_unassigned: true,
    }),
  );

/** A folder as objects, with every read recorded. */
function folder(objects: Record<string, Uint8Array>) {
  const reads: string[] = [];
  const backend: Backend = {
    async get(key) {
      reads.push(`get ${key}`);

      return objects[key];
    },
    async getRange(key, range) {
      reads.push(`range ${key}`);
      const o = objects[key];

      if (!o || !("suffixLength" in range)) return undefined;

      return o.slice(-range.suffixLength);
    },
    async open(key, signal) {
      reads.push(`open ${key}`);
      if (signal?.aborted) throw new DOMException("aborted", "AbortError");
      const o = objects[key];

      return (
        o &&
        new Response(o, { headers: { "content-length": String(o.length) } })
      );
    },
  };

  return { store: new ElementStore(backend), reads };
}

const gzipped = (c: Float32Array) => gzip(new Uint8Array(c.buffer));

describe("a molecule folder through the store", () => {
  it("names gene files as the processor did", () => {
    expect(sanitizeGeneName("Gad-1 (x)")).toBe("Gad_1_x");
    expect(folderGeneKey("Gad-1 (x)")).toBe("/genes/Gad_1_x.bin.gz");
    expect(folderGeneKey("A", true)).toBe("/genes/A_uuuuuuuuuu.bin.gz");
  });

  it("opens from the manifest only, then loads a gene whole with progress", async () => {
    const a = coords(5, 1);
    const { store, reads } = folder({
      "/manifest.json.gz": manifest({ A: 5, B: 2 }),
      "/genes/A.bin.gz": gzipped(a),
    });
    const ds = await SingleMoleculeDataset.fromElementStore(store, {
      id: "ds",
      geneKey: folderGeneKey,
    });

    expect(ds.uniqueGenes).toEqual(["A", "B"]);
    expect(ds.scalingFactor).toBe(2);
    expect(reads).toEqual(["get /manifest.json.gz"]);

    const progress = vi.fn();
    const got = await ds.getCoordinatesByGene("A", { onProgress: progress });

    expect(got).toEqual(a);
    expect(progress).toHaveBeenCalled();
    expect(progress.mock.lastCall![1]).toBe(gzipped(a).length);
    // Cached: no second read
    expect(await ds.getCoordinatesByGene("A")).toBe(got);
    expect(reads.filter((r) => r.startsWith("open"))).toEqual([
      "open /genes/A.bin.gz",
    ]);
    await expect(ds.getCoordinatesByGene("B")).rejects.toThrow(
      /No molecule file/,
    );
    await expect(ds.getCoordinatesByGene("Z")).rejects.toThrow(/not found/);
  });

  it("streams a gene in batches and denormalizes an old folder", async () => {
    const a = coords(600_000, 0);
    const { store } = folder({
      "/manifest.json.gz": manifest({ A: 600_000 }, false),
      "/genes/A.bin.gz": gzipped(a),
    });
    const ds = await SingleMoleculeDataset.fromElementStore(store, {
      id: "ds",
      geneKey: folderGeneKey,
    });
    const batches: number[] = [];
    const all = await ds.streamCoordinatesByGene("A", {
      onBatch: (b) => {
        batches.push(b.coords.length / 3);
      },
    });

    expect(batches.length).toBeGreaterThan(1);
    expect(batches.reduce((x, y) => x + y, 0)).toBe(600_000);
    // Normalized coordinates come back scaled by the folder's factor
    expect(all[3]).toBeCloseTo(a[3] * 2);
    expect(ds.scalingFactor).toBe(1);
  });

  it("counts a gene from the gzip trailer when the manifest does not say", async () => {
    const u = coords(7);
    const { store, reads } = folder({
      "/manifest.json.gz": manifest({ A: 3 }),
      "/genes/A_uuuuuuuuuu.bin.gz": gzipped(u),
    });
    const ds = await SingleMoleculeDataset.fromElementStore(store, {
      id: "ds",
      geneKey: folderGeneKey,
    });

    expect(await ds.resolveMoleculeCount("A")).toBe(3);
    expect(await ds.resolveMoleculeCount("A", "unassigned")).toBe(7);
    expect(await ds.resolveMoleculeCount("A", "unassigned")).toBe(7);
    expect(reads.filter((r) => r.startsWith("range"))).toEqual([
      "range /genes/A_uuuuuuuuuu.bin.gz",
    ]);
    expect(await ds.getUnassignedCoordinatesByGene("A")).toEqual(u);
    // A gene with no unassigned file has none
    expect(await ds.resolveMoleculeCount("Z", "unassigned")).toBeNull();
  });

  it("gives an empty unassigned set for a missing file, and a cancelled read fails", async () => {
    const { store } = folder({ "/manifest.json.gz": manifest({ A: 1 }) });
    const ds = await SingleMoleculeDataset.fromElementStore(store, {
      id: "ds",
      geneKey: folderGeneKey,
    });

    expect((await ds.getUnassignedCoordinatesByGene("A")).length).toBe(0);
    const batches: unknown[] = [];

    expect(
      (
        await ds.streamUnassignedCoordinatesByGene("A", {
          onBatch: (b) => {
            batches.push(b);
          },
        })
      ).length,
    ).toBe(0);
    expect(batches.length).toBe(1);

    const ctl = new AbortController();

    ctl.abort();
    await expect(
      ds.getCoordinatesByGene("A", { signal: ctl.signal }),
    ).rejects.toMatchObject({ name: "AbortError" });
  });
});
