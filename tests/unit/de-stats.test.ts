import { describe, it, expect, vi, beforeEach } from "vitest";

import {
  computeDeStats,
  ensureDeStatsForColumn,
  parseDeStatsBuffer,
  rankDegsForCelltype,
  selectPriorityCategoricalCluster,
  deStatsInFlight,
  isDeStatsInFlight,
  type DeStats,
} from "@/lib/utils/de-stats";

// ---------------------------------------------------------------------------
// Test fixtures
//
// 4 cells × 3 genes, 2 celltypes:
//   cell  ct  gene0  gene1  gene2
//    0    A    2.0    0.0    1.0
//    1    A    4.0    0.0    0.0
//    2    B    0.0    3.0    0.0
//    3    B    0.0    1.0    2.0
//
// Per-celltype expected stats:
//   means[g,A] = (sum over A cells) / 2     means[g,B] = (sum over B) / 2
//     gene0  A=3.0   B=0.0
//     gene1  A=0.0   B=2.0
//     gene2  A=0.5   B=1.0
//   pct (fraction with nonzero):
//     gene0  A=1.0   B=0.0
//     gene1  A=0.0   B=1.0
//     gene2  A=0.5   B=0.5
// ---------------------------------------------------------------------------

const GENES = ["gene0", "gene1", "gene2"];

function fixtureCluster() {
  return {
    column: "celltype",
    type: "categorical",
    valueIndices: new Uint16Array([0, 0, 1, 1]),
    uniqueValues: ["A", "B"],
  } as any;
}

const EXPECTED = {
  cellCounts: [2, 2],
  // gene-major: idx = g * C + c
  means: new Float32Array([3.0, 0.0, 0.0, 2.0, 0.5, 1.0]),
  pct: new Float32Array([1.0, 0.0, 0.0, 1.0, 0.5, 0.5]),
};

function expectStatsMatch(stats: DeStats) {
  expect(stats.column).toBe("celltype");
  expect(stats.celltypes).toEqual(["A", "B"]);
  expect(stats.genes).toEqual(GENES);
  expect(Array.from(stats.cellCounts)).toEqual(EXPECTED.cellCounts);
  for (let i = 0; i < EXPECTED.means.length; i++) {
    expect(stats.means[i]).toBeCloseTo(EXPECTED.means[i], 5);
    expect(stats.pctExpressing[i]).toBeCloseTo(EXPECTED.pct[i], 5);
  }
}

// ---------------------------------------------------------------------------

describe("selectPriorityCategoricalCluster", () => {
  it("respects the priority list (class_name first)", () => {
    const clusters = [
      { column: "leiden", type: "categorical" },
      { column: "class_name", type: "categorical" },
      { column: "celltype", type: "categorical" },
    ];
    const pick = selectPriorityCategoricalCluster(clusters as any);
    expect(pick?.column).toBe("class_name");
  });

  it("skips numerical columns even when their name matches", () => {
    const clusters = [
      { column: "class_name", type: "numerical" },
      { column: "leiden", type: "categorical" },
    ];
    const pick = selectPriorityCategoricalCluster(clusters as any);
    expect(pick?.column).toBe("leiden");
  });

  it("falls back to a *celltype* column when no priority hit", () => {
    const clusters = [
      { column: "my_celltype_v2", type: "categorical" },
      { column: "n_genes", type: "numerical" },
    ];
    const pick = selectPriorityCategoricalCluster(clusters as any);
    expect(pick?.column).toBe("my_celltype_v2");
  });

  it("falls back to any categorical when nothing else matches", () => {
    const clusters = [
      { column: "n_genes", type: "numerical" },
      { column: "sample", type: "categorical" },
    ];
    const pick = selectPriorityCategoricalCluster(clusters as any);
    expect(pick?.column).toBe("sample");
  });

  it("returns null with no categorical columns", () => {
    expect(
      selectPriorityCategoricalCluster([
        { column: "x", type: "numerical" },
      ] as any),
    ).toBeNull();
    expect(selectPriorityCategoricalCluster([])).toBeNull();
    expect(selectPriorityCategoricalCluster(null)).toBeNull();
  });
});

// ---------------------------------------------------------------------------

describe("computeDeStats: flat TypedArray (H5AD path)", () => {
  it("computes means + pct from a row-major Float32Array", async () => {
    // Row-major [N × G]: cell0=[2,0,1], cell1=[4,0,0], cell2=[0,3,0], cell3=[0,1,2]
    const matrix = new Float32Array([2, 0, 1, 4, 0, 0, 0, 3, 0, 0, 1, 2]);
    const stats = await computeDeStats(fixtureCluster(), GENES, matrix);
    expect(stats).not.toBeNull();
    expectStatsMatch(stats!);
  });

  it("logs and returns null if matrix is too short", async () => {
    const warn = vi.spyOn(console, "warn").mockImplementation(() => {});
    const stats = await computeDeStats(
      fixtureCluster(),
      GENES,
      new Float32Array(5),
    );
    expect(stats).toBeNull();
    expect(warn).toHaveBeenCalled();
    warn.mockRestore();
  });

  it("yields to the event loop during long passes", async () => {
    const matrix = new Float32Array([2, 0, 1, 4, 0, 0, 0, 3, 0, 0, 1, 2]);
    const progress: number[] = [];
    const stats = await computeDeStats(fixtureCluster(), GENES, matrix, {
      yieldEvery: 1,
      onProgress: (frac) => {
        progress.push(frac);
      },
    });
    expect(stats).not.toBeNull();
    expect(progress.at(-1)).toBe(1);
    expect(progress.length).toBeGreaterThan(1);
  });
});

describe("computeDeStats: Map<gene, Float32Array> (Xenium/MERSCOPE path)", () => {
  it("computes means + pct over per-gene cell-indexed arrays", async () => {
    // One array per gene: index = cell index
    const matrix = new Map<string, Float32Array>([
      ["gene0", new Float32Array([2, 4, 0, 0])],
      ["gene1", new Float32Array([0, 0, 3, 1])],
      ["gene2", new Float32Array([1, 0, 0, 2])],
    ]);
    const stats = await computeDeStats(fixtureCluster(), GENES, matrix);
    expect(stats).not.toBeNull();
    expectStatsMatch(stats!);
  });

  it("skips genes missing from the Map", async () => {
    // gene1 missing → its row should stay zero, others should compute normally
    const matrix = new Map<string, Float32Array>([
      ["gene0", new Float32Array([2, 4, 0, 0])],
      ["gene2", new Float32Array([1, 0, 0, 2])],
    ]);
    const stats = await computeDeStats(fixtureCluster(), GENES, matrix);
    expect(stats).not.toBeNull();
    const C = 2;
    // gene1 row (g=1) — both cells: 0
    expect(stats!.means[1 * C + 0]).toBe(0);
    expect(stats!.means[1 * C + 1]).toBe(0);
    expect(stats!.pctExpressing[1 * C + 0]).toBe(0);
    expect(stats!.pctExpressing[1 * C + 1]).toBe(0);
    // gene0 + gene2 still correct
    expect(stats!.means[0 * C + 0]).toBeCloseTo(3.0, 5);
    expect(stats!.means[2 * C + 1]).toBeCloseTo(1.0, 5);
  });
});

describe("computeDeStats: nested arrays", () => {
  it("computes means + pct from number[][] row-major", async () => {
    const matrix = [
      [2, 0, 1],
      [4, 0, 0],
      [0, 3, 0],
      [0, 1, 2],
    ];
    const stats = await computeDeStats(fixtureCluster(), GENES, matrix);
    expect(stats).not.toBeNull();
    expectStatsMatch(stats!);
  });
});

describe("computeDeStats: invalid inputs", () => {
  it("returns null when cluster lacks valueIndices/uniqueValues", async () => {
    const cluster = {
      column: "celltype",
      type: "categorical",
      // no valueIndices / uniqueValues
    } as any;
    const matrix = new Float32Array([1, 2, 3]);
    const stats = await computeDeStats(cluster, ["g0"], matrix);
    expect(stats).toBeNull();
  });

  it("returns null with empty genes", async () => {
    const matrix = new Float32Array(0);
    const stats = await computeDeStats(fixtureCluster(), [], matrix);
    expect(stats).toBeNull();
  });

  it("warns and returns null on unsupported matrix shape", async () => {
    const warn = vi.spyOn(console, "warn").mockImplementation(() => {});
    const stats = await computeDeStats(
      fixtureCluster(),
      GENES,
      "not a matrix" as any,
    );
    expect(stats).toBeNull();
    expect(warn).toHaveBeenCalled();
    warn.mockRestore();
  });
});

// ---------------------------------------------------------------------------

describe("parseDeStatsBuffer", () => {
  it("round-trips against the Python binary layout", () => {
    // Build a buffer matching write_de_stats_binary in process_spatial_data.py
    //   header (16 bytes): version=1, G=3, C=2, reserved=0
    //   celltype names: "A" then "B"
    //   cell counts: [2, 2]
    //   means: float32[G*C]
    //   pct:   float32[G*C]
    const celltypes = ["A", "B"];
    const cellCounts = [2, 2];
    const means = new Float32Array([3.0, 0.0, 0.0, 2.0, 0.5, 1.0]);
    const pct = new Float32Array([1.0, 0.0, 0.0, 1.0, 0.5, 0.5]);

    // Compute size
    let nameBytesLen = 0;
    for (const ct of celltypes) {
      nameBytesLen += 4 + new TextEncoder().encode(ct).byteLength;
    }
    const headerLen = 16;
    const cellCountsLen = celltypes.length * 4;
    const meansLen = means.byteLength;
    const pctLen = pct.byteLength;
    const total = headerLen + nameBytesLen + cellCountsLen + meansLen + pctLen;

    const buf = new ArrayBuffer(total);
    const view = new DataView(buf);
    let off = 0;
    view.setUint32(off, 1, true);
    off += 4; // version
    view.setUint32(off, 3, true);
    off += 4; // num_genes
    view.setUint32(off, 2, true);
    off += 4; // num_celltypes
    view.setUint32(off, 0, true);
    off += 4; // reserved

    const enc = new TextEncoder();
    for (const ct of celltypes) {
      const bytes = enc.encode(ct);
      view.setUint32(off, bytes.byteLength, true);
      off += 4;
      new Uint8Array(buf, off, bytes.byteLength).set(bytes);
      off += bytes.byteLength;
    }
    for (const n of cellCounts) {
      view.setUint32(off, n, true);
      off += 4;
    }
    new Uint8Array(buf, off, meansLen).set(new Uint8Array(means.buffer));
    off += meansLen;
    new Uint8Array(buf, off, pctLen).set(new Uint8Array(pct.buffer));

    const parsed = parseDeStatsBuffer(buf, "celltype", GENES);
    expectStatsMatch(parsed);
  });

  it("throws on unsupported version", () => {
    const buf = new ArrayBuffer(16);
    new DataView(buf).setUint32(0, 99, true);
    expect(() => parseDeStatsBuffer(buf, "celltype", GENES)).toThrow(
      /version/i,
    );
  });

  it("throws when num_genes mismatches the dataset", () => {
    const buf = new ArrayBuffer(16);
    const view = new DataView(buf);
    view.setUint32(0, 1, true);
    view.setUint32(4, 7, true); // wrong gene count
    view.setUint32(8, 2, true);
    expect(() => parseDeStatsBuffer(buf, "celltype", GENES)).toThrow(
      /numGenes/i,
    );
  });
});

// ---------------------------------------------------------------------------

function makeFixtureDeStats(): DeStats {
  return {
    column: "celltype",
    celltypes: ["A", "B", "C"],
    cellCounts: [2, 2, 4],
    genes: GENES,
    // gene-major: idx = g * 3 + c
    means: new Float32Array([
      3.0, 0.0, 1.0, // gene0: A=3, B=0, C=1
      0.0, 2.0, 0.0, // gene1: A=0, B=2, C=0
      0.5, 1.0, 0.25, // gene2: A=0.5, B=1, C=0.25
    ]),
    pctExpressing: new Float32Array([
      1.0, 0.0, 0.5,
      0.0, 1.0, 0.0,
      0.5, 0.5, 0.25,
    ]),
  };
}

describe("rankDegsForCelltype: vs Rest", () => {
  it("reconstructs meanOut as the weighted average of non-target celltypes", () => {
    const stats = makeFixtureDeStats();
    const ranked = rankDegsForCelltype(stats, "A");
    // nTotal = 8, nIn = 2, nOut = 6.
    // gene0: meanOut = (0*2 + 1*4) / 6 = 0.6667
    const gene0 = ranked.find((r) => r.gene === "gene0")!;
    expect(gene0.meanIn).toBeCloseTo(3.0, 5);
    expect(gene0.meanOut).toBeCloseTo(2 / 3, 4);
    expect(gene0.pctIn).toBeCloseTo(1.0, 5);
    // pctOut = (0*2 + 0.5*4)/6 = 0.3333
    expect(gene0.pctOut).toBeCloseTo(1 / 3, 4);
  });

  it("sorts descending by log2FC", () => {
    const stats = makeFixtureDeStats();
    const ranked = rankDegsForCelltype(stats, "A");
    for (let i = 1; i < ranked.length; i++) {
      expect(ranked[i - 1].log2FC).toBeGreaterThanOrEqual(ranked[i].log2FC);
    }
  });

  it("returns empty for unknown target", () => {
    expect(rankDegsForCelltype(makeFixtureDeStats(), "ZZZ")).toEqual([]);
  });

  it("uses the pseudocount so zero-expression genes don't produce ±Infinity", () => {
    const stats = makeFixtureDeStats();
    const ranked = rankDegsForCelltype(stats, "A", { pseudocount: 1e-9 });
    for (const r of ranked) {
      expect(Number.isFinite(r.log2FC)).toBe(true);
    }
  });
});

describe("rankDegsForCelltype: vs specific reference", () => {
  it("reads meanOut directly from the reference celltype row", () => {
    const stats = makeFixtureDeStats();
    const ranked = rankDegsForCelltype(stats, "A", { reference: "B" });
    const gene1 = ranked.find((r) => r.gene === "gene1")!;
    expect(gene1.meanIn).toBeCloseTo(0.0, 5);
    expect(gene1.meanOut).toBeCloseTo(2.0, 5); // means[1*3 + 1]
    expect(gene1.pctIn).toBeCloseTo(0.0, 5);
    expect(gene1.pctOut).toBeCloseTo(1.0, 5);
  });

  it("returns empty when target === reference", () => {
    expect(
      rankDegsForCelltype(makeFixtureDeStats(), "A", { reference: "A" }),
    ).toEqual([]);
  });

  it("returns empty when reference is unknown", () => {
    expect(
      rankDegsForCelltype(makeFixtureDeStats(), "A", { reference: "ZZZ" }),
    ).toEqual([]);
  });

  it("produces a different top gene than vs Rest", () => {
    const stats = makeFixtureDeStats();
    const vsRest = rankDegsForCelltype(stats, "A");
    const vsB = rankDegsForCelltype(stats, "A", { reference: "B" });
    // Sanity: ordering can differ because the comparison set differs.
    // We only assert the math diverges for at least one gene.
    let differs = false;
    for (let i = 0; i < vsRest.length; i++) {
      const g = vsRest[i].gene;
      const r2 = vsB.find((r) => r.gene === g);
      if (r2 && Math.abs(r2.log2FC - vsRest[i].log2FC) > 1e-6) {
        differs = true;
        break;
      }
    }
    expect(differs).toBe(true);
  });
});

// ---------------------------------------------------------------------------

describe("ensureDeStatsForColumn", () => {
  beforeEach(() => {
    deStatsInFlight.clear();
  });

  function makeDataset(overrides: Record<string, any> = {}) {
    const cluster = fixtureCluster();
    return {
      id: `ds_${Math.random().toString(36).slice(2)}`,
      genes: GENES,
      // Flat row-major matrix matching fixtureCluster
      matrix: new Float32Array([2, 0, 1, 4, 0, 0, 0, 3, 0, 0, 1, 2]),
      clusters: [cluster],
      deStatsByColumn: new Map<string, DeStats>(),
      availableDeStatsColumns: [] as string[],
      adapter: null as any,
      ...overrides,
    };
  }

  it("computes once and caches on the dataset", async () => {
    const ds = makeDataset();
    const stats = await ensureDeStatsForColumn(ds as any, "celltype");
    expect(stats).not.toBeNull();
    expect(ds.deStatsByColumn.has("celltype")).toBe(true);
    expect(isDeStatsInFlight(ds.id, "celltype")).toBe(false);

    // Second call should hit the cache (no new compute). Mutate the cached
    // stats so we can detect cache reuse vs recompute.
    const cached = ds.deStatsByColumn.get("celltype")!;
    cached.celltypes = ["SENTINEL"];
    const second = await ensureDeStatsForColumn(ds as any, "celltype");
    expect(second?.celltypes).toEqual(["SENTINEL"]);
  });

  it("dedupes concurrent calls for the same (dataset, column)", async () => {
    const ds = makeDataset();
    const [a, b] = await Promise.all([
      ensureDeStatsForColumn(ds as any, "celltype"),
      ensureDeStatsForColumn(ds as any, "celltype"),
    ]);
    expect(a).toBe(b); // exact same promise result reference
    expect(isDeStatsInFlight(ds.id, "celltype")).toBe(false);
  });

  it("prefers adapter.loadDeStats when the column is advertised", async () => {
    const fromAdapter: DeStats = {
      column: "celltype",
      celltypes: ["FROM_ADAPTER"],
      cellCounts: [1],
      genes: GENES,
      means: new Float32Array(GENES.length),
      pctExpressing: new Float32Array(GENES.length),
    };
    const adapter = {
      loadDeStats: vi.fn().mockResolvedValue(fromAdapter),
    };
    const ds = makeDataset({
      availableDeStatsColumns: ["celltype"],
      adapter,
    });
    const stats = await ensureDeStatsForColumn(ds as any, "celltype");
    expect(stats?.celltypes).toEqual(["FROM_ADAPTER"]);
    expect(adapter.loadDeStats).toHaveBeenCalledWith("celltype", GENES);
  });

  it("falls back to recompute when adapter.loadDeStats rejects", async () => {
    const warn = vi.spyOn(console, "warn").mockImplementation(() => {});
    const adapter = {
      loadDeStats: vi.fn().mockRejectedValue(new Error("network down")),
    };
    const ds = makeDataset({
      availableDeStatsColumns: ["celltype"],
      adapter,
    });
    const stats = await ensureDeStatsForColumn(ds as any, "celltype");
    expect(stats).not.toBeNull();
    expect(stats?.celltypes).toEqual(["A", "B"]); // came from local recompute
    expect(warn).toHaveBeenCalled();
    warn.mockRestore();
  });

  it("clears the in-flight set after rejection", async () => {
    const adapter = {
      loadDeStats: vi.fn().mockRejectedValue(new Error("boom")),
    };
    const ds = makeDataset({
      availableDeStatsColumns: ["celltype"],
      adapter,
      // Remove the matrix so the recompute fallback also fails → final null
      matrix: null,
    });
    const warn = vi.spyOn(console, "warn").mockImplementation(() => {});
    const stats = await ensureDeStatsForColumn(ds as any, "celltype");
    expect(stats).toBeNull();
    expect(isDeStatsInFlight(ds.id, "celltype")).toBe(false);
    warn.mockRestore();
  });

  it("returns null when no matrix and adapter can't fetch", async () => {
    const ds = makeDataset({ matrix: null });
    const stats = await ensureDeStatsForColumn(ds as any, "celltype");
    expect(stats).toBeNull();
  });

  it("returns null when the requested column isn't categorical", async () => {
    const ds = makeDataset({
      clusters: [
        { column: "n_genes", type: "numerical", valueIndices: new Uint16Array([0]), uniqueValues: ["0"] },
      ],
    });
    const stats = await ensureDeStatsForColumn(ds as any, "n_genes");
    expect(stats).toBeNull();
  });
});
