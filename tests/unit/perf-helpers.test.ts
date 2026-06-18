import { describe, it, expect, vi, afterEach } from "vitest";

import { median } from "@/e2e/helpers/perf";
import {
  pctChange,
  classify,
  recordKey,
  baselinePath,
} from "@/scripts/perf/compare-baseline";

describe("perf.median", () => {
  it("returns the middle value for odd-length input", () => {
    expect(median([3, 1, 2])).toBe(2);
    expect(median([5])).toBe(5);
  });

  it("averages the two middle values for even-length input", () => {
    expect(median([1, 2, 3, 4])).toBe(2.5);
    expect(median([10, 20])).toBe(15);
  });

  it("is order-independent", () => {
    expect(median([4, 2, 1, 3])).toBe(median([1, 2, 3, 4]));
  });

  it("returns NaN for empty input", () => {
    expect(Number.isNaN(median([]))).toBe(true);
  });
});

describe("baseline math", () => {
  it("pctChange computes signed percent change vs baseline", () => {
    expect(pctChange(100, 120)).toBeCloseTo(20);
    expect(pctChange(100, 80)).toBeCloseTo(-20);
    expect(pctChange(50, 75)).toBeCloseTo(50);
    expect(pctChange(100, 100)).toBeCloseTo(0);
  });

  it("classify maps percent change to status by threshold", () => {
    // warn > 10%, fail > 20%
    expect(classify(0, 10, 20)).toBe("OK");
    expect(classify(15, 10, 20)).toBe("WARN");
    expect(classify(25, 10, 20)).toBe("FAIL");
    expect(classify(-15, 10, 20)).toBe("FAST");
  });

  it("classify thresholds are strict (boundary values)", () => {
    expect(classify(10, 10, 20)).toBe("OK"); // not > warn
    expect(classify(20, 10, 20)).toBe("WARN"); // > warn, not > fail
    expect(classify(-10, 10, 20)).toBe("OK"); // not < -warn
  });
});

describe("baseline keys / paths", () => {
  it("recordKey namespaces by browser/dataset/metric", () => {
    expect(
      recordKey({ browser: "chromium", datasetId: "h5ad-tiny", metric: "loadToRender" }),
    ).toBe("chromium/h5ad-tiny/loadToRender");
  });

  it("recordKey falls back to 'unknown' browser", () => {
    expect(recordKey({ datasetId: "sm-tiny", metric: "geneSelectionLatency" })).toBe(
      "unknown/sm-tiny/geneSelectionLatency",
    );
  });

  it("baselinePath is per-platform", () => {
    expect(baselinePath("linux").endsWith("linux.baseline.json")).toBe(true);
    expect(baselinePath("darwin").endsWith("darwin.baseline.json")).toBe(true);
  });
});

describe("test-hooks (window.__merfish)", () => {
  it("is inert when NEXT_PUBLIC_E2E is not set", async () => {
    vi.resetModules();
    delete process.env.NEXT_PUBLIC_E2E;
    (globalThis as any).window = {};
    const m = await import("@/lib/utils/test-hooks");
    m.resetHooks();
    m.markRenderComplete({ pointCount: 1 });
    expect((globalThis as any).window.__merfish).toBeUndefined();
  });

  it("records timing + stats when enabled", async () => {
    vi.resetModules();
    process.env.NEXT_PUBLIC_E2E = "1";
    (globalThis as any).window = {};
    const m = await import("@/lib/utils/test-hooks");

    m.resetHooks();
    m.markSceneReady();
    m.markRenderComplete({ pointCount: 42, geneCount: 7, dataType: "single_cell" });
    m.markRenderComplete({ pointCount: 42 }); // second render bumps the counter
    m.markGeneRendered("gene_0", 5);

    const h = (globalThis as any).window.__merfish;
    expect(h.enabled).toBe(true);
    expect(h.sceneReady).toBeGreaterThanOrEqual(0);
    expect(h.renderComplete).toBeGreaterThanOrEqual(0);
    expect(h.renderCount).toBe(2);
    expect(h.stats.pointCount).toBe(42);
    expect(h.stats.geneCount).toBe(7);
    expect(h.lastGene).toEqual(
      expect.objectContaining({ gene: "gene_0", count: 5 }),
    );

    m.resetHooks();
    expect(h.renderCount).toBe(0);
    expect(h.renderComplete).toBeNull();
    expect(h.lastGene).toBeNull();
  });
});

// Keep the global clean for other suites.
afterEach(() => {
  delete (globalThis as any).window;
  delete process.env.NEXT_PUBLIC_E2E;
});
