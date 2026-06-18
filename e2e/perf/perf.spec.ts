import {
  test,
  expect,
  datasetsByTier,
  isMaterialized,
  annotate,
  loadLocalDataset,
  selectGeneAndMeasure,
  recordSamples,
  type DatasetEntry,
} from "../fixtures/merfish-fixtures";

/**
 * Performance suite — tagged @perf. Takes N samples per metric and records the
 * median to perf/results/latest.json (via the perf reporter) for baseline
 * comparison. Includes tiny by default; medium/large are included only when
 * materialized on disk (run `npm run test:data:fetch -- --tier medium`).
 *
 * Each sample runs in a FRESH browser context: very large datasets (e.g.
 * sm-large, 20M molecules) accumulate GPU/JS memory across loads, and reusing
 * one tab eventually crashes the renderer. A clean context per load isolates
 * memory and keeps timings independent.
 */

const SAMPLES = Number(process.env.PERF_SAMPLES ?? 3);
// browser.newContext() does not inherit `use.baseURL`, so set it explicitly to
// match playwright.config.ts.
const PORT = process.env.E2E_PORT || "3100";
const BASE_URL = process.env.E2E_BASE_URL || `http://localhost:${PORT}`;

const candidates: DatasetEntry[] = [
  ...datasetsByTier("tiny"),
  ...datasetsByTier("medium", "large"),
];

for (const ds of candidates) {
  test(`@perf timings: ${ds.id} (${ds.format})`, async ({ browser }, testInfo) => {
    annotate(testInfo, ds.id, ds.format);

    if (!isMaterialized(ds)) {
      test.skip(true, `${ds.id} not on disk — run \`npm run test:data:fetch -- ${ds.id}\``);
    }

    // Run one action in a throwaway context (fresh tab → memory reset each time).
    const inFreshContext = async <T>(
      fn: (page: import("@playwright/test").Page) => Promise<T>,
    ): Promise<T> => {
      const context = await browser.newContext({ baseURL: BASE_URL });
      try {
        return await fn(await context.newPage());
      } finally {
        await context.close();
      }
    };

    // Warmup (uncounted) so route compilation / first-paint cost in a dev server
    // doesn't skew the first measured sample.
    await inFreshContext((page) => loadLocalDataset(page, ds));

    // Sample load → render time, each in its own fresh context.
    const loadSamples: number[] = [];
    for (let i = 0; i < SAMPLES; i++) {
      const ms = await inFreshContext(async (page) => {
        const timing = await loadLocalDataset(page, ds);
        return timing.loadToRender;
      });
      loadSamples.push(ms);
    }
    await recordSamples(
      testInfo,
      { datasetId: ds.id, format: ds.format, metric: "loadToRender", tier: ds.tier },
      loadSamples,
    );
    expect(Math.min(...loadSamples)).toBeGreaterThan(0);

    // Sample gene-selection latency — load then select in the same fresh context.
    const gene = ds.expected.knownGenes[0];
    const geneSamples: number[] = [];
    for (let i = 0; i < SAMPLES; i++) {
      const ms = await inFreshContext(async (page) => {
        await loadLocalDataset(page, ds);
        return selectGeneAndMeasure(page, ds, gene);
      });
      geneSamples.push(ms);
    }
    await recordSamples(
      testInfo,
      {
        datasetId: ds.id,
        format: ds.format,
        metric: "geneSelectionLatency",
        tier: ds.tier,
      },
      geneSamples,
    );
  });
}
