import {
  test,
  expect,
  datasetsByTier,
  annotate,
  assertCount,
  assertEqual,
  loadLocalDataset,
  waitForRender,
  readHooks,
  assertCanvasRendered,
  selectGeneAndMeasure,
  dragDropFile,
  type DatasetEntry,
} from "../../fixtures/merfish-fixtures";

/**
 * Local viewing tests — one per committed tiny dataset (every supported
 * format). Tagged @quick so they form the fast PR suite. Each step is wrapped
 * in test.step() so a failure points at the exact stage; assertions carry the
 * dataset/format/expected/actual context.
 */

const datasets = datasetsByTier("tiny");

for (const ds of datasets) {
  test(`@quick local view: ${ds.id} (${ds.format})`, async ({ page }, testInfo) => {
    annotate(testInfo, ds.id, ds.format);

    await test.step("load dataset locally + measure render time", async () => {
      const timing = await loadLocalDataset(page, ds);
      expect(
        timing.loadToRender,
        `${ds.id}: loadToRender should be positive`,
      ).toBeGreaterThan(0);
      testInfo.annotations.push({
        type: "timing",
        description: `loadToRender=${timing.loadToRender.toFixed(0)}ms loadToScene=${
          timing.loadToScene?.toFixed(0) ?? "n/a"
        }ms`,
      });
    });

    await test.step("verify dataset size / counts", async () => {
      const hooks = await readHooks(page);
      if (!hooks?.enabled || !hooks.stats) {
        test.info().annotations.push({
          type: "warning",
          description: "hooks disabled — skipping count assertions",
        });
        return;
      }
      const ctx = { dataset: ds.id, format: ds.format, step: "verify counts" };
      assertEqual(hooks.stats.dataType, ds.dataType, ctx);
      assertEqual(hooks.stats.dimensions, ds.expected.dimensions, ctx);
      assertCount(hooks.stats.geneCount ?? -1, ds.expected.genes, ctx);
      if (ds.dataType === "single_cell" && ds.expected.cells != null) {
        assertCount(hooks.stats.pointCount ?? -1, ds.expected.cells, ctx);
      }
    });

    await test.step("verify spatial coordinates render (non-blank canvas)", async () => {
      await assertCanvasRendered(page, ds);
    });

    await test.step("select known genes + measure selection latency", async () => {
      const genes = ds.expected.knownGenes.slice(0, 3);
      for (const gene of genes) {
        const latency = await selectGeneAndMeasure(page, ds, gene);
        expect(
          latency,
          `${ds.id}: gene '${gene}' selection latency should be finite`,
        ).toBeGreaterThanOrEqual(0);
        testInfo.annotations.push({
          type: "timing",
          description: `gene ${gene} latency=${latency.toFixed(0)}ms`,
        });
      }
    });
  });
}

// Drag-and-drop is only meaningful for single-file formats.
const dragDatasets: DatasetEntry[] = datasets.filter(
  (d) => d.format === "h5ad" || d.format === "single-molecule",
);

for (const ds of dragDatasets) {
  test(`@quick drag-drop load: ${ds.id} (${ds.format})`, async ({ page }, testInfo) => {
    annotate(testInfo, ds.id, ds.format);
    const ms = await dragDropFile(page, ds);
    expect(ms, `${ds.id}: drag-drop-to-render should be positive`).toBeGreaterThan(0);
    await waitForRender(page, ds);
    await assertCanvasRendered(page, ds);
    testInfo.annotations.push({
      type: "timing",
      description: `dragDropToRender=${ms.toFixed(0)}ms`,
    });
  });
}
