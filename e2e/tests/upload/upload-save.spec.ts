import {
  test,
  expect,
  getDataset,
  annotate,
  assertCount,
  assertEqual,
  loadLocalDataset,
  uploadAndSave,
  openSavedDataset,
} from "../../fixtures/merfish-fixtures";

/**
 * Upload → save → reload tests. Tagged @upload; these need a reachable Postgres
 * + S3 (MinIO). Run via `npm run test:upload` after `npm run test:infra:up`.
 * Skipped unless E2E_UPLOAD=1 so the default/quick suite never hits services.
 */

const UPLOAD_ENABLED = process.env.E2E_UPLOAD === "1";

// One representative dataset per upload path.
const UPLOAD_TARGETS = ["h5ad-tiny", "sm-tiny"];

test.describe("@upload upload + save + reload", () => {
  test.skip(!UPLOAD_ENABLED, "Set E2E_UPLOAD=1 (and start docker infra) to run upload tests.");

  for (const id of UPLOAD_TARGETS) {
    const ds = getDataset(id);

    test(`upload/save/reload: ${ds.id} (${ds.format})`, async ({ page }, testInfo) => {
      annotate(testInfo, ds.id, ds.format);

      await test.step("load dataset locally", async () => {
        await loadLocalDataset(page, ds);
      });

      let datasetId = "";
      await test.step("upload + save + measure time", async () => {
        const res = await uploadAndSave(page, ds, { name: `e2e-${ds.id}` });
        datasetId = res.datasetId;
        expect(datasetId, "saved dataset id").toBeTruthy();
        testInfo.annotations.push({
          type: "timing",
          description: `uploadSaveTotal=${res.totalMs}ms`,
        });
      });

      await test.step("verify saved dataset is COMPLETE via API", async () => {
        const seg = ds.dataType === "single_molecule" ? "single-molecule" : "datasets";
        const resp = await page.request.get(`/api/${seg}/${datasetId}`);
        expect(resp.ok(), `GET /api/${seg}/${datasetId} -> ${resp.status()}`).toBeTruthy();
        const body = await resp.json();
        // Status field name varies by route; accept either shape.
        const status = body.status ?? body.dataset?.status;
        if (status) assertEqual(status, "COMPLETE", {
          dataset: ds.id,
          format: ds.format,
          step: "verify status COMPLETE",
        });
      });

      await test.step("reload saved dataset from S3 + verify counts", async () => {
        const hooks = await openSavedDataset(page, ds, datasetId);
        if (!hooks?.enabled || !hooks.stats) return;
        const ctx = { dataset: ds.id, format: ds.format, step: "reload verify counts" };
        assertEqual(hooks.stats.dimensions, ds.expected.dimensions, ctx);
        assertCount(hooks.stats.geneCount ?? -1, ds.expected.genes, ctx);
        if (ds.dataType === "single_cell" && ds.expected.cells != null) {
          assertCount(hooks.stats.pointCount ?? -1, ds.expected.cells, ctx);
        }
      });
    });
  }
});
