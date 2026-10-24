import fs from "node:fs";
import path from "node:path";

import {
  test,
  expect,
  readHooks,
  recordMetric,
  SELECTORS,
} from "../fixtures/merfish-fixtures";

/**
 * Real SpatialData stores dropped for a browser preview: how long until the
 * viewer opens and first draws, how much memory that takes, and whether each
 * element the store holds is shown. The stores are on local disk, so this
 * runs only when they are named:
 *
 *   MERFISHEYES_REAL_STORES=/path/a.zarr,/path/b.zarr npx playwright test e2e/perf/real-drop.spec.ts
 *
 * Opened through the folder picker, which hands the browser the files
 * themselves (the drop zone's walk is covered by drop-matrix.spec.ts on the
 * tiny stores). Tagged @perf.
 */

const STORES = (process.env.MERFISHEYES_REAL_STORES ?? "")
  .split(",")
  .filter(Boolean);

test.skip(STORES.length === 0, "Set MERFISHEYES_REAL_STORES to run.");

for (const store of STORES) {
  // Named by its whole path: stores are often all called setN/setN.zarr
  const name = store;

  test(`@perf real drop: ${name}`, async ({ page }, info) => {
    test.setTimeout(900_000);
    const errors: string[] = [];

    page.on("pageerror", (e) => errors.push(e.message));
    await page.goto("/", { waitUntil: "domcontentloaded" });

    const start = Date.now();

    await page.locator("#sc-file-input-folder").setInputFiles(store);
    await page.waitForURL(/\/(viewer|sm-viewer|image-viewer)\/from-local/, {
      timeout: 600_000,
    });
    const toViewer = Date.now() - start;

    await page.waitForFunction(
      () => (window as any).__merfish?.renderComplete != null,
      undefined,
      { timeout: 600_000 },
    );
    const toRender = Date.now() - start;

    // The image and molecules draw after the cells; wait until each has
    const kinds = ["images", "points"].filter((k) =>
      fs.existsSync(path.join(store, k)),
    );

    for (const kind of kinds) {
      const layer = kind === "images" ? "image" : "points";

      await expect.soft
        .poll(async () => (await readHooks(page))?.layers[layer]?.drawn ?? 0, {
          message: `${layer} drawn`,
          timeout: 120_000,
        })
        .toBeGreaterThan(0);
    }
    const toAll = Date.now() - start;

    await recordMetric(info, {
      datasetId: name,
      format: "spatialdata-drop",
      metric: "dropToAllLayers",
      ms: toAll,
    });
    const hooks = (await readHooks(page))!;
    const record = (metric: string, ms: number) =>
      recordMetric(info, {
        datasetId: name,
        format: "spatialdata-drop",
        metric,
        ms,
      });

    await record("dropToViewer", toViewer);
    await record("dropToRender", toRender);
    await info.attach("summary", {
      contentType: "application/json",
      body: JSON.stringify({
        store,
        toViewer,
        toRender,
        toAll,
        stats: hooks.stats,
        jsHeapMB: (hooks.memory?.jsHeap ?? 0) / 2 ** 20,
        layers: hooks.layers,
        errors,
      }),
    });
    expect.soft(errors, "page errors").toEqual([]);
    expect.soft(hooks.stats?.pointCount ?? 0, "cells").toBeGreaterThan(0);
    await expect.soft(page.getByText(/Skipped:/)).toHaveCount(0);
    // Whatever else the store holds should reach the screen
    if (fs.existsSync(path.join(store, "images"))) {
      await expect
        .soft(page.locator(SELECTORS.imageButton), "the image")
        .toBeVisible({ timeout: 60_000 });
    }
    if (fs.existsSync(path.join(store, "points"))) {
      await expect
        .soft(
          page
            .getByTestId("points-status")
            .or(page.locator(SELECTORS.smLayerButton))
            .first(),
          "the molecules",
        )
        .toBeVisible({ timeout: 60_000 });
    }
  });
}
