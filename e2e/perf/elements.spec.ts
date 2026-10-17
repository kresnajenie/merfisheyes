import fs from "node:fs";
import path from "node:path";

import {
  test,
  expect,
  startSpatialDataServer,
  storeUrl,
  recordMetric,
  REPO_ROOT,
  SELECTORS,
} from "../fixtures/merfish-fixtures";
import {
  firstContent,
  frames,
  layer,
  memory,
  nudgeSlider,
  pageNow,
  pan,
  resetFrames,
  settledAfter,
  sweep,
  zoom,
} from "../helpers/budgets";
import { type FixtureServer } from "../helpers/fixture-server";

/**
 * The budgets of docs/specs/BUDGETS.md, measured — tagged @perf, written to
 * perf/results/latest.json like the other perf metrics and gated against
 * the baseline. Every store is served by the local stand-in for a bucket at
 * a capped 12 MB/s (a 100 Mbit/s connection), so the numbers mean the same
 * in CI and on a laptop.
 *
 * Stores:
 *   tiny    the committed fixtures (30 cells, 3,000 molecules, an image,
 *           30 polygons) — always
 *   medium  a generated one-section store (scripts/testdata/make-medium-store.py)
 *           when it is on disk
 *
 * Metrics (ms unless noted):
 *   openToCells     navigation → the cells drawn
 *   firstMolecules  navigation → the first molecules on screen
 *   viewSettled     a zoom → every molecule the view wants drawn
 *   bytesPerView    bytes of points read for that zoom (reported in ms
 *                   slot as kilobytes, so the baseline gate applies)
 *   frameOver16Pct  share of frames over 20 ms during a one-second pan, in %
 *   sliderToFrame   a dot-size step → the next frame drawn
 *   hoverFrameMs    worst frame while the pointer sweeps the canvas
 */

const BANDWIDTH = 12_000_000;
const MEDIUM = path.join(
  REPO_ROOT,
  "tests",
  "data",
  "generated",
  "medium",
  "points-store",
  "store.zarr",
);

interface Store {
  id: string;
  url: (server: FixtureServer) => string;
  present: () => boolean;
  /** A gene heavy enough that the point budget bites, selected after open. */
  heavyGene?: string;
}

const stores: Store[] = [
  {
    id: "spatialdata-tiny",
    url: (server) =>
      storeUrl(server, {
        cells: true,
        points: true,
        image: true,
        shapes: true,
      }),
    present: () => true,
  },
  {
    id: "spatialdata-medium",
    url: (server) => `${server.url}/medium/store.zarr`,
    present: () => fs.existsSync(MEDIUM),
    heavyGene: "Gene000",
  },
];

// Millions of molecules on a software GPU make slow frames; give the
// scripted interactions room.
test.setTimeout(600_000);

let server: FixtureServer;

/** One molecule dot-size step → the frame after it, in ms. */
async function sliderStep(
  page: import("@playwright/test").Page,
): Promise<number> {
  const before = await frames(page);
  const at = await pageNow(page);

  await nudgeSlider(page, "Dot size");
  // The frame after the change (a software GPU can take seconds per frame)
  await expect
    .poll(async () => (await frames(page)).count, { timeout: 120_000 })
    .toBeGreaterThan(before.count + 1);

  return (await pageNow(page)) - at;
}

/** Bytes and requests the server sent for the points element since `from`. */
function pointsTraffic(from: number): { bytes: number; requests: number } {
  const reads = server.requests
    .slice(from)
    .filter((r) => r.path.includes("/points/"));

  return {
    bytes: reads.reduce((sum, r) => sum + r.bytes, 0),
    requests: reads.length,
  };
}

test.beforeAll(async () => {
  server = await startSpatialDataServer({ bytesPerSecond: BANDWIDTH });
});
test.afterAll(() => server.close());

for (const store of stores) {
  test(`@perf budgets: ${store.id}`, async ({ page }, testInfo) => {
    if (!store.present()) {
      test.skip(true, `${store.id} is not on disk`);
    }
    const record = (metric: string, ms: number) =>
      recordMetric(testInfo, {
        datasetId: store.id,
        format: "spatialdata",
        metric,
        ms: Math.round(ms * 10) / 10,
        tier: store.id.endsWith("tiny") ? "tiny" : "medium",
      });

    await page.goto(
      `/spatialdata-viewer/from-s3?url=${encodeURIComponent(store.url(server))}`,
      { waitUntil: "commit" },
    );
    const t0 = await pageNow(page);

    // Open → cells
    await expect(page.locator(SELECTORS.scCanvas).first()).toBeVisible({
      timeout: 120_000,
    });
    const cellsAt = await firstContent(page, "table");

    await record("openToCells", cellsAt - t0);

    // → first molecules (the default genes of the points element)
    const moleculesAt = await firstContent(page, "points");

    await record("firstMolecules", moleculesAt - t0);
    await settledAfter(page, "points", moleculesAt - 1);

    // Hover across the canvas and one dot-size step, on the light scene
    // the store opens with
    await resetFrames(page);
    await sweep(page);
    await record("hoverFrameMs", (await frames(page)).worstMs);
    await record("sliderToFrame", await sliderStep(page));

    // A heavy gene selected by hand: time to its first molecules, to a
    // settled view, and the bytes it cost
    if (store.heavyGene) {
      const gene = store.heavyGene;

      await page.locator('[aria-label="Single molecule"]').first().click();
      await page.locator('button:has-text("Gene")').first().click();
      await page.locator('input[placeholder="Search genes..."]').fill(gene);
      const geneFrom = server.requests.length;
      const selectAt = await pageNow(page);

      await page
        .locator("label", { hasText: new RegExp(`^${gene}`) })
        .first()
        .click();
      const geneChanged = await expect
        .poll(async () => (await layer(page, "points"))?.changedAt ?? 0, {
          timeout: 30_000,
        })
        .toBeGreaterThan(selectAt)
        .then(async () => (await layer(page, "points"))!.changedAt!);
      const geneSettled = await settledAfter(page, "points", geneChanged);

      await record("geneToSettled", geneSettled - selectAt);
      await record("geneBytesKB", pointsTraffic(geneFrom).bytes / 1024);
      await record("geneRequests", pointsTraffic(geneFrom).requests);
      // Close the picker so the canvas is free for the pointer
      await page.keyboard.press("Escape");
    }

    // A zoom in: the view changes, molecules are read for it
    const zoomFrom = server.requests.length;
    const zoomAt = await pageNow(page);

    // Deep enough that a gene drawn from a coarser level comes down to
    // every molecule in the region: a fifth of the width
    await zoom(page, -32);
    const changed = await expect
      .poll(async () => (await layer(page, "points"))?.changedAt ?? 0, {
        timeout: 10_000,
      })
      .toBeGreaterThan(zoomAt)
      .then(async () => (await layer(page, "points"))!.changedAt!);
    // Damping keeps the camera moving for many frames after the wheel; on
    // a software GPU those frames are slow, so allow the view a long time
    const settledAt = await settledAfter(page, "points", changed, 300_000);

    await record("viewSettled", settledAt - changed);
    await record("bytesPerViewKB", pointsTraffic(zoomFrom).bytes / 1024);
    await record("requestsPerView", pointsTraffic(zoomFrom).requests);

    // A one-second pan: frames
    await resetFrames(page);
    await pan(page, 0.3, store.heavyGene ? 10 : 30);
    const panned = await frames(page);

    await record(
      "frameOver16Pct",
      panned.count ? (100 * panned.over16ms) / panned.count : 0,
    );
    await record("worstFrameMs", panned.worstMs);
    // CPU time inside render per frame: the part the code controls even on
    // a software GPU (headless Chromium), where the frame interval is not
    await record(
      "avgRenderMs",
      panned.count ? panned.renderMs / panned.count : 0,
    );

    // The dot-size step again, now with the heavy gene drawn
    if (store.heavyGene) {
      await record("sliderToFrameHeavy", await sliderStep(page));
    }

    // Memory at the end, in MB, for the record
    const mem = await memory(page);

    if (mem.jsHeap !== null) await record("jsHeapMB", mem.jsHeap / 1048576);
    await record("drawCalls", mem.drawCalls);
  });
}
