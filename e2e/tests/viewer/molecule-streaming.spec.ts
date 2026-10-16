import type { Page } from "@playwright/test";

import {
  test,
  expect,
  readHooks,
  startSpatialDataServer,
  storeUrl,
  cellsUrl,
  largeMoleculesUrl,
  LARGE_MOLECULES,
  SELECTORS,
  type MerfishGeneLoad,
} from "../../fixtures/merfish-fixtures";
import { type FixtureServer } from "../../helpers/fixture-server";

/**
 * Molecule gene loading is one code path (lib/webgl/MoleculeLayer.ts) on
 * every surface that draws molecules. These open each surface on a generated
 * molecule folder whose two biggest genes are above the streaming threshold
 * and check, through the `window.__merfish.geneLoads` hooks, that those genes
 * stream in batch by batch while the small one loads whole — and that
 * closing a gene's loading toast cancels its download. Tagged @quick: no
 * database or S3, the folder is served by a local stand-in for a bucket.
 */

let server: FixtureServer;

test.beforeAll(async () => {
  server = await startSpatialDataServer();
});
test.afterAll(() => server.close());

const geneLoad = (page: Page, gene: string) =>
  page.evaluate(
    (g) =>
      (window as any).__merfish?.geneLoads?.[g] as MerfishGeneLoad | undefined,
    gene,
  );

/** Wait until a gene's load has finished, and return how it went. */
async function finishedLoad(
  page: Page,
  gene: string,
): Promise<MerfishGeneLoad> {
  await expect
    .poll(async () => (await geneLoad(page, gene))?.done, {
      message: `${gene} should finish loading`,
      timeout: 60_000,
    })
    .toBe(true);

  return (await geneLoad(page, gene))!;
}

/**
 * The two big genes streamed in several batches; the small one loaded whole.
 * Each gene was downloaded once — unless the scene is rebuilt while loading
 * (a split view opening), which restarts what was in flight.
 */
async function expectStreamedDefaults(
  page: Page,
  tag: string,
  { countDownloads = true } = {},
): Promise<void> {
  const [big, huge, small] = LARGE_MOLECULES.genes(tag);

  for (const gene of [big, huge]) {
    const load = await finishedLoad(page, gene);

    expect(load, gene).toMatchObject({
      mode: "stream",
      expected: LARGE_MOLECULES.large,
      loaded: LARGE_MOLECULES.large,
      aborted: false,
      failed: false,
    });
    expect(load.batches, `${gene} batches`).toBeGreaterThanOrEqual(2);
  }
  expect(await finishedLoad(page, small)).toMatchObject({
    mode: "whole",
    loaded: LARGE_MOLECULES.small,
    batches: 1,
    aborted: false,
    failed: false,
  });
  // The fourth gene is not a default; nothing was fetched for it.
  expect(await geneLoad(page, LARGE_MOLECULES.genes(tag)[3])).toBeUndefined();
  if (countDownloads) {
    expect(
      server.requests.filter((r) =>
        r.path.startsWith(`/molecules-large/${tag}/genes/`),
      ).length,
      "one download per shown gene",
    ).toBe(3);
  }
}

/** Close the toast of a gene whose download is still in progress. */
async function cancelMidStream(page: Page, gene: string): Promise<void> {
  await expect
    .poll(async () => (await geneLoad(page, gene))?.batches ?? 0, {
      message: `${gene} should have started streaming`,
      timeout: 60_000,
    })
    .toBeGreaterThanOrEqual(1);
  const load = (await geneLoad(page, gene))!;

  expect(load.done, "the download should still be in progress").toBe(false);
  await page
    .locator(`[id="loading-gene-${gene}"] .Toastify__close-button`)
    .click();
  await expect
    .poll(async () => (await geneLoad(page, gene))?.aborted, {
      timeout: 30_000,
    })
    .toBe(true);
  expect((await geneLoad(page, gene))!.loaded).toBeLessThan(
    LARGE_MOLECULES.large,
  );
}

test("@quick standalone molecule viewer streams its big genes", async ({
  page,
}) => {
  await page.goto(
    `/sm-viewer/from-s3?url=${encodeURIComponent(largeMoleculesUrl(server, "a"))}`,
  );
  await expect(page.locator(SELECTORS.smCanvas).first()).toBeVisible({
    timeout: 120_000,
  });
  await expectStreamedDefaults(page, "a");
  await page.waitForFunction(
    () => (window as any).__merfish?.renderComplete != null,
    undefined,
    { timeout: 60_000 },
  );
  expect((await readHooks(page))!.stats!.pointClouds).toBe(3);
});

test("@quick molecules over cells stream the same way", async ({ page }) => {
  await page.goto(
    `/viewer/from-s3?url=${encodeURIComponent(cellsUrl(server, "b"))}`,
  );
  await expect(page.locator(SELECTORS.scCanvas).first()).toBeVisible({
    timeout: 120_000,
  });
  // The cells scene draws them; there is no second scene.
  await expect(page.locator(SELECTORS.smLayerButton)).toBeVisible({
    timeout: 60_000,
  });
  await expect(page.locator(SELECTORS.smCanvas)).toHaveCount(0);
  await expectStreamedDefaults(page, "b");
});

test("@quick a SpatialData store with cells streams its linked molecules", async ({
  page,
}) => {
  await page.goto(
    `/spatialdata-viewer/from-s3?url=${encodeURIComponent(
      storeUrl(server, { cells: true, molecules: "large" }),
    )}`,
  );
  await expect(page.locator(SELECTORS.scCanvas).first()).toBeVisible({
    timeout: 120_000,
  });
  await expect(page.locator(SELECTORS.smLayerButton)).toBeVisible({
    timeout: 60_000,
  });
  await expectStreamedDefaults(page, "store");
});

test("@quick a SpatialData store with only molecules streams them", async ({
  page,
}) => {
  await page.goto(
    `/spatialdata-viewer/from-s3?url=${encodeURIComponent(
      storeUrl(server, { molecules: "large" }),
    )}`,
  );
  await expect(page.locator(SELECTORS.smCanvas).first()).toBeVisible({
    timeout: 120_000,
  });
  // The same scene as the standalone viewer, whose test counts downloads;
  // here a dev-server refresh mid-download has been seen to rebuild it.
  await expectStreamedDefaults(page, "store", { countDownloads: false });
});

test("@quick both panels of a split view stream independently", async ({
  page,
}) => {
  const left = largeMoleculesUrl(server, "a");
  const right = largeMoleculesUrl(server, "b");

  await page.goto(
    `/sm-viewer/from-s3?url=${encodeURIComponent(left)}&splitS3Url=${encodeURIComponent(
      encodeURIComponent(right),
    )}&splitType=sm`,
  );
  await expect(page.locator(SELECTORS.smCanvas)).toHaveCount(2, {
    timeout: 120_000,
  });
  await expectStreamedDefaults(page, "a", { countDownloads: false });
  await expectStreamedDefaults(page, "b", { countDownloads: false });
});

test("@quick closing a gene's toast cancels its download in the molecule viewer", async ({
  page,
}) => {
  await page.goto(
    `/sm-viewer/from-s3?url=${encodeURIComponent(largeMoleculesUrl(server, "c"))}`,
  );
  await expect(page.locator(SELECTORS.smCanvas).first()).toBeVisible({
    timeout: 120_000,
  });
  await cancelMidStream(page, "Big_c");
  // The gene left the selection: two remain, and their loads still finish.
  await expect(page.getByText(/Selected Genes \(2\)/)).toBeVisible();
  await finishedLoad(page, "Huge_c");
  await expect
    .poll(async () => (await readHooks(page))!.stats?.pointClouds)
    .toBe(2);
});

test("@quick closing a gene's toast cancels its download over cells too", async ({
  page,
}) => {
  await page.goto(
    `/viewer/from-s3?url=${encodeURIComponent(cellsUrl(server, "d"))}`,
  );
  await expect(page.locator(SELECTORS.smLayerButton)).toBeVisible({
    timeout: 120_000,
  });
  await cancelMidStream(page, "Big_d");
  await expect(page.getByText(/Selected Genes \(2\)/)).toBeVisible();
  await finishedLoad(page, "Huge_d");
  expect((await geneLoad(page, "Huge_d"))!.aborted).toBe(false);
});
