import type { Page } from "@playwright/test";

import {
  test,
  expect,
  startSpatialDataServer,
  storeUrl,
  cellsUrl,
  SELECTORS,
} from "../../fixtures/merfish-fixtures";
import { type FixtureServer } from "../../helpers/fixture-server";

/**
 * Split-screen sync between the viewers as they are now: a cell folder, a
 * SpatialData store, an image. The button is offered when both panels have
 * something in common (lib/stores/splitScreenStore.ts syncCapabilities) and
 * what it offers works. It used to compare the viewers instead, so a cell
 * folder beside a store showed no button, and two stores showed one that
 * synced nothing. Tagged @quick: committed tiny fixtures on a local stand-in
 * for the bucket.
 */

let server: FixtureServer;

test.beforeAll(async () => {
  server = await startSpatialDataServer();
});
test.afterAll(() => server.close());

/** The left page with the right panel opened from a URL, as a shared link does. */
function split(left: string, right: string, type: string, sync = false) {
  const url = new URL(left, "http://x");

  url.searchParams.set("splitS3Url", encodeURIComponent(right));
  url.searchParams.set("splitType", type);
  if (sync) url.searchParams.set("sync", "1");

  return url.pathname + url.search;
}

const store = (parts: Parameters<typeof storeUrl>[1]) =>
  `/spatialdata-viewer/from-s3?url=${encodeURIComponent(storeUrl(server, parts))}`;

const syncButton = (page: Page) =>
  page.getByRole("button", { name: /^(Enable|Disable) sync$/ });

/** Choose a gene in the left panel's gene list. */
async function chooseGeneOnLeft(page: Page, gene: string) {
  await page.locator(SELECTORS.geneModeButton).first().click();
  const search = page.locator(SELECTORS.geneSearchInput).first();

  await expect(search).toBeVisible();
  await search.fill(gene);
  await page.getByRole("radio", { name: gene, exact: true }).first().click();
}

const badges = (page: Page, gene: string) =>
  page.locator(`${SELECTORS.selectedGeneBadge}[data-gene="${gene}"]`);

test("@quick a cell folder beside a SpatialData store with cells offers sync", async ({
  page,
}) => {
  await page.goto(
    split(
      `/viewer/from-s3?url=${encodeURIComponent(cellsUrl(server, "split-sync"))}`,
      storeUrl(server, { cells: true }),
      "spatialdata",
    ),
  );

  await expect(syncButton(page)).toBeVisible({ timeout: 120_000 });
  await expect(syncButton(page)).toHaveAttribute("aria-pressed", "false");
});

test("@quick two SpatialData stores sync a gene chosen on the left", async ({
  page,
}) => {
  await page.goto(
    split(
      store({ cells: true }),
      storeUrl(server, { cells: true }),
      "spatialdata",
    ),
  );

  const button = syncButton(page);

  await expect(button).toBeVisible({ timeout: 120_000 });
  await button.click();
  await expect(button).toHaveAttribute("aria-pressed", "true");

  await chooseGeneOnLeft(page, "gene2");
  // The left panel's legend and the right one's both show it
  await expect(badges(page, "gene2")).toHaveCount(2, { timeout: 30_000 });
  await expect(page).toHaveURL(/[?&]sync=1/);
});

test("@quick a shared link with sync=1 opens with sync on, and it syncs", async ({
  page,
}) => {
  await page.goto(
    split(
      store({ cells: true }),
      storeUrl(server, { cells: true }),
      "spatialdata",
      true,
    ),
  );

  const button = syncButton(page);

  await expect(button).toBeVisible({ timeout: 120_000 });
  await expect(button).toHaveAttribute("aria-pressed", "true");
  // A restored sync lets both panels apply their own link first: the hooks
  // only propagate 3 s after both have loaded (useSyncVisualization)
  await page.waitForTimeout(4_000);
  await chooseGeneOnLeft(page, "gene1");
  await expect(badges(page, "gene1")).toHaveCount(2, { timeout: 30_000 });
});

test("@quick two images offer no sync, as nothing in an image syncs yet", async ({
  page,
}) => {
  await page.goto(
    split(
      store({ image: true }),
      storeUrl(server, { image: true }),
      "spatialdata",
    ),
  );

  // Both panels are up before judging the button's absence
  await expect(page.locator(SELECTORS.imageCanvas)).toHaveCount(2, {
    timeout: 120_000,
  });
  await expect(syncButton(page)).toHaveCount(0);
});
