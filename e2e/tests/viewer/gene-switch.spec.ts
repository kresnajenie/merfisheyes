import {
  expect,
  SELECTORS,
  startSpatialDataServer,
  storeUrl,
  test,
} from "../../fixtures/merfish-fixtures";
import { type FixtureServer } from "../../helpers/fixture-server";

/**
 * Picking a second gene while the first is still downloading
 * (components/three-scene.tsx): each gene's "Loading expression" toast goes
 * when its own read ends. Tagged @quick: the tiny store, its expression
 * reads held back so that the genes' reads overlap.
 */

let server: FixtureServer;

test.beforeAll(async () => {
  server = await startSpatialDataServer();
});
test.afterAll(() => server.close());

test("@quick no loading toast is left after switching genes mid-download", async ({
  page,
}) => {
  await page.goto(
    `/spatialdata-viewer/from-s3?url=${encodeURIComponent(
      storeUrl(server, { cells: true }),
    )}`,
  );
  await expect(page.locator(SELECTORS.scCanvas).first()).toBeVisible({
    timeout: 120_000,
  });
  // Every read of the expression matrix takes a second and a half
  await page.route(/\/tables\/[^/]+\/(X|layers)\//, async (route) => {
    await new Promise((resolve) => setTimeout(resolve, 1500));
    await route.continue();
  });
  await page.locator(SELECTORS.geneModeButton).first().click();
  await expect(page.locator(SELECTORS.geneSearchInput)).toBeVisible();

  const genes = page.getByRole("radio");
  const loading = page.locator(".Toastify__toast", {
    hasText: "Loading expression",
  });

  await genes.nth(0).click({ force: true });
  await expect(loading).toHaveCount(1);
  await genes.nth(1).click({ force: true });
  await genes.nth(2).click({ force: true });

  await expect(page.locator(SELECTORS.selectedGeneBadge)).toBeVisible();
  await expect(loading).toHaveCount(0, { timeout: 60_000 });
});
