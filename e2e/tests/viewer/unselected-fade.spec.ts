import {
  expect,
  SELECTORS,
  startSpatialDataServer,
  storeUrl,
  test,
} from "../../fixtures/merfish-fixtures";
import { type FixtureServer } from "../../helpers/fixture-server";

/**
 * Cells outside the selected cell types fade as the view zooms out
 * (lib/webgl/shaders.ts). How far in they are opaque again is "Unselected
 * fade" in Advanced settings, a uniform of the cell cloud. Tagged @quick.
 */

let server: FixtureServer;

test.beforeAll(async () => {
  server = await startSpatialDataServer();
});
test.afterAll(() => server.close());

test("@quick Unselected fade is set in Advanced settings and kept in the link", async ({
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
  // The gear at the top right, found by its icon
  const advanced = page
    .locator("button")
    .filter({ has: page.locator('path[d^="M9.594 3.94"]') });

  await advanced.click();

  const slider = page.getByRole("slider", {
    name: "Unselected fade (zoom ×, 1 = off)",
  });

  // The default: opaque again four times closer than the whole dataset
  await expect(slider).toHaveValue("4");
  // The least turns the fade off
  await slider.focus();
  await page.keyboard.press("Home");
  await expect(slider).toHaveValue("1");
  // In the link's `v` (base64url JSON), with the other advanced settings
  await expect
    .poll(() => {
      const v = new URL(page.url()).searchParams.get("v");

      return v ? Buffer.from(v, "base64url").toString() : "";
    })
    .toContain('"unselectedFadeZoom":1');

  // The link opens with it
  await page.reload();
  await expect(page.locator(SELECTORS.scCanvas).first()).toBeVisible({
    timeout: 120_000,
  });
  await advanced.click();
  await expect(slider).toHaveValue("1");
});
