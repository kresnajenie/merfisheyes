import type { Page } from "@playwright/test";

import {
  test,
  expect,
  startSpatialDataServer,
  storeUrl,
  largeMoleculesUrl,
  SELECTORS,
} from "../../fixtures/merfish-fixtures";
import { type FixtureServer } from "../../helpers/fixture-server";

/**
 * The molecule opacity slider: both molecule layers — a store's tiled points
 * element and a molecule folder's per-gene clouds — dim when it goes down and
 * come back when it goes up. Measured on the drawn canvas. Tagged @quick.
 */

let server: FixtureServer;

test.beforeAll(async () => {
  server = await startSpatialDataServer();
});
test.afterAll(() => server.close());

/**
 * How much light the middle of the canvas gives off: the mean of each
 * pixel's brightest channel. Only the middle, where the molecules are: a
 * screenshot of the canvas also holds the controls and banners drawn over it.
 */
async function drawnBrightness(page: Page): Promise<number> {
  const shot = await page.locator(SELECTORS.smCanvas).first().screenshot();

  return page.evaluate(async (png) => {
    const image = await createImageBitmap(
      await (await fetch(`data:image/png;base64,${png}`)).blob(),
    );
    const context = new OffscreenCanvas(image.width, image.height).getContext("2d")!;

    context.drawImage(image, 0, 0);
    const { data } = context.getImageData(
      Math.floor(image.width * 0.33),
      Math.floor(image.height * 0.33),
      Math.floor(image.width * 0.3),
      Math.floor(image.height * 0.25),
    );
    let sum = 0;

    for (let i = 0; i < data.length; i += 4) {
      sum += Math.max(data[i], data[i + 1], data[i + 2]);
    }

    return sum / (data.length / 4);
  }, shot.toString("base64"));
}

async function expectOpacityToDim(page: Page): Promise<void> {
  await expect(page.locator(SELECTORS.smCanvas).first()).toBeVisible({
    timeout: 120_000,
  });
  // Something drawn first
  await expect.poll(() => drawnBrightness(page), { timeout: 60_000 }).toBeGreaterThan(0.5);
  const opaque = await drawnBrightness(page);
  const thumb = page.getByRole("slider", { name: "Molecule opacity" });

  await thumb.focus();
  await page.keyboard.press("Home");
  await expect
    .poll(() => drawnBrightness(page), { timeout: 15_000 })
    .toBeLessThan(opaque * 0.3);

  await page.keyboard.press("End");
  await expect
    .poll(() => drawnBrightness(page), { timeout: 15_000 })
    .toBeGreaterThan(opaque * 0.9);
}

test("@quick a store's tiled molecules fade with the opacity slider", async ({ page }) => {
  await page.goto(
    `/spatialdata-viewer/from-s3?url=${encodeURIComponent(storeUrl(server, { points: true }))}`,
  );
  await expectOpacityToDim(page);
});

test("@quick a molecule folder's genes fade with the opacity slider", async ({ page }) => {
  await page.goto(
    `/sm-viewer/from-s3?url=${encodeURIComponent(largeMoleculesUrl(server, "opacity"))}`,
  );
  await expectOpacityToDim(page);
});
