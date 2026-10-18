import {
  expect,
  SELECTORS,
  startSpatialDataServer,
  storeUrl,
  test,
} from "../../fixtures/merfish-fixtures";
import { firstContent, layer, settledAfter } from "../../helpers/budgets";
import { type FixtureServer } from "../../helpers/fixture-server";

/**
 * A store holding every element kind, opened through the one dataset
 * loader: each kind is found by the element registry, reaches its layer,
 * draws and settles, and its controls appear where the registry puts them.
 * Tagged @quick: the tiny fixtures, served by a local stand-in for the
 * bucket, so no database or S3 is needed.
 */

let server: FixtureServer;

test.beforeAll(async () => {
  server = await startSpatialDataServer();
});
test.afterAll(() => server.close());

test("@quick every element kind of a store reaches its layer", async ({
  page,
}) => {
  await page.goto(
    `/spatialdata-viewer/from-s3?url=${encodeURIComponent(
      storeUrl(server, {
        cells: true,
        points: true,
        shapes: true,
        image: true,
      }),
    )}`,
  );
  await expect(page.locator(SELECTORS.scCanvas).first()).toBeVisible({
    timeout: 120_000,
  });
  await firstContent(page, "table");

  // The tiny store is about 100 units across, under the width at which
  // polygons are drawn, so every kind draws in the opening view
  for (const kind of ["image", "points", "shapes"]) {
    await test.step(`${kind}: first content, then a settled view`, async () => {
      const at = await firstContent(page, kind);

      await settledAfter(page, kind, at - 1);
      expect((await layer(page, kind))!.drawn, `${kind} drawn`).toBeGreaterThan(
        0,
      );
    });
  }

  await test.step("each layer's controls are where the registry puts them", async () => {
    // The image controls and the molecule status: panels inside the scene
    await expect(page.locator(SELECTORS.imageButton)).toBeVisible();
    await expect(page.getByTestId("points-status").first()).toBeVisible();
    // The cell shapes: a button in the viewer's rail
    await expect(page.getByTestId("shapes-button")).toBeVisible();
    await page.getByTestId("shapes-button").click();
    await expect(page.getByLabel("Show cell shapes")).toBeChecked();
  });
});
