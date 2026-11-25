import {
  expect,
  SELECTORS,
  startSpatialDataServer,
  storeUrl,
  test,
} from "../../fixtures/merfish-fixtures";
import { type FixtureServer } from "../../helpers/fixture-server";

/**
 * The download bar (components/load-meter-bar.tsx, lib/net/load-meter.ts):
 * while data is being read it says how much has arrived of how much was
 * announced. Tagged @quick: the tiny store, served slowly by a local stand-in
 * for the bucket so that a read lasts long enough to be seen.
 */

let server: FixtureServer;

test.beforeAll(async () => {
  server = await startSpatialDataServer({ bytesPerSecond: 60_000 });
});
test.afterAll(() => server.close());

test.describe("@quick the download bar", () => {
  test("shows megabytes received of megabytes announced, then goes", async ({
    page,
  }) => {
    await page.goto(
      `/spatialdata-viewer/from-s3?url=${encodeURIComponent(
        storeUrl(server, { cells: true, image: true }),
      )}`,
    );

    // "412 KB / 1.3 MB · 3 of 9 files"
    const amounts =
      /\d+(\.\d+)? (KB|MB|GB) \/ \d+(\.\d+)? (KB|MB|GB) · \d+ of \d+ files?/;
    const meter = page.getByTestId("load-meter");

    // While the store opens, the loading screen says it and the bar stays away
    await expect(page.getByTestId("load-meter-text")).toContainText(amounts, {
      timeout: 60_000,
    });
    await expect(meter).toBeHidden();

    // In the viewer, what is still being read shows on the bar
    await expect(page.locator(SELECTORS.scCanvas).first()).toBeVisible({
      timeout: 120_000,
    });
    await expect(meter).toContainText(amounts);
    await expect(meter.getByRole("progressbar")).toBeVisible();
    // Once nothing is in flight it is gone
    await expect(meter).toBeHidden({ timeout: 120_000 });
  });
});
