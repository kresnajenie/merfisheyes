import fs from "node:fs";
import os from "node:os";
import path from "node:path";

import {
  test,
  expect,
  readHooks,
  REPO_ROOT,
  SELECTORS,
} from "../../fixtures/merfish-fixtures";

/**
 * A SpatialData store whose table carries its DE stats (varm / uns, written
 * by scripts/de_stats.py) opens with the DEG panel; the same store without
 * them has none to offer — the viewer reads a store's table a gene at a time
 * and cannot compute them itself. Dropped as a folder for a browser preview.
 * Tagged @quick.
 */

const ZARR = path.join(REPO_ROOT, "tests", "data", "tiny", "zarr");

function dropFolder(store: string): string {
  const folder = path.join(
    fs.mkdtempSync(path.join(os.tmpdir(), "deg-")),
    "store.zarr",
  );

  fs.cpSync(path.join(ZARR, store), folder, { recursive: true });

  return folder;
}

for (const { store, de } of [
  { store: "spatialdata_v3_de.zarr", de: true },
  { store: "spatialdata_v3.zarr", de: false },
]) {
  test(`@quick DEG from the table: ${store}`, async ({ page }) => {
    await page.goto("/", { waitUntil: "domcontentloaded" });
    await page
      .locator("#sc-file-input-folder")
      .setInputFiles(dropFolder(store));
    await page.waitForURL(/\/viewer\/from-local/, { timeout: 120_000 });
    await expect
      .poll(async () => (await readHooks(page))?.stats?.pointCount, {
        timeout: 60_000,
      })
      .toBe(30);

    const button = page.locator(SELECTORS.degButton);

    if (!de) {
      await expect(button).toHaveCount(0);

      return;
    }
    await expect(button).toBeVisible({ timeout: 30_000 });
    await button.click();
    const panel = page.locator(SELECTORS.degPanel);

    await expect(panel).toBeVisible();
    const auto = panel.locator('button:has-text("Auto")').first();

    if (await auto.isVisible()) await auto.click();
    await panel.getByLabel("Target celltype", { exact: false }).click();
    await page.getByRole("option", { name: "neuron" }).click();
    await expect(
      page.locator(`${SELECTORS.degPanel} [data-testid="deg-row"]`).first(),
    ).toBeVisible({ timeout: 10_000 });
  });
}
