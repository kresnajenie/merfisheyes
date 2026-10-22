import fs from "node:fs";
import os from "node:os";
import path from "node:path";

import {
  test,
  expect,
  readHooks,
  litFraction,
  dropFolders,
  makeSeparateFolders,
  writeStore,
  SELECTORS,
  type StoreSpec,
} from "../../fixtures/merfish-fixtures";
import { signInDev } from "../../helpers/auth";

/**
 * Already-processed SpatialData dropped on the homepage, on the two paths
 * that keep it: the browser preview's "Upload & Save" (every part uploaded
 * as it is, then reopened from the saved link), and the "Upload & process on
 * server" card, which should take processed data as it is rather than look
 * for raw input. The same layouts as e2e/tests/viewer/drop-matrix.spec.ts.
 * Tagged @upload (needs Postgres + S3).
 */

const UPLOAD_ENABLED = process.env.E2E_UPLOAD === "1";

type Page = import("@playwright/test").Page;

const tmp = (prefix: string) => fs.mkdtempSync(path.join(os.tmpdir(), prefix));

interface Case {
  name: string;
  spec: StoreSpec;
  /** The folders dropped, as one item or several. */
  folders: () => string[];
}

const inStore = (spec: StoreSpec) => () => {
  const store = path.join(tmp("store-"), "sample.zarr");

  writeStore(store, spec);

  return [store];
};
const separate = (spec: StoreSpec) => () => [makeSeparateFolders(spec)];
/** The same folders, dragged in at once as several items. */
const together = (spec: StoreSpec) => () => {
  const folder = makeSeparateFolders(spec);

  return fs.readdirSync(folder).map((name) => path.join(folder, name));
};

const CASES: Case[] = [
  {
    name: "store: table (CSR X + CSC layer) + image + points",
    spec: { table: "csr", image: true, points: true },
    folders: inStore({ table: "csr", image: true, points: true }),
  },
  {
    name: "store: table (CSC X) + image",
    spec: { table: "csc_x", image: true },
    folders: inStore({ table: "csc_x", image: true }),
  },
  {
    name: "store: image + points, no table",
    spec: { image: true, points: true },
    folders: inStore({ image: true, points: true }),
  },
  {
    name: "store: an image only",
    spec: { image: true },
    folders: inStore({ image: true }),
  },
  {
    name: "separate, in a folder: table + image + points",
    spec: { table: "csr", image: true, points: true },
    folders: separate({ table: "csr", image: true, points: true }),
  },
  {
    name: "separate, dragged together: table + image + points",
    spec: { table: "csr", image: true, points: true },
    folders: together({ table: "csr", image: true, points: true }),
  },
];

const dropZone = (page: Page) =>
  page.locator("#sc-file-input-folder").locator("..");

/** What a reopened dataset must show: its cells, image and molecules. */
async function showsParts(page: Page, spec: StoreSpec): Promise<void> {
  if (spec.table) {
    await expect.soft
      .poll(async () => (await readHooks(page))?.stats?.pointCount, {
        message: "cells",
        timeout: 120_000,
      })
      .toBe(30);
  } else if (spec.points) {
    await expect
      .soft(page.locator(SELECTORS.smCanvas).first())
      .toBeVisible({ timeout: 120_000 });
  }
  if (spec.image && !spec.table && !spec.points) {
    // An image alone is its own scene
    await expect.soft
      .poll(() => litFraction(page, SELECTORS.imageCanvas), {
        message: "the image drawn",
        timeout: 60_000,
      })
      .toBeGreaterThan(0.2);
  } else if (spec.image) {
    await expect
      .soft(page.locator(SELECTORS.imageButton))
      .toBeVisible({ timeout: 60_000 });
  }
  if (spec.points && spec.table) {
    await expect
      .soft(
        page
          .getByTestId("points-status")
          .or(page.locator(SELECTORS.smLayerButton))
          .first(),
        "a molecule layer over the cells",
      )
      .toBeVisible({ timeout: 60_000 });
  }
}

/** Name it in the save dialog, upload, and open what was saved. */
async function saveAndReopen(
  page: Page,
  title: string,
  spec: StoreSpec,
): Promise<void> {
  await page.getByLabel(/Project name|Dataset name/).fill(title);
  await page.getByRole("button", { name: "Upload", exact: true }).click();
  await expect(page.getByTestId("bundle-saved")).toBeVisible({
    timeout: 180_000,
  });
  await page.getByText("Open the saved dataset").click();
  await page.waitForURL(/\/ds_/, { timeout: 60_000 });
  await showsParts(page, spec);
}

test.describe("@upload drop matrix: Upload & Save from the preview", () => {
  test.skip(
    !UPLOAD_ENABLED,
    "Set E2E_UPLOAD=1 (and start docker infra) to run upload tests.",
  );

  for (const { name, spec, folders } of CASES) {
    test(name, async ({ page }) => {
      const title = `e2e-drop-${Date.now()}`;

      await signInDev(page, `${title}@example.com`);
      await page.goto("/", { waitUntil: "domcontentloaded" });
      await dropFolders(page, dropZone(page), folders());
      await page.waitForURL(/\/(viewer|sm-viewer|image-viewer)\/from-local/, {
        timeout: 120_000,
      });

      await test.step("the preview shows every part", () =>
        showsParts(page, spec));

      await test.step("saved as it is, and reopened from the link", async () => {
        await page.locator(SELECTORS.uploadSaveButton).first().click();
        await saveAndReopen(page, title, spec);
      });
    });
  }
});

test.describe("@upload drop matrix: Upload & process on server", () => {
  test.skip(
    !UPLOAD_ENABLED,
    "Set E2E_UPLOAD=1 (and start docker infra) to run upload tests.",
  );

  for (const { name, spec, folders } of CASES) {
    test(name, async ({ page }) => {
      const title = `e2e-server-drop-${Date.now()}`;

      await signInDev(page, `${title}@example.com`);
      await page.goto("/", { waitUntil: "domcontentloaded" });
      await page
        .getByRole("radio", { name: "Upload & process on server" })
        .click();
      await dropFolders(page, dropZone(page), folders());

      // Processed data needs no processing: it is taken as it is, never
      // turned away for not being raw input
      const dialog = page.getByLabel(/Project name|Dataset name/);
      const refusal = page.getByText(/No single-cell data found/);

      await expect(dialog.or(refusal).first()).toBeVisible({ timeout: 30_000 });
      expect(await refusal.isVisible(), "turned away").toBe(false);

      await test.step("saved as it is, and opened from the link", () =>
        saveAndReopen(page, title, spec));
    });
  }
});
