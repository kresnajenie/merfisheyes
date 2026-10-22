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

/**
 * Already-processed SpatialData dropped on the homepage's single-cell Folder
 * card for a browser preview, in every combination of a table (X as CSR with
 * a CSC layer, or X itself CSC), an OME-Zarr image and a tiled points element,
 * and every way of laying them out:
 *
 *   store                  the dropped folder is one store holding them all
 *   store in a folder      a folder holding that store
 *   separate, in a folder  a folder holding table.zarr, image.zarr, points.zarr
 *   separate, together     those three dragged in at once, as three items
 *
 * Dropped through the drop zone (dropFolders), not the folder picker, so the
 * walk of the dropped entries is part of what is tested. Every check is soft:
 * a failing case reports all that is wrong with it. Tagged @quick.
 */

type Layout =
  | "store"
  | "store in a folder"
  | "separate, in a folder"
  | "separate, together";

const LAYOUTS: Layout[] = [
  "store",
  "store in a folder",
  "separate, in a folder",
  "separate, together",
];

const SPECS: StoreSpec[] = [];

for (const table of [undefined, "csr", "csc_x"] as const) {
  for (const image of [false, true]) {
    for (const points of [false, true]) {
      if (table || image || points) SPECS.push({ table, image, points });
    }
  }
}

const label = (spec: StoreSpec) =>
  [
    spec.table &&
      `table (${spec.table === "csr" ? "CSR X + CSC layer" : "CSC X"})`,
    spec.image && "image",
    spec.points && "points",
  ]
    .filter(Boolean)
    .join(" + ");

const tmp = (prefix: string) => fs.mkdtempSync(path.join(os.tmpdir(), prefix));

/** The folders to drop, as one item or several. */
function folders(layout: Layout, spec: StoreSpec): string[] {
  if (layout === "store") {
    const store = path.join(tmp("store-"), "sample.zarr");

    writeStore(store, spec);

    return [store];
  }
  if (layout === "store in a folder") {
    const folder = path.join(tmp("bundle-"), "sample");

    writeStore(path.join(folder, "store.zarr"), spec);

    return [folder];
  }
  const folder = makeSeparateFolders(spec);

  return layout === "separate, in a folder"
    ? [folder]
    : fs.readdirSync(folder).map((name) => path.join(folder, name));
}

// Compile the pages once, so no case spends its time budget on the dev server
test.beforeAll(async ({ browser }) => {
  test.setTimeout(300_000);
  const page = await browser.newPage();

  for (const route of [
    "/",
    "/viewer/from-local",
    "/sm-viewer/from-local",
    "/image-viewer/from-local",
  ]) {
    await page.goto(route, { waitUntil: "load", timeout: 240_000 });
  }
  await page.close();
});

const TINY_CELLS = 30;
const TINY_GENES = 8;
const POINTS_GENES = 4;

for (const layout of LAYOUTS) {
  for (const spec of SPECS) {
    test(`@quick drop matrix: ${layout}: ${label(spec)}`, async ({ page }) => {
      const errors: string[] = [];

      page.on("pageerror", (e) => errors.push(e.message));
      await page.goto("/", { waitUntil: "domcontentloaded" });
      await dropFolders(
        page,
        page.locator("#sc-file-input-folder").locator(".."),
        folders(layout, spec),
      );

      // Cells are the scene when there are any; then molecules; then the image
      const viewer = spec.table
        ? /\/viewer\/from-local/
        : spec.points
          ? /\/sm-viewer\/from-local/
          : /\/image-viewer\/from-local/;
      const arrived = await page
        .waitForURL(viewer, { timeout: 30_000 })
        .then(() => true)
        .catch(() => false);

      expect.soft(arrived, `opens ${viewer} (at ${page.url()})`).toBe(true);
      if (!arrived) {
        expect.soft(errors, "page errors").toEqual([]);

        return;
      }

      if (spec.table) {
        await test.step("reads the table, and a gene from it", async () => {
          await page
            .waitForFunction(
              () => (window as any).__merfish?.renderComplete != null,
              undefined,
              { timeout: 30_000 },
            )
            .catch(() => undefined);
          const stats = (await readHooks(page))?.stats;

          expect.soft(stats?.pointCount, "cells").toBe(TINY_CELLS);
          expect.soft(stats?.geneCount, "genes").toBe(TINY_GENES);

          const search = page.locator(SELECTORS.geneSearchInput);

          if (!(await search.isVisible())) {
            await page.locator(SELECTORS.geneModeButton).first().click();
          }
          await search.fill("gene1");
          await page.getByRole("radio", { name: "gene1", exact: true }).click();
          await expect
            .soft(page.locator(SELECTORS.selectedGeneBadge))
            .toHaveAttribute("data-gene", "gene1", { timeout: 15_000 });
        });
      }

      if (spec.points) {
        await test.step("shows the points' molecules", async () => {
          if (spec.table) {
            await expect
              .soft(
                page
                  .getByTestId("points-status")
                  .or(page.locator(SELECTORS.smLayerButton))
                  .first(),
                "a molecule layer over the cells",
              )
              .toBeVisible({ timeout: 15_000 });
          } else {
            await expect
              .soft(page.locator(SELECTORS.smCanvas).first())
              .toBeVisible({ timeout: 15_000 });
            await expect.soft
              .poll(async () => (await readHooks(page))?.stats?.geneCount, {
                message: "the points' genes",
                timeout: 15_000,
              })
              .toBe(POINTS_GENES);
          }
        });
      } else {
        await expect
          .soft(page.getByTestId("points-status"), "no molecule layer")
          .toHaveCount(0);
        await expect
          .soft(page.locator(SELECTORS.smLayerButton), "no molecule layer")
          .toHaveCount(0);
      }

      await test.step("shows the image when there is one", async () => {
        if (!spec.image) {
          await expect.soft(page.locator(SELECTORS.imageButton)).toHaveCount(0);
        } else if (spec.table || spec.points) {
          await expect
            .soft(page.locator(SELECTORS.imageButton))
            .toBeVisible({ timeout: 15_000 });
        } else {
          await expect.soft
            .poll(() => litFraction(page, SELECTORS.imageCanvas), {
              message: "the image drawn",
              timeout: 15_000,
            })
            .toBeGreaterThan(0.2);
        }
      });

      await expect
        .soft(page.getByText(/Skipped:/), "nothing skipped")
        .toHaveCount(0);
      expect.soft(errors, "page errors").toEqual([]);
    });
  }
}
