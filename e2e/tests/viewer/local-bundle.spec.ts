import {
  test,
  expect,
  readHooks,
  litFraction,
  makeBundle,
  makeImageFolder,
  SELECTORS,
  type BundleParts,
} from "../../fixtures/merfish-fixtures";

/**
 * Folders of already-processed data dropped on the homepage for a browser
 * preview: an OME-Zarr image on its own, and cells, molecules and an image
 * together, overlaid. Nothing is converted and nothing is uploaded, so these
 * run without a database or S3. Tagged @quick.
 */

async function drop(
  page: import("@playwright/test").Page,
  folder: string,
): Promise<void> {
  await page.goto("/", { waitUntil: "domcontentloaded" });
  await page.locator("#sc-file-input-folder").setInputFiles(folder);
}

async function cellsRendered(page: import("@playwright/test").Page) {
  await page.waitForURL(/\/viewer\/from-local/, { timeout: 120_000 });
  await expect(page.locator(SELECTORS.scCanvas).first()).toBeVisible();
  await page.waitForFunction(
    () => (window as any).__merfish?.renderComplete != null,
    undefined,
    { timeout: 120_000 },
  );

  return (await readHooks(page))!.stats!;
}

test("@quick local preview: an OME-Zarr image folder on its own", async ({
  page,
}) => {
  await drop(page, makeImageFolder());
  await page.waitForURL(/\/image-viewer\/from-local/, { timeout: 120_000 });
  await expect(page.locator(SELECTORS.imageButton)).toBeVisible();
  await expect
    .poll(() => litFraction(page, SELECTORS.imageCanvas), {
      message: "the image should be drawn from the dropped files",
      timeout: 60_000,
    })
    .toBeGreaterThan(0.2);
});

/** merfisheyes outputs side by side, and a SpatialData store. */
const WITH_CELLS: { name: string; parts: BundleParts; cells: number }[] = [
  {
    name: "cells + molecules + image (merfisheyes outputs)",
    parts: { cells: "chunked", molecules: true, images: ["dapi.zarr"] },
    cells: 300,
  },
  {
    name: "cells + image (merfisheyes outputs)",
    parts: { cells: "chunked", images: ["dapi.zarr"] },
    cells: 300,
  },
  {
    name: "cells + molecules (merfisheyes outputs)",
    parts: { cells: "chunked", molecules: true },
    cells: 300,
  },
  {
    name: "SpatialData store (cells + its image) + molecules",
    parts: { cells: "spatialdata", molecules: true },
    cells: 30,
  },
  {
    name: "SpatialData store alone",
    parts: { cells: "spatialdata" },
    cells: 30,
  },
];

for (const { name, parts, cells } of WITH_CELLS) {
  test(`@quick local preview: ${name}`, async ({ page }) => {
    await drop(page, makeBundle(parts));
    const stats = await cellsRendered(page);

    expect(stats.pointCount, "cells").toBe(cells);

    // A store carries its own image; the other form brings one beside it
    const hasImage = parts.cells === "spatialdata" || !!parts.images?.length;

    if (hasImage) {
      await expect(page.locator(SELECTORS.imageButton)).toBeVisible();
    } else {
      await expect(page.locator(SELECTORS.imageButton)).toHaveCount(0);
    }
    if (parts.cells === "chunked") {
      // The image is stretched under all the cells, so it fills the view;
      // the cells on their own are a sprinkle of dots
      const lit = expect.poll(() => litFraction(page, SELECTORS.scCanvas), {
        timeout: 60_000,
      });

      await (hasImage ? lit.toBeGreaterThan(0.08) : lit.toBeLessThan(0.05));
    }

    if (parts.molecules) {
      await expect(page.locator(SELECTORS.smLayerButton)).toBeVisible();
    } else {
      await expect(page.locator(SELECTORS.smLayerButton)).toHaveCount(0);
    }
  });
}

test("@quick local preview: molecules + image, no cells", async ({ page }) => {
  await drop(page, makeBundle({ molecules: true, images: ["dapi.zarr"] }));
  await page.waitForURL(/\/sm-viewer\/from-local/, { timeout: 120_000 });
  await expect(page.locator(SELECTORS.smCanvas).first()).toBeVisible();
  await expect(page.locator(SELECTORS.imageButton)).toBeVisible();
  await expect
    .poll(() => litFraction(page, SELECTORS.smCanvas), { timeout: 60_000 })
    .toBeGreaterThan(0.2);
});

test("@quick local preview: the first of several images is used, and says so", async ({
  page,
}) => {
  await drop(
    page,
    makeBundle({ cells: "chunked", images: ["b.zarr", "a.zarr"] }),
  );
  await cellsRendered(page);
  await expect(page.getByText(/Skipped: image b\.zarr/)).toBeVisible();
  await expect(page.locator(SELECTORS.imageButton)).toBeVisible();
});

test("@quick local preview: raw cells beside processed parts show the cells only", async ({
  page,
}) => {
  await drop(
    page,
    makeBundle({ cells: "h5ad", molecules: true, images: ["dapi.zarr"] }),
  );
  await expect(
    page.getByText(/The browser preview shows the cells only/),
  ).toBeVisible({ timeout: 60_000 });
  const stats = await cellsRendered(page);

  expect(stats.pointCount).toBe(300);
  await expect(page.locator(SELECTORS.imageButton)).toHaveCount(0);
  await expect(page.locator(SELECTORS.smLayerButton)).toHaveCount(0);
});

test("@quick local preview: an image is not drawn under normalized cells, and says so", async ({
  page,
}) => {
  await drop(page, makeBundle({ cells: "normalized", images: ["dapi.zarr"] }));
  await cellsRendered(page);
  await expect(page.getByText(/The image is not drawn/)).toBeVisible();
  await expect(page.locator(SELECTORS.imageButton)).toHaveCount(0);
});
