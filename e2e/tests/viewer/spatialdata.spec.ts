import {
  test,
  expect,
  readHooks,
  litFraction,
  startSpatialDataServer,
  storeUrl,
  imageUrl,
  SELECTORS,
  type StoreLayout,
  type StoreParts,
} from "../../fixtures/merfish-fixtures";
import { type FixtureServer } from "../../helpers/fixture-server";

/**
 * SpatialData and image viewers, opened by URL the way a public bucket is
 * (/spatialdata-viewer/from-s3 and /image-viewer/from-s3). Tagged @quick: the
 * stores are the committed tiny fixtures, served by a local stand-in for the
 * bucket, so no database or S3 is needed.
 */

const CELLS = 30;
const GENES = 8;

const SCENES = {
  cells: SELECTORS.scCanvas,
  molecules: SELECTORS.smCanvas,
  image: SELECTORS.imageCanvas,
};

/** A store may hold any of cells, molecules and an image. */
const STORES: { parts: StoreParts; layout?: StoreLayout }[] = [
  { parts: { cells: true } },
  { parts: { cells: true, image: true } },
  { parts: { cells: true, molecules: true } },
  { parts: { cells: true, molecules: true, image: true } },
  { parts: { molecules: true } },
  { parts: { molecules: true, image: true } },
  { parts: { image: true } },
  // The table as anndata >= 0.12 writes it: sharded arrays, nullable strings
  { parts: { cells: true, image: true }, layout: "spatialdata_v3_sharded" },
];

let server: FixtureServer;

test.beforeAll(async () => {
  server = await startSpatialDataServer();
});
test.afterAll(() => server.close());
test.beforeEach(() => {
  server.requests.length = 0;
});

const served = (pattern: RegExp) =>
  server.requests.some((r) => r.status < 300 && pattern.test(r.path));
const names = (parts: StoreParts) =>
  Object.keys(parts).join(" + ") || "nothing";

for (const { parts, layout } of STORES) {
  const label = `${names(parts)}${layout ? ` (${layout})` : ""}`;

  test(`@quick SpatialData store with ${label}`, async ({ page }) => {
    await page.goto(
      `/spatialdata-viewer/from-s3?url=${encodeURIComponent(storeUrl(server, parts, layout))}`,
    );

    // Cells are the scene when there are any; then molecules; then the image
    const scene = parts.cells
      ? "cells"
      : parts.molecules
        ? "molecules"
        : "image";

    await test.step(`shows the ${scene} scene, and only that one`, async () => {
      await expect(page.locator(SCENES[scene]).first()).toBeVisible({
        timeout: 120_000,
      });
      for (const [other, selector] of Object.entries(SCENES)) {
        if (other !== scene)
          await expect(page.locator(selector)).toHaveCount(0);
      }
    });

    if (parts.cells) {
      await test.step("reads the table", async () => {
        await page.waitForFunction(
          () => (window as any).__merfish?.renderComplete != null,
          undefined,
          { timeout: 120_000 },
        );
        const stats = (await readHooks(page))!.stats!;

        expect(stats.pointCount, "cells").toBe(CELLS);
        expect(stats.geneCount, "genes").toBe(GENES);
      });
    }

    await test.step("draws the image only when the store has one", async () => {
      const button = page.locator(SELECTORS.imageButton);

      if (parts.image) {
        await expect(button).toBeVisible();
        await expect
          .poll(() => served(/\/images\/dapi\/\d\/c\//), {
            message: "an image tile should be fetched",
          })
          .toBe(true);
      } else {
        await expect(button).toHaveCount(0);
        expect(served(/\/images\//), "no image should be read").toBe(false);
      }
    });

    if (parts.cells) {
      await test.step("draws molecules over the cells only when linked", async () => {
        // The scene looks for the link once it is up
        await expect
          .poll(() =>
            server.requests.some((r) => r.path.endsWith("/mapping.json")),
          )
          .toBe(true);
        if (parts.molecules) {
          await expect(page.locator(SELECTORS.smLayerButton)).toBeVisible();
          await expect
            .poll(() => served(/\/molecules\/manifest\.json\.gz$/))
            .toBe(true);
        } else {
          await expect(page.locator(SELECTORS.smLayerButton)).toHaveCount(0);
        }
      });
    } else if (parts.molecules) {
      await test.step("lists the molecules' genes", async () => {
        await expect
          .poll(() => served(/\/molecules\/manifest\.json\.gz$/))
          .toBe(true);
        await expect(
          page.locator(SELECTORS.smGenesButton).first(),
        ).toBeVisible();
      });
    }

    if (layout === "spatialdata_v3_sharded") {
      await test.step("reads shards by byte range", async () => {
        expect(
          server.requests.some(
            (r) => r.status === 206 && r.path.includes("/tables/"),
          ),
          "a ranged read of the table",
        ).toBe(true);
      });
    }
  });
}

test("@quick SpatialData store with nothing to show says so", async ({
  page,
}) => {
  await page.goto(
    `/spatialdata-viewer/from-s3?url=${encodeURIComponent(storeUrl(server, {}))}`,
  );
  await expect(
    page.getByText(
      "This SpatialData store has no cells, molecules or image to show",
    ),
  ).toBeVisible({ timeout: 120_000 });
});

test("@quick image viewer draws an image on its own", async ({ page }) => {
  await page.goto(
    `/image-viewer/from-s3?url=${encodeURIComponent(imageUrl(server, 50))}`,
  );

  await expect(page.locator(SELECTORS.imageCanvas)).toBeVisible({
    timeout: 120_000,
  });
  await expect(page.locator(SELECTORS.imageButton)).toBeVisible();
  // Framed to fill the view, so most of the canvas is image once tiles land
  await expect
    .poll(() => litFraction(page, SELECTORS.imageCanvas), {
      message: "the image should be drawn",
      timeout: 60_000,
    })
    .toBeGreaterThan(0.2);

  await test.step("hiding the image clears the canvas", async () => {
    await page.locator(SELECTORS.imageButton).click();
    await page.getByLabel("Show image").uncheck();
    // What stays lit is the controls drawn over the canvas
    await expect
      .poll(() => litFraction(page, SELECTORS.imageCanvas), { timeout: 30_000 })
      .toBeLessThan(0.05);
  });
});

test("@quick image viewer explains a URL it can't read", async ({ page }) => {
  await page.goto(
    `/image-viewer/from-s3?url=${encodeURIComponent(`${server.url}/missing.zarr`)}`,
  );
  await expect(
    page.getByText("Couldn't read an image at that URL"),
  ).toBeVisible({
    timeout: 120_000,
  });
});
