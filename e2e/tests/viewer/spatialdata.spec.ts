import {
  test,
  expect,
  readHooks,
  litFraction,
  startSpatialDataServer,
  storeUrl,
  storeFolderUrl,
  imageUrl,
  visiumUrl,
  visiumHdUrl,
  libraryStoreUrl,
  selectGeneAndMeasure,
  SELECTORS,
  type DatasetEntry,
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

test("@quick Visium store: the full image under the spots, not the downscaled ones", async ({
  page,
}) => {
  // No `global` system: the spots' system is the target, and the full image
  // is identity into it, so it spans the spots and no more: a square in the
  // middle of the view framed on them. Multiplied by every transform of the
  // spots (as before), it was drawn 64x too large and ran off the view.
  await page.goto(
    `/spatialdata-viewer/from-s3?url=${encodeURIComponent(visiumUrl(server))}`,
  );
  await expect(page.locator(SCENES.cells).first()).toBeVisible({
    timeout: 120_000,
  });
  await page.waitForFunction(
    () => (window as any).__merfish?.renderComplete != null,
    undefined,
    { timeout: 120_000 },
  );
  expect((await readHooks(page))!.stats!.pointCount, "spots").toBe(36);
  await expect(page.locator(SELECTORS.imageButton)).toBeVisible();
  await expect
    .poll(() => served(/\/images\/tiny_full_image\/s\d\/c\//), {
      message: "the full image should be read",
    })
    .toBe(true);
  await expect
    .poll(() => litFraction(page, SCENES.cells), {
      message: "the image should be drawn",
      timeout: 60_000,
    })
    .toBeGreaterThan(0.08);
  // A band along the right edge, clear of the buttons: dark beside the
  // spots, lit when the image runs off past them
  expect(
    await litFraction(page, SCENES.cells, {
      x0: 0.92,
      y0: 0.3,
      x1: 0.99,
      y1: 0.7,
    }),
    "the image should cover the spots, not run past them",
  ).toBeLessThan(0.05);
  expect(
    served(/\/images\/tiny_(hires|lowres)_image\//),
    "the downscaled images should not be read",
  ).toBe(false);
});

test("@quick Visium HD store: the bins drawn as squares of their real size", async ({
  page,
}) => {
  // Square bins of side 16 on a 16 grid: at their real size they tile the
  // section (1 in 7 of the fixture's bins is missing: about 86% of it lit
  // with the image hidden); as dots, almost nothing is lit
  await page.goto(
    `/spatialdata-viewer/from-s3?url=${encodeURIComponent(visiumHdUrl(server))}`,
  );
  await expect(page.locator(SCENES.cells).first()).toBeVisible({
    timeout: 120_000,
  });
  await page.waitForFunction(
    () => (window as any).__merfish?.renderComplete != null,
    undefined,
    { timeout: 120_000 },
  );
  expect((await readHooks(page))!.stats!.pointCount, "bins").toBe(219);
  const middle = { x0: 0.45, y0: 0.45, x1: 0.55, y1: 0.55 };
  const binsLit = () => litFraction(page, SCENES.cells, middle);

  // The bins alone: no column to colour by, so they are translucent
  // neutral (NEUTRAL_POINT_ALPHA), and the image is hidden
  await page.locator(SELECTORS.imageButton).click();
  await page.getByLabel("Show image").uncheck();
  // Close the image panel: it sits over the rail's buttons
  await page.locator(SELECTORS.imageButton).click();
  await expect
    .poll(binsLit, {
      message: "the bins should tile the section",
      timeout: 60_000,
    })
    .toBeGreaterThan(0.7);

  await test.step("coloured by a gene, every bin keeps its size", async () => {
    // Low expression shrank a dot; a bin stays a bin (colour carries it)
    await selectGeneAndMeasure(
      page,
      { dataType: "single_cell" } as DatasetEntry,
      "gene0",
    );
    await expect
      .poll(binsLit, {
        message: "the bins should still tile the section",
        timeout: 30_000,
      })
      .toBeGreaterThan(0.7);
  });

  await test.step("the SC button hides them (an image is under them)", async () => {
    await page.getByRole("button", { name: "SC", exact: true }).click();
    await expect
      .poll(binsLit, { message: "the bins should be hidden", timeout: 30_000 })
      .toBeLessThan(0.05);
  });
});

test("@quick CODEX store: a dense X, a 6-channel image shown one channel at a time", async ({
  page,
}) => {
  await page.goto(
    `/spatialdata-viewer/from-s3?url=${encodeURIComponent(libraryStoreUrl(server, "codex"))}`,
  );
  await expect(page.locator(SCENES.cells).first()).toBeVisible({
    timeout: 120_000,
  });
  await page.waitForFunction(
    () => (window as any).__merfish?.renderComplete != null,
    undefined,
    { timeout: 120_000 },
  );
  expect((await readHooks(page))!.stats!.pointCount, "cells").toBe(36);
  await expect(page.locator(SELECTORS.imageButton)).toBeVisible();
  await expect
    .poll(() => served(/\/images\/reg001_image\/s\d\/c\/0\//), {
      message: "the first channel should be read",
    })
    .toBe(true);
  await expect
    .poll(() => litFraction(page, SCENES.cells), {
      message: "the image should be drawn",
      timeout: 60_000,
    })
    .toBeGreaterThan(0.08);
  // Only the shown channel is read; the panel lists every one
  expect(served(/\/images\/reg001_image\/s\d\/c\/[1-5]\//)).toBe(false);
  await page.locator(SELECTORS.imageButton).click();
  await expect(page.locator('input[aria-label$=" colour"]')).toHaveCount(6);
  await page.getByLabel("Show CD3", { exact: true }).check();
  await expect
    .poll(() => served(/\/images\/reg001_image\/s\d\/c\/1\//), {
      message: "a channel switched on should be read",
    })
    .toBe(true);
  // A marker's intensities colour the cells
  await page.locator(SELECTORS.imageButton).click();
  await selectGeneAndMeasure(
    page,
    { dataType: "single_cell" } as DatasetEntry,
    "CD3",
  );
});

test("@quick MALDI-MSI store: ions as features, the optical image placed by a rotating affine", async ({
  page,
}) => {
  await page.goto(
    `/spatialdata-viewer/from-s3?url=${encodeURIComponent(libraryStoreUrl(server, "maldi"))}`,
  );
  await expect(page.locator(SCENES.cells).first()).toBeVisible({
    timeout: 120_000,
  });
  await page.waitForFunction(
    () => (window as any).__merfish?.renderComplete != null,
    undefined,
    { timeout: 120_000 },
  );
  expect((await readHooks(page))!.stats!.pointCount, "spots").toBe(64);
  // obsm/spatial is the ion image's pixel index: the spots are placed from
  // their own Point element, read whole
  expect(
    server.requests.some(
      (r) =>
        r.status < 300 &&
        !r.range &&
        r.path.endsWith("/shapes/maldi_spots/shapes.parquet"),
    ),
  ).toBe(true);
  await expect(page.locator(SELECTORS.imageButton)).toBeVisible();
  await expect
    .poll(() => served(/\/images\/optical_image\/s\d\/c\//), {
      message: "the optical image should be read",
    })
    .toBe(true);
  await expect
    .poll(() => litFraction(page, SCENES.cells), {
      message: "the rotated image should be drawn",
      timeout: 60_000,
    })
    .toBeGreaterThan(0.08);
  await selectGeneAndMeasure(
    page,
    { dataType: "single_cell" } as DatasetEntry,
    "C6H12O6+Na",
  );
});

test("@quick CosMx store: cells placed from their shapes when the table has no obsm/spatial", async ({
  page,
}) => {
  await page.goto(
    `/spatialdata-viewer/from-s3?url=${encodeURIComponent(libraryStoreUrl(server, "cosmx"))}`,
  );
  await expect(page.locator(SCENES.cells).first()).toBeVisible({
    timeout: 120_000,
  });
  await page.waitForFunction(
    () => (window as any).__merfish?.renderComplete != null,
    undefined,
    { timeout: 120_000 },
  );
  expect((await readHooks(page))!.stats!.pointCount, "cells").toBe(25);
  // The Point element is read whole; of the polygons (265 MB in the real
  // store) only the footer, by byte range, for the cells' footprint
  const parquet = (name: string) =>
    server.requests.filter(
      (r) =>
        r.status < 300 && r.path.endsWith(`/shapes/${name}/shapes.parquet`),
    );

  expect(parquet("cells_polygons_for_vitessce").some((r) => !r.range)).toBe(
    true,
  );
  expect(parquet("cells_polygons").every((r) => r.range)).toBe(true);
  await expect(page.locator(SELECTORS.imageButton)).toHaveCount(0);
  await selectGeneAndMeasure(
    page,
    { dataType: "single_cell" } as DatasetEntry,
    "CD45",
  );
});

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

/** Paste a URL into the homepage's "Load from S3" box. */
async function loadFromS3(
  page: import("@playwright/test").Page,
  url: string,
): Promise<void> {
  await page.goto("/", { waitUntil: "domcontentloaded" });
  await page.getByText("Load from S3").first().click();
  await page.getByLabel("Dataset Folder URL").fill(url);
  await page.getByRole("button", { name: "Load Dataset" }).click();
}

test("@quick Load from S3: a folder holding a store opens overlaid", async ({
  page,
}) => {
  const parts = { cells: true, molecules: true, image: true };

  // Only the folder is given; its mapping.json names the store inside it
  await loadFromS3(page, storeFolderUrl(server, parts));
  await page.waitForURL(/\/spatialdata-viewer\/from-s3\?url=/, {
    timeout: 60_000,
  });
  await expect(page.locator(SELECTORS.scCanvas).first()).toBeVisible({
    timeout: 120_000,
  });
  await page.waitForFunction(
    () => (window as any).__merfish?.renderComplete != null,
    undefined,
    { timeout: 120_000 },
  );
  expect((await readHooks(page))!.stats!.pointCount).toBe(CELLS);
  await expect(page.locator(SELECTORS.imageButton)).toBeVisible();
  await expect(page.locator(SELECTORS.smLayerButton)).toBeVisible();

  await test.step("the link reopens all of it", async () => {
    await page.reload();
    await expect(page.locator(SELECTORS.scCanvas).first()).toBeVisible({
      timeout: 120_000,
    });
    await expect(page.locator(SELECTORS.imageButton)).toBeVisible();
    await expect(page.locator(SELECTORS.smLayerButton)).toBeVisible();
  });
});

test("@quick Load from S3: a store and an image open in their viewers", async ({
  page,
}) => {
  await loadFromS3(page, storeUrl(server, { cells: true }));
  await page.waitForURL(/\/spatialdata-viewer\/from-s3\?url=/, {
    timeout: 60_000,
  });
  await expect(page.locator(SELECTORS.scCanvas).first()).toBeVisible({
    timeout: 120_000,
  });

  await loadFromS3(page, imageUrl(server, 50));
  await page.waitForURL(/\/image-viewer\/from-s3\?url=/, { timeout: 60_000 });
  await expect(page.locator(SELECTORS.imageCanvas)).toBeVisible({
    timeout: 120_000,
  });
});

test("@quick Load from S3: says so when there is nothing to open", async ({
  page,
}) => {
  await loadFromS3(page, `${server.url}/nothing-here`);
  await expect(page.getByText(/Nothing to open at that URL/)).toBeVisible({
    timeout: 60_000,
  });
});

test("@quick image settings are kept in the link", async ({ page }) => {
  await page.goto(
    `/image-viewer/from-s3?url=${encodeURIComponent(imageUrl(server, 50))}`,
  );
  await expect(page.locator(SELECTORS.imageCanvas)).toBeVisible({
    timeout: 120_000,
  });
  await page.locator(SELECTORS.imageButton).click();

  const colour = page.getByLabel(/colour$/).first();
  const maximum = page.getByLabel(/maximum$/).first();

  // An untouched image leaves the link alone
  expect(new URL(page.url()).searchParams.has("img")).toBe(false);

  const img = () => new URL(page.url()).searchParams.get("img");

  await colour.fill("#ff8800");
  // The colour alone writes the settings (with the auto-contrast window):
  // wait for that, then for the link to change again with the new maximum,
  // or the reload can race the debounced write
  await expect.poll(img).not.toBeNull();
  const afterColour = img();

  await maximum.fill("120");
  await expect.poll(img).not.toBe(afterColour);

  await page.reload();
  await page.locator(SELECTORS.imageButton).click();
  await expect(colour).toHaveValue("#ff8800");
  await expect(maximum).toHaveValue("120");
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
