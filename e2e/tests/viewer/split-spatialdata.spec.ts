import {
  test,
  expect,
  startSpatialDataServer,
  storeUrl,
  imageUrl,
  SELECTORS,
} from "../../fixtures/merfish-fixtures";
import { type FixtureServer } from "../../helpers/fixture-server";

/**
 * A split view's right panel opens what its viewer would on a page of its
 * own: a SpatialData store (cells, or molecules alone) and an image, by a
 * shared split link or a link pasted into the picker; and a SpatialData page
 * restores its split from the URL. Tagged @quick: the stores are served by a
 * local stand-in for a bucket.
 */

let server: FixtureServer;

test.beforeAll(async () => {
  server = await startSpatialDataServer();
});
test.afterAll(() => server.close());

const fromS3 = (viewer: string, url: string) =>
  `/${viewer}/from-s3?url=${encodeURIComponent(url)}`;

test("@quick a split link opens a molecules store beside a cells store", async ({
  page,
}) => {
  const cells = storeUrl(server, { cells: true });
  const molecules = storeUrl(server, { points: true });

  await page.goto(
    `${fromS3("spatialdata-viewer", cells)}&splitS3Url=${encodeURIComponent(
      encodeURIComponent(molecules),
    )}&splitType=spatialdata`,
  );
  // Cells on the left, the store's molecules alone on the right
  await expect(page.locator(SELECTORS.scCanvas).first()).toBeVisible({
    timeout: 120_000,
  });
  await expect(page.locator(SELECTORS.smCanvas).first()).toBeVisible({
    timeout: 120_000,
  });
  // Drawn by the tiled points layer, which says what it drew
  await expect(page.getByText(/[\d.]+ ?k? drawn/).first()).toBeVisible({
    timeout: 60_000,
  });
  // The split stays in the link
  expect(new URL(page.url()).searchParams.get("splitType")).toBe("spatialdata");
});

test("@quick the picker opens the same store, a pasted store and a pasted image", async ({
  page,
}) => {
  await page.goto(fromS3("spatialdata-viewer", storeUrl(server, { cells: true })));
  await expect(page.locator(SELECTORS.scCanvas).first()).toBeVisible({
    timeout: 120_000,
  });

  // "Same Dataset": the store again, as a store
  await page.getByRole("button", { name: "Split screen" }).click();
  await page.getByText("Same Dataset").click();
  await expect(page.locator(SELECTORS.scCanvas)).toHaveCount(2, {
    timeout: 120_000,
  });
  expect(new URL(page.url()).searchParams.get("splitType")).toBe("spatialdata");
  await page.getByRole("button", { name: "Close split view" }).click();

  // A SpatialData store's link: its cells in the panel
  await page.getByRole("button", { name: "Split screen" }).click();
  await page.getByText("Paste Link").click();
  await page
    .getByPlaceholder("Paste dataset URL or ID...")
    .fill(`${server.url}${fromS3("spatialdata-viewer", storeUrl(server, { cells: true, image: true }))}`);
  await page.getByRole("button", { name: "Load" }).click();
  await expect(page.locator(SELECTORS.scCanvas)).toHaveCount(2, {
    timeout: 120_000,
  });
  expect(new URL(page.url()).searchParams.get("splitType")).toBe("spatialdata");

  // An image's link: the image in the panel, inside it (not over the page)
  await page.getByRole("button", { name: "Close split view" }).click();
  await page.getByRole("button", { name: "Split screen" }).click();
  await page.getByText("Paste Link").click();
  await page
    .getByPlaceholder("Paste dataset URL or ID...")
    .fill(`${server.url}${fromS3("image-viewer", imageUrl(server, 1000))}`);
  await page.getByRole("button", { name: "Load" }).click();
  const image = page.locator(SELECTORS.imageCanvas).first();

  await expect(image).toBeVisible({ timeout: 120_000 });
  await expect
    .poll(() => new URL(page.url()).searchParams.get("splitType"))
    .toBe("image");
  const panel = await image.boundingBox();
  const viewport = page.viewportSize()!;

  expect(panel!.x).toBeGreaterThan(viewport.width * 0.3);
  // The cells are still on the left
  await expect(page.locator(SELECTORS.scCanvas).first()).toBeVisible();
});
