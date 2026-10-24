import {
  test,
  expect,
  getDataset,
  annotate,
  loadLocalDataset,
  uploadAndSave,
  openSavedDataset,
  startSpatialDataServer,
  imageUrl,
  storeUrl,
  SELECTORS,
} from "../../fixtures/merfish-fixtures";
import { signInDev } from "../../helpers/auth";
import { type FixtureServer } from "../../helpers/fixture-server";

/**
 * An image is a dataset: opened by URL it can be added to an account, and the
 * owner of a single cell dataset then picks it from their images to draw it
 * under the cells. Tagged @upload (needs Postgres + S3, see
 * upload-save.spec.ts); the images are served by a local stand-in for a
 * public bucket.
 */

const UPLOAD_ENABLED = process.env.E2E_UPLOAD === "1";

test.describe("@upload image overlay", () => {
  test.skip(
    !UPLOAD_ENABLED,
    "Set E2E_UPLOAD=1 (and start docker infra) to run upload tests.",
  );

  let server: FixtureServer;

  test.beforeAll(async () => {
    server = await startSpatialDataServer();
  });
  test.afterAll(() => server?.close());

  test("a SpatialData store opened by URL can be added to the account", async ({
    page,
  }) => {
    const store = storeUrl(server, { cells: true, image: true });

    await signInDev(page, `e2e-store-url-${Date.now()}@example.com`);
    await page.goto(
      `/spatialdata-viewer/from-s3?url=${encodeURIComponent(store)}`,
    );
    await expect(page.locator(SELECTORS.scCanvas).first()).toBeVisible({
      timeout: 120_000,
    });
    await page.getByRole("button", { name: "Add to my account" }).click();
    await expect(page.getByText("Dataset added to your account")).toBeVisible();

    const mine = (
      await (await page.request.get("/api/ingest/mine?type=spatialdata")).json()
    ).datasets as { title: string; viewerUrl: string }[];

    // Listed with the user's datasets, and opened by the same link
    expect(mine.map((d) => d.title)).toEqual(["store"]);
    expect(mine[0].viewerUrl).toBe(
      `/spatialdata-viewer/from-s3?url=${encodeURIComponent(store)}`,
    );
  });

  test("owner adds an image to their account and draws it under the cells", async ({
    page,
    browser,
    baseURL,
  }, testInfo) => {
    // Its cells keep their raw coordinates (about -1000..1000), which is what
    // an image is placed in.
    const ds = getDataset("h5ad-tiny");
    const image = imageUrl(server, 1000);
    const second = imageUrl(server, 500);
    const tileServed = () =>
      server.requests.some(
        (r) => r.status === 200 && /\/image\/1000\.zarr\/\d\/c\//.test(r.path),
      );
    const myImages = async () =>
      (await (await page.request.get("/api/ingest/mine?type=image")).json())
        .datasets as { id: string; title: string }[];

    annotate(testInfo, ds.id, ds.format);

    let datasetId = "";
    let imageId = "";

    await test.step("sign in, upload and open the saved dataset", async () => {
      // A fresh account, so "my images" starts empty
      await signInDev(page, `e2e-image-${Date.now()}@example.com`);
      await loadLocalDataset(page, ds);
      datasetId = (await uploadAndSave(page, ds, { name: "e2e-image-overlay" }))
        .datasetId;
      expect(await myImages()).toHaveLength(0);
    });

    const overlay = `/api/datasets/${datasetId}/overlay`;
    const stored = async () => (await page.request.get(overlay)).json();

    await test.step("an image opened by URL can be added to the account", async () => {
      await page.goto(`/image-viewer/from-s3?url=${encodeURIComponent(image)}`);
      await expect(page.locator(SELECTORS.imageCanvas)).toBeVisible();
      await page.getByRole("button", { name: "Add to my account" }).click();
      await expect(
        page.getByText("Dataset added to your account"),
      ).toBeVisible();
      await expect(
        page.getByRole("button", { name: "Add to my account" }),
      ).toHaveCount(0);

      const mine = await myImages();

      expect(mine.map((d) => d.title)).toEqual(["image/1000"]);
      imageId = mine[0].id;
    });

    await test.step("the overlay dialog lists it; picking it attaches it", async () => {
      await openSavedDataset(page, ds, datasetId);
      await expect(page.locator(SELECTORS.imageButton)).toHaveCount(0);
      await page.locator(SELECTORS.overlayButton).click();
      await page
        .getByTestId("image-overlay-list")
        .getByRole("button", { name: /image\/1000/ })
        .click();
      await expect(page.getByText("Image overlay attached")).toBeVisible();
      await expect(page.getByText("Current image")).toBeVisible();
      expect(await stored()).toMatchObject({ image, imageDatasetId: imageId });
    });

    await test.step("the viewer draws it under the cells", async () => {
      await openSavedDataset(page, ds, datasetId);
      await expect(page.locator(SELECTORS.imageButton)).toBeVisible();
      await expect
        .poll(tileServed, { message: "an image tile should be fetched" })
        .toBe(true);
    });

    await test.step("each channel has a lower and an upper bound", async () => {
      await page.locator(SELECTORS.imageButton).click();
      const min = page.getByLabel("DAPI minimum", { exact: true });
      const max = page.getByLabel("DAPI maximum", { exact: true });

      await expect(min).toBeVisible();
      await max.fill("40000");
      await min.fill("1000");
      const minValue = page.getByLabel("DAPI minimum value", { exact: true });
      const maxValue = page.getByLabel("DAPI maximum value", { exact: true });

      await expect(minValue).toHaveValue("1000");
      await expect(maxValue).toHaveValue("40000");
      // The bounds can't cross
      await min.fill("65535");
      await expect(minValue).toHaveValue("39999");
      await expect(maxValue).toHaveValue("40000");
    });

    await test.step("an image added by URL joins the account and replaces it", async () => {
      await page.locator(SELECTORS.overlayButton).click();
      await page.getByLabel("Image URL").fill(`${server.url}/missing.zarr`);
      await page.getByRole("button", { name: "Add", exact: true }).click();
      await expect(
        page.getByText("Couldn't read an image at that URL"),
      ).toBeVisible();
      expect((await stored()).image).toBe(image);

      await page.getByLabel("Image URL").fill(second);
      await page.getByRole("button", { name: "Add", exact: true }).click();
      await expect(
        page.getByText("Image added to your account and attached"),
      ).toBeVisible();
      expect((await stored()).image).toBe(second);
      expect((await myImages()).map((d) => d.title).sort()).toEqual([
        "image/1000",
        "image/500",
      ]);
    });

    await test.step("only the owner can change it, and only to their own image", async () => {
      const anonymous = await browser.newContext({ baseURL });

      expect(
        (
          await anonymous.request.patch(overlay, { data: { image: null } })
        ).status(),
      ).toBe(401);
      await anonymous.close();

      const other = await browser.newContext({ baseURL });

      await signInDev(await other.newPage(), "e2e-other@example.com");
      expect(
        (
          await other.request.patch(overlay, { data: { image: null } })
        ).status(),
        "someone else's dataset",
      ).toBe(404);
      await other.close();

      expect(
        (
          await page.request.patch(overlay, {
            data: { imageDatasetId: datasetId },
          })
        ).status(),
        "a dataset that is not an image",
      ).toBe(400);
      expect(
        (
          await page.request.patch(overlay, {
            data: { image: "s3://bucket/image.zarr" },
          })
        ).status(),
        "only http(s) URLs are accepted",
      ).toBe(400);
      expect((await stored()).image).toBe(second);
    });

    await test.step("deleting the image stops the cells from drawing it", async () => {
      const secondId = (await stored()).imageDatasetId;

      expect(
        (await page.request.delete(`/api/ingest/${secondId}`)).ok(),
      ).toBeTruthy();
      expect(await stored()).toMatchObject({
        image: null,
        imageDatasetId: null,
      });
      await openSavedDataset(page, ds, datasetId);
      await expect(page.locator(SELECTORS.imageButton)).toHaveCount(0);
    });
  });
});
