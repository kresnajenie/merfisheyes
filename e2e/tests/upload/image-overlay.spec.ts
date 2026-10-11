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
  SELECTORS,
} from "../../fixtures/merfish-fixtures";
import { signInDev } from "../../helpers/auth";
import { type FixtureServer } from "../../helpers/fixture-server";

/**
 * Image overlay on a saved single cell dataset: its owner attaches an image
 * from the overlay button, the viewer draws it under the cells, and only the
 * owner can change it. Tagged @upload (needs Postgres + S3, see
 * upload-save.spec.ts); the image is served by a local stand-in for a public
 * bucket.
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

  test("owner attaches, sees and removes an image", async ({
    page,
    browser,
    baseURL,
  }, testInfo) => {
    // Its cells keep their raw coordinates (about -1000..1000), which is what
    // an image is placed in.
    const ds = getDataset("h5ad-tiny");
    const image = imageUrl(server, 1000);
    const tileServed = () =>
      server.requests.some(
        (r) => r.status === 200 && /\/image\/1000\.zarr\/\d\/c\//.test(r.path),
      );

    annotate(testInfo, ds.id, ds.format);

    let datasetId = "";

    await test.step("sign in, upload and open the saved dataset", async () => {
      await signInDev(page);
      await loadLocalDataset(page, ds);
      datasetId = (await uploadAndSave(page, ds, { name: "e2e-image-overlay" }))
        .datasetId;
      await openSavedDataset(page, ds, datasetId);
      await expect(page.locator(SELECTORS.imageButton)).toHaveCount(0);
    });

    const overlay = `/api/datasets/${datasetId}/overlay`;

    await test.step("a URL with no image behind it is refused", async () => {
      await page.locator(SELECTORS.overlayButton).click();
      await page.getByLabel("Image URL").fill(`${server.url}/missing.zarr`);
      await page.getByRole("button", { name: "Attach" }).click();
      await expect(
        page.getByText("Couldn't read an image at that URL"),
      ).toBeVisible();
      expect((await (await page.request.get(overlay)).json()).image).toBeNull();
    });

    await test.step("attach the image", async () => {
      await page.getByLabel("Image URL").fill(image);
      await page.getByRole("button", { name: "Attach" }).click();
      await expect(page.getByText("Image overlay attached")).toBeVisible();
      await expect(page.getByText("Current image")).toBeVisible();
      expect((await (await page.request.get(overlay)).json()).image).toBe(
        image,
      );
    });

    await test.step("the viewer draws it under the cells", async () => {
      await openSavedDataset(page, ds, datasetId);
      await expect(page.locator(SELECTORS.imageButton)).toBeVisible();
      await expect
        .poll(tileServed, { message: "an image tile should be fetched" })
        .toBe(true);
    });

    await test.step("only the owner can change it", async () => {
      const change = { data: { image: null } };
      const anonymous = await browser.newContext({ baseURL });

      expect((await anonymous.request.patch(overlay, change)).status()).toBe(
        401,
      );
      await anonymous.close();

      const other = await browser.newContext({ baseURL });

      await signInDev(await other.newPage(), "e2e-other@example.com");
      expect((await other.request.patch(overlay, change)).status()).toBe(404);
      await other.close();

      expect(
        (
          await page.request.patch(overlay, {
            data: { image: "s3://bucket/image.zarr" },
          })
        ).status(),
        "only http(s) URLs are accepted",
      ).toBe(400);
      expect((await (await page.request.get(overlay)).json()).image).toBe(
        image,
      );
    });

    await test.step("remove the image", async () => {
      page.once("dialog", (dialog) => dialog.accept());
      await page.locator(SELECTORS.overlayButton).click();
      await page.getByRole("button", { name: "Remove" }).click();
      await expect(page.getByText("Image overlay removed")).toBeVisible();
      expect((await (await page.request.get(overlay)).json()).image).toBeNull();

      await openSavedDataset(page, ds, datasetId);
      await expect(page.locator(SELECTORS.imageButton)).toHaveCount(0);
    });
  });
});
