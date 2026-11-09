import {
  test,
  expect,
  litFraction,
  makeBundle,
  makeImageFolder,
  makeTableFolder,
  SELECTORS,
} from "../../fixtures/merfish-fixtures";
import { signInDev } from "../../helpers/auth";

/**
 * "Upload & Save" from a local preview of already-processed data: every part
 * is uploaded as it is, and the saved link opens with the same overlays the
 * preview had; the save then asks for one "ready" email for all its parts
 * (/api/uploads/ready, D6 of docs/plans/PREVIEW-SAVE.md). Tagged @upload
 * (needs Postgres + S3).
 */

const UPLOAD_ENABLED = process.env.E2E_UPLOAD === "1";

type Page = import("@playwright/test").Page;
type Request = import("@playwright/test").Request;

/** Every "ready" email the page asks for, from now on. */
function watchReady(page: Page): Request[] {
  const asked: Request[] = [];

  page.on("request", (r) => {
    if (r.method() === "POST" && r.url().endsWith("/api/uploads/ready")) {
      asked.push(r);
    }
  });

  return asked;
}

/** The one "ready" email a save asked for: its datasets, and the answer. */
async function oneReady(asked: Request[]) {
  await expect.poll(() => asked.length).toBe(1);
  const response = await asked[0].response();

  expect(response?.status()).toBe(200);

  return {
    datasetIds: asked[0].postDataJSON().datasetIds as string[],
    // No SES in a test run, so nothing really goes out
    answer: await response!.json(),
  };
}

/** Drop a folder, wait for its preview, save it, and open the saved link. */
async function previewAndSave(
  page: Page,
  folder: string,
  name: string,
  previewUrl: RegExp,
): Promise<void> {
  await page.goto("/", { waitUntil: "domcontentloaded" });
  await page.locator("#sc-file-input-folder").setInputFiles(folder);
  await page.waitForURL(previewUrl, { timeout: 120_000 });
  await page.locator(SELECTORS.uploadSaveButton).first().click();
  await page.getByLabel(/Project name|Dataset name/).fill(name);
  await page.getByRole("button", { name: "Upload", exact: true }).click();
  await expect(page.getByTestId("bundle-saved")).toBeVisible({
    timeout: 180_000,
  });
  await page.getByText("Open the saved dataset").click();
}

const mine = async (page: Page, type: string, name: string) =>
  (
    (await (await page.request.get(`/api/ingest/mine?type=${type}`)).json())
      .datasets as { id: string; title: string; viewerUrl: string | null }[]
  ).filter((d) => d.title.startsWith(name));

test.describe("@upload save a local preview as it is", () => {
  test.skip(
    !UPLOAD_ENABLED,
    "Set E2E_UPLOAD=1 (and start docker infra) to run upload tests.",
  );

  test("merfisheyes cells + molecules + image stay overlaid", async ({
    page,
  }) => {
    const name = `e2e-bundle-${Date.now()}`;
    const ready = watchReady(page);

    await signInDev(page, `e2e-bundle-${Date.now()}@example.com`);
    await previewAndSave(
      page,
      makeBundle({ cells: "chunked", molecules: true, images: ["dapi.zarr"] }),
      name,
      /\/viewer\/from-local/,
    );

    await test.step("the saved cells open with both overlays", async () => {
      await page.waitForURL(/\/viewer\/ds_/, { timeout: 60_000 });
      await expect(page.locator(SELECTORS.scCanvas).first()).toBeVisible();
      await expect(page.locator(SELECTORS.imageButton)).toBeVisible();
      await expect(page.locator(SELECTORS.smLayerButton)).toBeVisible();
      await expect
        .poll(() => litFraction(page, SELECTORS.scCanvas), { timeout: 60_000 })
        .toBeGreaterThan(0.08);
    });

    await test.step("three datasets, linked, in one project", async () => {
      const [cells] = await mine(page, "single_cell", name);
      const [molecules] = await mine(page, "single_molecule", name);
      const [image] = await mine(page, "image", name);

      expect(cells.title).toBe(`${name} (cells)`);
      expect(molecules.title).toBe(`${name} (molecules)`);
      expect(image.title).toBe(`${name} (image)`);
      expect(
        await (
          await page.request.get(`/api/datasets/${cells.id}/mapping`)
        ).json(),
      ).toMatchObject({
        links: { __all__: molecules.id },
        source: "app",
        image: image.id,
        imageSource: "app",
      });
    });

    await test.step("one ready email for the three, the cells first", async () => {
      const [cells] = await mine(page, "single_cell", name);
      const [molecules] = await mine(page, "single_molecule", name);
      const [image] = await mine(page, "image", name);
      const { datasetIds, answer } = await oneReady(ready);

      expect(datasetIds).toEqual([cells.id, molecules.id, image.id]);
      // All three are the user's and finished: the email was for all of them
      expect(answer).toMatchObject({ success: true, counted: 3 });
    });
  });

  test("a SpatialData store + molecules stay overlaid", async ({ page }) => {
    const name = `e2e-store-${Date.now()}`;

    await signInDev(page, `e2e-store-${Date.now()}@example.com`);
    await previewAndSave(
      page,
      makeBundle({ cells: "spatialdata", molecules: true }),
      name,
      /\/viewer\/from-local/,
    );

    // The store is saved as it is and read by the SpatialData viewer
    await page.waitForURL(/\/spatialdata-viewer\/ds_/, { timeout: 60_000 });
    await expect(page.locator(SELECTORS.scCanvas).first()).toBeVisible({
      timeout: 120_000,
    });
    await expect(page.locator(SELECTORS.imageButton)).toBeVisible();
    await expect(page.locator(SELECTORS.smLayerButton)).toBeVisible();

    const [store] = await mine(page, "spatialdata", name);

    expect(store.viewerUrl).toBe(`/spatialdata-viewer/${store.id}`);
  });

  test("a table on its own is saved as a store", async ({ page }) => {
    const name = `e2e-table-${Date.now()}`;

    await signInDev(page, `e2e-table-${Date.now()}@example.com`);
    await previewAndSave(
      page,
      makeTableFolder("table"),
      name,
      /\/viewer\/from-local/,
    );

    // Saved where a store keeps its table, with the root metadata it lacked
    await page.waitForURL(/\/spatialdata-viewer\/ds_/, { timeout: 60_000 });
    await expect(page.locator(SELECTORS.scCanvas).first()).toBeVisible({
      timeout: 120_000,
    });
    await page.waitForFunction(
      () => (window as any).__merfish?.stats?.pointCount === 30,
      undefined,
      { timeout: 120_000 },
    );
  });

  test("an image on its own", async ({ page }) => {
    const name = `e2e-image-${Date.now()}`;
    const ready = watchReady(page);

    await signInDev(page, `e2e-local-image-${Date.now()}@example.com`);
    await previewAndSave(
      page,
      makeImageFolder(),
      name,
      /\/image-viewer\/from-local/,
    );
    await page.waitForURL(/\/image-viewer\/ds_/, { timeout: 60_000 });
    await expect
      .poll(() => litFraction(page, SELECTORS.imageCanvas), { timeout: 60_000 })
      .toBeGreaterThan(0.2);
    const [image] = await mine(page, "image", name);

    expect(image.title).toBe(name);
    const { datasetIds, answer } = await oneReady(ready);

    expect(datasetIds).toEqual([image.id]);
    expect(answer).toMatchObject({ success: true, counted: 1 });
  });
});
