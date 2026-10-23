import path from "node:path";

import {
  test,
  expect,
  REPO_ROOT,
  SELECTORS,
} from "../../fixtures/merfish-fixtures";
import { signInDev } from "../../helpers/auth";

/**
 * A pre-chunked molecule folder (manifest.json.gz + genes/*.bin.gz) opened in
 * single molecule mode and saved with Upload & Save — the path of
 * kresnajenie/merfisheyes#166. Tagged @upload.
 */

const UPLOAD_ENABLED = process.env.E2E_UPLOAD === "1";
const CHUNKED = path.join(REPO_ROOT, "tests", "data", "tiny", "single-molecule", "chunked");

test.describe("@upload pre-chunked molecules: Upload & Save", () => {
  test.skip(!UPLOAD_ENABLED, "Set E2E_UPLOAD=1 (and start docker infra) to run upload tests.");

  test("the folder is saved and opens from the app's storage", async ({ page }) => {
    const name = `e2e-sm-chunked-${Date.now()}`;

    await signInDev(page, `e2e-sm-chunked-${Date.now()}@example.com`);
    await page.goto("/", { waitUntil: "domcontentloaded" });
    await page.getByRole("button", { name: /single molecule/i }).first().click();
    await page.locator("#sm-file-input-chunked").setInputFiles(CHUNKED);
    await page.waitForURL("**/sm-viewer/from-local", { timeout: 60_000 });
    await expect(page.locator(SELECTORS.smCanvas).first()).toBeVisible({ timeout: 60_000 });

    await page.locator(SELECTORS.uploadSaveButton).first().click();
    await page.getByPlaceholder("Enter dataset name").fill(name);
    await page.locator(SELECTORS.processUploadButton).click();

    await expect(page.locator(SELECTORS.uploadSuccessHeading).first()).toBeVisible({
      timeout: 120_000,
    });
    const href = await page.locator('a[href*="/sm-viewer/"]').first().getAttribute("href");
    const id = href?.match(/\/sm-viewer\/([^/?#]+)/)?.[1];

    expect(id).toBeTruthy();
    await page.goto(`/sm-viewer/${id}`);
    await expect(page.locator(SELECTORS.smCanvas).first()).toBeVisible({ timeout: 60_000 });
  });

  test("a failed save keeps its message on screen until it is closed", async ({ page }) => {
    await signInDev(page, `e2e-sm-chunked-fail-${Date.now()}@example.com`);
    await page.goto("/", { waitUntil: "domcontentloaded" });
    await page.getByRole("button", { name: /single molecule/i }).first().click();
    await page.locator("#sm-file-input-chunked").setInputFiles(CHUNKED);
    await page.waitForURL("**/sm-viewer/from-local", { timeout: 60_000 });

    // The server turns the upload down, as it might for a real reason
    await page.route("**/api/single-molecule/initiate", (route) =>
      route.fulfill({ status: 500, json: { error: "Storage is unavailable" } }),
    );
    const title = `e2e-sm-fail-${Date.now()}`;
    const reported = page.waitForRequest("**/api/uploads/failed");

    await page.locator(SELECTORS.uploadSaveButton).first().click();
    await page.getByPlaceholder("Enter dataset name").fill(title);
    await page.locator(SELECTORS.processUploadButton).click();

    // The reason goes by email too; no dataset was made to mark
    expect((await reported).postDataJSON()).toEqual({
      title,
      error: "Storage is unavailable",
    });

    const message = page.getByText("Failed to process dataset: Storage is unavailable");

    await expect(message).toBeVisible();
    // Toasts close after 3 s by default; an error waits to be read
    await page.waitForTimeout(10_000);
    await expect(message).toBeVisible();
    // One click closes it, and the dialog behind it stays open to try again
    await page
      .locator(".Toastify__toast", { hasText: "Storage is unavailable" })
      .locator(".Toastify__close-button")
      .click();
    await expect(message).toBeHidden();
    await expect(page.getByPlaceholder("Enter dataset name")).toBeVisible();
  });

  test("a save that fails after its dataset was made marks it FAILED", async ({ page }) => {
    await signInDev(page, `e2e-sm-chunked-late-${Date.now()}@example.com`);
    await page.goto("/", { waitUntil: "domcontentloaded" });
    await page.getByRole("button", { name: /single molecule/i }).first().click();
    await page.locator("#sm-file-input-chunked").setInputFiles(CHUNKED);
    await page.waitForURL("**/sm-viewer/from-local", { timeout: 60_000 });

    // The files go up; finishing the upload is what fails
    await page.route("**/api/single-molecule/*/complete", (route) =>
      route.fulfill({ status: 500, json: { error: "Could not finish the upload" } }),
    );
    const reported = page.waitForResponse("**/api/uploads/failed");

    await page.locator(SELECTORS.uploadSaveButton).first().click();
    await page.getByPlaceholder("Enter dataset name").fill(`e2e-sm-late-${Date.now()}`);
    await page.locator(SELECTORS.processUploadButton).click();

    const response = await reported;
    const { datasetId } = response.request().postDataJSON();

    expect(datasetId).toMatch(/^sm_/);
    expect(await response.json()).toMatchObject({ marked: true });
    const dataset = await (await page.request.get(`/api/single-molecule/${datasetId}`)).json();

    expect(dataset.status ?? dataset.dataset?.status).toBe("FAILED");
  });
});
