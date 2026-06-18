import { type Page, expect } from "@playwright/test";

import { type DatasetEntry } from "./databank";
import { SELECTORS } from "./selectors";
import { waitForRender, readHooks, type MerfishHooks } from "./viewer";

const UPLOAD_TIMEOUT = 180_000;

export interface UploadResult {
  datasetId: string;
  /** Total time from confirm-click to success state (ms, wall clock). */
  totalMs: number;
}

/**
 * Drive the "Upload & Save" modal for a dataset already open in the local
 * viewer. Works for both single-cell and single-molecule paths. Requires the
 * app to have a reachable DB + S3 (see docker-compose.test.yml).
 */
export async function uploadAndSave(
  page: Page,
  ds: DatasetEntry,
  opts: { name?: string; email?: string } = {},
): Promise<UploadResult> {
  const name = opts.name ?? `e2e-${ds.id}-${Date.now()}`;
  const email = opts.email ?? "e2e@example.com";

  await page.locator(SELECTORS.uploadSaveButton).first().click();

  await page.getByPlaceholder("Enter dataset name").fill(name);
  await page.getByPlaceholder("Enter your email").fill(email);

  const confirm =
    ds.dataType === "single_molecule"
      ? page.locator(SELECTORS.processUploadButton)
      : page.locator(SELECTORS.processSaveButton);
  await expect(confirm).toBeEnabled();

  const t0 = Date.now();
  await confirm.click();

  // Success surfaces as a heading ("Upload Successful" / "Upload Complete!").
  await expect(page.locator(SELECTORS.uploadSuccessHeading).first()).toBeVisible({
    timeout: UPLOAD_TIMEOUT,
  });
  const totalMs = Date.now() - t0;

  const datasetId = await extractDatasetId(page, ds);
  return { datasetId, totalMs };
}

/** Pull the saved dataset id from the success link in the modal. */
async function extractDatasetId(page: Page, ds: DatasetEntry): Promise<string> {
  const seg = ds.dataType === "single_molecule" ? "sm-viewer" : "viewer";
  const link = page.locator(`a[href*="/${seg}/"]`).first();
  await expect(link).toBeVisible({ timeout: 10_000 });
  const href = (await link.getAttribute("href")) ?? "";
  const m = href.match(new RegExp(`/${seg}/([^/?#]+)`));
  if (!m) throw new Error(`Could not parse dataset id from href: ${href}`);
  return m[1];
}

/** Open a previously-saved dataset by id (loads from S3) and wait for render. */
export async function openSavedDataset(
  page: Page,
  ds: DatasetEntry,
  datasetId: string,
): Promise<MerfishHooks | null> {
  const seg = ds.dataType === "single_molecule" ? "sm-viewer" : "viewer";
  await page.goto(`/${seg}/${datasetId}`, { waitUntil: "domcontentloaded" });
  await waitForRender(page, ds);
  return readHooks(page);
}
