import fs from "node:fs";
import os from "node:os";
import path from "node:path";

import {
  test,
  expect,
  litFraction,
  REPO_ROOT,
  SELECTORS,
} from "../../fixtures/merfish-fixtures";
import { signInDev } from "../../helpers/auth";

/**
 * One folder holding cells, molecules and an image, uploaded for server
 * processing: each becomes its own dataset in one project, and the cells'
 * dataset is told to overlay the other two. Tagged @upload.
 *
 * The test stack has no processing queue, so the cells and molecules stay
 * queued; what is checked is everything that happens before that — and the
 * image, which needs no processing and is viewable straight away.
 */

const UPLOAD_ENABLED = process.env.E2E_UPLOAD === "1";
const TINY = path.join(REPO_ROOT, "tests", "data", "tiny");

test.describe("@upload combined upload", () => {
  test.skip(
    !UPLOAD_ENABLED,
    "Set E2E_UPLOAD=1 (and start docker infra) to run upload tests.",
  );
  // Under the one pipeline the folder is one dataset (prepare-upload.spec.ts)
  test.skip(
    process.env.NEXT_PUBLIC_INGEST_PIPELINE === "prepare",
    "The separate-datasets upload is not made under NEXT_PUBLIC_INGEST_PIPELINE=prepare.",
  );

  let folder: string;

  test.beforeAll(() => {
    folder = path.join(
      fs.mkdtempSync(path.join(os.tmpdir(), "e2e-")),
      "sample",
    );
    fs.mkdirSync(folder);
    fs.copyFileSync(
      path.join(TINY, "h5ad", "tiny.h5ad"),
      path.join(folder, "cells.h5ad"),
    );
    fs.copyFileSync(
      path.join(TINY, "single-molecule", "tiny.parquet"),
      path.join(folder, "transcripts.parquet"),
    );
    fs.cpSync(
      path.join(TINY, "zarr", "spatialdata_v3.zarr", "images", "dapi"),
      path.join(folder, "dapi.zarr"),
      { recursive: true },
    );
  });
  test.afterAll(() => fs.rmSync(path.dirname(folder), { recursive: true }));

  test("cells, molecules and an image from one folder", async ({ page }) => {
    const name = `e2e-combined-${Date.now()}`;
    const mine = async (type: string) =>
      (
        (await (await page.request.get(`/api/ingest/mine?type=${type}`)).json())
          .datasets as { id: string; title: string; status: string }[]
      ).filter((d) => d.title.startsWith(name));

    await signInDev(page, `e2e-combined-${Date.now()}@example.com`);

    await test.step("the folder is recognised as all three", async () => {
      await page.goto("/", { waitUntil: "domcontentloaded" });
      await page
        .getByRole("radio", { name: "Upload & process on server" })
        .click();
      await page.locator("#sc-file-input-folder").setInputFiles(folder);

      await expect(page.getByText("Name your dataset")).toBeVisible();
      await expect(
        page.getByText("Also upload the single-molecule data"),
      ).toBeVisible();
      await expect(page.getByText("Also upload the image")).toBeVisible();
      await expect(page.getByText("dapi.zarr")).toBeVisible();
    });

    await test.step("upload all three under one name", async () => {
      await page.getByLabel("Project name").fill(name);
      await page.getByRole("button", { name: "Upload all three" }).click();
      await page.waitForURL("**/explore");
      // The image needs no processing: its bar ends on a link to it
      await expect(page.getByText("Ready to view")).toBeVisible({
        timeout: 120_000,
      });
    });

    let imageId = "";
    let cellsId = "";

    await test.step("three datasets, grouped in one project", async () => {
      await expect
        .poll(async () => (await mine("single_cell")).length, {
          timeout: 60_000,
        })
        .toBe(1);
      const [image] = await mine("image");
      const [cells] = await mine("single_cell");
      const [molecules] = await mine("single_molecule");

      expect(image).toMatchObject({
        title: `${name} (image)`,
        status: "COMPLETE",
      });
      expect(cells.title).toBe(`${name} (cells)`);
      expect(molecules.title).toBe(`${name} (molecules)`);
      imageId = image.id;
      cellsId = cells.id;

      await expect
        .poll(async () => {
          const projects = (
            await (await page.request.get("/api/projects")).json()
          ).projects as { title: string; datasetCount?: number; id: string }[];

          return projects.filter((p) => p.title === name).length;
        })
        .toBe(1);

      // The cells' dataset names the other two as its overlays
      const mapping = await (
        await page.request.get(`/api/datasets/${cellsId}/mapping`)
      ).json();

      expect(mapping).toMatchObject({
        linkColumn: "__all__",
        links: { __all__: molecules.id },
        source: "app",
        image: imageId,
        imageSource: "app",
        imageDatasetId: imageId,
      });
    });

    await test.step("the uploaded image opens from the app's own storage", async () => {
      await page.goto(`/image-viewer/${imageId}`);
      await expect(page.locator(SELECTORS.imageCanvas)).toBeVisible();
      await expect
        .poll(() => litFraction(page, SELECTORS.imageCanvas), {
          message: "the image should be drawn",
          timeout: 60_000,
        })
        .toBeGreaterThan(0.2);
    });
  });

  test("molecules beside the cells are found by their columns, not their name", async ({
    page,
  }) => {
    // A lab's own naming: no "transcript" anywhere
    const named = path.join(
      fs.mkdtempSync(path.join(os.tmpdir(), "e2e-")),
      "lab",
    );

    fs.mkdirSync(named);
    fs.copyFileSync(
      path.join(TINY, "h5ad", "tiny.h5ad"),
      path.join(named, "cells.h5ad"),
    );
    fs.copyFileSync(
      path.join(TINY, "single-molecule", "tiny.parquet"),
      path.join(named, "set1__single_molecules.parquet"),
    );

    await signInDev(page, `e2e-combined-named-${Date.now()}@example.com`);
    await page.goto("/", { waitUntil: "domcontentloaded" });
    await page
      .getByRole("radio", { name: "Upload & process on server" })
      .click();
    await page.locator("#sc-file-input-folder").setInputFiles(named);

    await expect(page.getByText("Name your dataset")).toBeVisible();
    await expect(
      page.getByText("Also upload the single-molecule data"),
    ).toBeVisible();
    fs.rmSync(path.dirname(named), { recursive: true });
  });
});
