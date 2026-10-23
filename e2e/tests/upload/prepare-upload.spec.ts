import { execFileSync } from "node:child_process";
import fs from "node:fs";
import os from "node:os";
import path from "node:path";

import { ListObjectsV2Command, S3Client } from "@aws-sdk/client-s3";
import type { Page } from "@playwright/test";

import {
  test,
  expect,
  REPO_ROOT,
  SELECTORS,
} from "../../fixtures/merfish-fixtures";
import { signInDev } from "../../helpers/auth";

/**
 * "Upload & process on server" under the one pipeline
 * (NEXT_PUBLIC_INGEST_PIPELINE=prepare, docs/plans/ONE-PIPELINE.md): a folder
 * of cells, molecules and an image becomes ONE dataset — the cells and
 * molecules always sent, the image only when ticked. Tagged @upload; runs
 * only when the app was started with the setting.
 *
 * With E2E_WORKER_IMAGE (a local build of worker/Dockerfile) and
 * CALLBACK_SECRET (the app's), the second test also runs that image on the
 * upload against the test MinIO, as AWS Batch would, signs its callback to
 * this app, and opens the store it wrote.
 */

const UPLOAD_ENABLED = process.env.E2E_UPLOAD === "1";
const PREPARE = process.env.NEXT_PUBLIC_INGEST_PIPELINE === "prepare";
const WORKER_IMAGE = process.env.E2E_WORKER_IMAGE;
const TINY = path.join(REPO_ROOT, "tests", "data", "tiny");

const s3 = () =>
  new S3Client({
    endpoint: process.env.AWS_S3_ENDPOINT,
    region: process.env.AWS_REGION || "us-east-1",
    forcePathStyle: true,
    credentials: {
      accessKeyId: process.env.AWS_ACCESS_KEY_ID ?? "",
      secretAccessKey: process.env.AWS_SECRET_ACCESS_KEY ?? "",
    },
  });

/** The keys uploaded under raw/<id>/, relative to it. */
async function rawKeys(id: string): Promise<string[]> {
  const res = await s3().send(
    new ListObjectsV2Command({
      Bucket: process.env.AWS_S3_BUCKET,
      Prefix: `raw/${id}/`,
    }),
  );

  return (res.Contents ?? [])
    .map((o) => o.Key!.slice(`raw/${id}/`.length))
    .sort();
}

const mine = async (page: Page, type: string, name: string) =>
  (
    (await (await page.request.get(`/api/ingest/mine?type=${type}`)).json())
      .datasets as {
      id: string;
      title: string;
      status: string;
      viewerUrl: string | null;
    }[]
  ).filter((d) => d.title === name);

test.describe("@upload one store from the server", () => {
  test.skip(
    !UPLOAD_ENABLED || !PREPARE,
    "Set E2E_UPLOAD=1 and NEXT_PUBLIC_INGEST_PIPELINE=prepare to run.",
  );

  let folder: string;

  test.beforeAll(() => {
    folder = path.join(fs.mkdtempSync(path.join(os.tmpdir(), "e2e-")), "lab");
    fs.mkdirSync(folder);
    fs.copyFileSync(
      path.join(TINY, "h5ad", "tiny.h5ad"),
      path.join(folder, "cells.h5ad"),
    );
    // A lab's own name: found by its columns
    fs.copyFileSync(
      path.join(TINY, "single-molecule", "tiny.parquet"),
      path.join(folder, "set1__single_molecules.parquet"),
    );
    fs.cpSync(
      path.join(TINY, "zarr", "spatialdata_v3.zarr", "images", "dapi"),
      path.join(folder, "dapi.zarr"),
      { recursive: true },
    );
  });
  test.afterAll(() => fs.rmSync(path.dirname(folder), { recursive: true }));

  async function dropFolder(page: Page, user: string) {
    await signInDev(page, `${user}-${Date.now()}@example.com`);
    await page.goto("/", { waitUntil: "domcontentloaded" });
    await page.getByRole("radio", { name: "Upload & process on server" }).click();
    await page.locator("#sc-file-input-folder").setInputFiles(folder);
    await expect(page.getByText("Name your dataset")).toBeVisible();
  }

  test("cells and molecules go up as one dataset; the image waits to be asked", async ({
    page,
  }) => {
    const name = `e2e-prepare-${Date.now()}`;

    await dropFolder(page, "e2e-prepare");

    await test.step("the dialog lists what goes into the store", async () => {
      await expect(page.getByText("into one SpatialData store")).toBeVisible();
      await expect(
        page.getByText("set1__single_molecules.parquet"),
      ).toBeVisible();
      const images = page.getByRole("checkbox", {
        name: /Also upload the images/,
      });

      await expect(images).not.toBeChecked();
      // No separate molecules or image dataset is offered
      await expect(
        page.getByText("Also upload the single-molecule data"),
      ).toHaveCount(0);
    });

    await test.step("one dataset, the image left behind", async () => {
      await page.getByLabel("Dataset name").fill(name);
      await page.getByRole("button", { name: "Start upload" }).click();
      await page.waitForURL("**/explore");

      await expect
        .poll(async () => (await mine(page, "single_cell", name)).length, {
          timeout: 60_000,
        })
        .toBe(1);
      const [cells] = await mine(page, "single_cell", name);

      expect(await mine(page, "single_molecule", name)).toHaveLength(0);
      expect(await mine(page, "image", name)).toHaveLength(0);
      expect(await rawKeys(cells.id)).toEqual([
        "cells.h5ad",
        "set1__single_molecules.parquet",
      ]);
    });
  });

  test("ticked, the image goes too; the worker makes one store that opens", async ({
    page,
  }) => {
    test.skip(
      !WORKER_IMAGE || !process.env.CALLBACK_SECRET,
      "Set E2E_WORKER_IMAGE (and CALLBACK_SECRET, as the app has it) to run the worker.",
    );
    test.setTimeout(600_000);
    const name = `e2e-prepare-worker-${Date.now()}`;

    // A dataset is fingerprinted by its raw files, and a re-upload of the same
    // ones is turned away as a duplicate: this run's are its own
    fs.writeFileSync(path.join(folder, "dapi.zarr", "e2e-run.txt"), name);
    await dropFolder(page, "e2e-prepare-worker");
    await page.getByRole("checkbox", { name: /Also upload the images/ }).check();
    await page.getByLabel("Dataset name").fill(name);
    await page.getByRole("button", { name: "Start upload" }).click();
    await page.waitForURL("**/explore");

    let id = "";

    await expect
      .poll(
        async () =>
          (await mine(page, "single_cell", name)).find(
            (d) => d.status === "QUEUED",
          )?.id ?? "",
        { timeout: 120_000 },
      )
      .not.toBe("");
    id = (await mine(page, "single_cell", name))[0].id;
    expect(await rawKeys(id)).toContain("dapi.zarr/zarr.json");

    await test.step("the worker image runs on it, as Batch would", async () => {
      const base = new URL(page.url()).origin;

      execFileSync(
        "docker",
        [
          "run",
          "--rm",
          "--network",
          "host",
          "-e",
          `DATASET_ID=${id}`,
          "-e",
          "DATASET_KIND=single_cell",
          "-e",
          `PROCESSING_PARAMS=${JSON.stringify({ kind: "single_cell", pipeline: "prepare", stages: { chunk: { chunkSize: 1 } } })}`,
          "-e",
          `AWS_S3_ENDPOINT=${process.env.AWS_S3_ENDPOINT}`,
          "-e",
          `AWS_S3_BUCKET=${process.env.AWS_S3_BUCKET}`,
          "-e",
          `AWS_ACCESS_KEY_ID=${process.env.AWS_ACCESS_KEY_ID}`,
          "-e",
          `AWS_SECRET_ACCESS_KEY=${process.env.AWS_SECRET_ACCESS_KEY}`,
          "-e",
          `CALLBACK_URL=${base}/api/ingest/${id}/callback`,
          "-e",
          `CALLBACK_SECRET=${process.env.CALLBACK_SECRET}`,
          WORKER_IMAGE!,
        ],
        { stdio: "inherit", timeout: 300_000 },
      );
    });

    await test.step("it is a SpatialData dataset and opens with all three", async () => {
      // The callback finishes it; on a cold dev server that takes a moment
      await expect
        .poll(async () => (await mine(page, "spatialdata", name)).length, {
          timeout: 120_000,
        })
        .toBe(1);
      const [store] = await mine(page, "spatialdata", name);

      expect(store).toMatchObject({
        status: "COMPLETE",
        viewerUrl: `/spatialdata-viewer/${id}`,
      });
      await page.goto(store.viewerUrl!);
      await expect(page.locator(SELECTORS.scCanvas).first()).toBeVisible({
        timeout: 120_000,
      });
      await page.waitForFunction(
        () => (window as any).__merfish?.stats?.pointCount === 300,
        undefined,
        { timeout: 120_000 },
      );
      await expect(page.locator(SELECTORS.imageButton)).toBeVisible();
      await expect(page.locator(SELECTORS.smLayerButton)).toBeVisible();
    });
  });
});

test.describe("@upload molecule columns on the server card", () => {
  test.skip(
    !UPLOAD_ENABLED,
    "Set E2E_UPLOAD=1 (and start docker infra) to run upload tests.",
  );

  test("a table whose columns are not recognised asks for them", async ({
    page,
  }) => {
    const dir = fs.mkdtempSync(path.join(os.tmpdir(), "e2e-"));
    const table = path.join(dir, "spots.csv");

    fs.writeFileSync(table, "probe,px,py\nGad1,1.5,2.5\nSst,3.0,4.0\n");

    await signInDev(page, `e2e-columns-${Date.now()}@example.com`);
    await page.goto("/", { waitUntil: "domcontentloaded" });
    await page
      .getByRole("button", { name: /single molecule/i })
      .first()
      .click();
    await page.getByRole("radio", { name: "Upload & process on server" }).click();
    await page.locator("#sm-file-input-file").setInputFiles(table);

    // The columns come from the file's own header, so the list is never blank
    await expect(page.getByText("Confirm columns")).toBeVisible();
    await page.getByLabel(/Gene \/ feature/).selectOption("probe");
    await page.getByLabel(/X coordinate/).selectOption("px");
    await page.getByLabel(/Y coordinate/).selectOption("py");
    await page.getByRole("button", { name: "Looks good — process" }).click();

    await expect(page.getByText("Name your dataset")).toBeVisible();
    const initiate = page.waitForRequest("**/api/ingest/initiate");

    await page.getByLabel("Dataset name").fill(`e2e-columns-${Date.now()}`);
    await page.getByRole("button", { name: "Start upload" }).click();

    expect((await initiate).postDataJSON().processingParams).toMatchObject({
      kind: "single_molecule",
      columnMapping: { gene: "probe", x: "px", y: "py" },
    });
    fs.rmSync(dir, { recursive: true });
  });
});
