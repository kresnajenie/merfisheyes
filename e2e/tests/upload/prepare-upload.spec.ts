import { execFileSync } from "node:child_process";
import fs from "node:fs";
import os from "node:os";
import path from "node:path";

import {
  GetObjectCommand,
  ListObjectsV2Command,
  S3Client,
} from "@aws-sdk/client-s3";
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

/** The axes of the points element the worker wrote. */
async function pointsAxes(id: string): Promise<string[]> {
  const res = await s3().send(
    new GetObjectCommand({
      Bucket: process.env.AWS_S3_BUCKET,
      Key: `datasets/${id}/${id}.zarr/points/molecules/zarr.json`,
    }),
  );

  return JSON.parse(await res.Body!.transformToString()).attributes.axes;
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
      // The molecules have a z column: offered, and off unless ticked
      await expect(
        page.getByRole("checkbox", { name: /Keep the molecules' Z/ }),
      ).not.toBeChecked();
      // No separate molecules or image dataset is offered
      await expect(
        page.getByText("Also upload the single-molecule data"),
      ).toHaveCount(0);
    });

    await test.step("one dataset, the image left behind", async () => {
      await page.getByLabel("Dataset name").fill(name);
      const initiate = page.waitForRequest("**/api/ingest/initiate");

      await page.getByRole("button", { name: "Start upload" }).click();
      const params = (await initiate).postDataJSON().processingParams;

      expect(params.pipeline).toBe("prepare");
      expect(params.keepZ).toBeUndefined();
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
    await page.getByRole("checkbox", { name: /Keep the molecules' Z/ }).check();
    await page.getByLabel("Dataset name").fill(name);
    const initiate = page.waitForRequest("**/api/ingest/initiate");

    await page.getByRole("button", { name: "Start upload" }).click();
    // What the app stored is what Batch would hand the worker
    const params = (await initiate).postDataJSON().processingParams;

    expect(params).toMatchObject({ pipeline: "prepare", keepZ: true });
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
          `PROCESSING_PARAMS=${JSON.stringify(params)}`,
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
      // Z was ticked: the molecules kept it
      expect(await pointsAxes(id)).toEqual(["x", "y", "z"]);
    });

    await test.step("the upload bar links its viewer, and a cell-viewer link opens it too", async () => {
      const status = await (
        await page.request.get(`/api/ingest/${id}/status`)
      ).json();

      expect(status.viewerUrl).toBe(`/spatialdata-viewer/${id}`);
      // The upload bar gave /viewer links before; one already shared opens
      await page.goto(`/viewer/${id}`);
      await page.waitForFunction(
        () => (window as any).__merfish?.stats?.pointCount === 300,
        undefined,
        { timeout: 120_000 },
      );
      await expect(page.locator(SELECTORS.smLayerButton)).toBeVisible();
    });
  });

  test("a molecule table alone becomes a store that opens", async ({ page }) => {
    test.skip(
      !WORKER_IMAGE || !process.env.CALLBACK_SECRET,
      "Set E2E_WORKER_IMAGE (and CALLBACK_SECRET, as the app has it) to run the worker.",
    );
    test.setTimeout(600_000);
    const name = `e2e-prepare-sm-${Date.now()}`;
    // Its own molecules, so it is never a duplicate of an earlier run's
    const table = path.join(fs.mkdtempSync(path.join(os.tmpdir(), "e2e-")), "spots.csv");
    const rows = Array.from(
      { length: 2000 },
      (_, i) => `gene_${i % 10},${(Math.random() * 500).toFixed(2)},${(Math.random() * 500).toFixed(2)}`,
    );

    fs.writeFileSync(table, ["feature_name,x_location,y_location", ...rows].join("\n"));

    await signInDev(page, `e2e-prepare-sm-${Date.now()}@example.com`);
    await page.goto("/", { waitUntil: "domcontentloaded" });
    await page.getByRole("button", { name: /single molecule/i }).first().click();
    await page.getByRole("radio", { name: "Upload & process on server" }).click();
    await page.locator("#sm-file-input-file").setInputFiles(table);
    await page.getByLabel("Dataset name").fill(name);
    const initiate = page.waitForRequest("**/api/ingest/initiate");

    await page.getByRole("button", { name: "Start upload" }).click();
    const params = (await initiate).postDataJSON().processingParams;

    expect(params).toMatchObject({ kind: "single_molecule", pipeline: "prepare" });
    await page.waitForURL("**/explore");

    let id = "";

    await expect
      .poll(
        async () =>
          (await mine(page, "single_molecule", name)).find((d) => d.status === "QUEUED")?.id ?? "",
        { timeout: 120_000 },
      )
      .toMatch(/^sm_/);
    id = (await mine(page, "single_molecule", name))[0].id;

    execFileSync(
      "docker",
      [
        "run", "--rm", "--network", "host",
        "-e", `DATASET_ID=${id}`,
        "-e", "DATASET_KIND=single_molecule",
        "-e", `PROCESSING_PARAMS=${JSON.stringify(params)}`,
        "-e", `AWS_S3_ENDPOINT=${process.env.AWS_S3_ENDPOINT}`,
        "-e", `AWS_S3_BUCKET=${process.env.AWS_S3_BUCKET}`,
        "-e", `AWS_ACCESS_KEY_ID=${process.env.AWS_ACCESS_KEY_ID}`,
        "-e", `AWS_SECRET_ACCESS_KEY=${process.env.AWS_SECRET_ACCESS_KEY}`,
        "-e", `CALLBACK_URL=${new URL(page.url()).origin}/api/ingest/${id}/callback`,
        "-e", `CALLBACK_SECRET=${process.env.CALLBACK_SECRET}`,
        WORKER_IMAGE!,
      ],
      { stdio: "inherit", timeout: 300_000 },
    );
    fs.rmSync(path.dirname(table), { recursive: true });

    // Still its sm_ id, now a SpatialData store, opened by that id
    await expect
      .poll(async () => (await mine(page, "spatialdata", name)).length, {
        timeout: 120_000,
      })
      .toBe(1);
    const [store] = await mine(page, "spatialdata", name);

    expect(store).toMatchObject({
      id,
      status: "COMPLETE",
      viewerUrl: `/spatialdata-viewer/${id}`,
    });
    await page.goto(store.viewerUrl!);
    await expect(page.locator(SELECTORS.smCanvas).first()).toBeVisible({
      timeout: 120_000,
    });
    await expect(page.getByText(/Dataset fetch failed|Invalid dataset ID/)).toHaveCount(0);
  });

  test("a molecules store overlays a cells store, attached in the dialog", async ({
    page,
  }) => {
    test.skip(
      !WORKER_IMAGE || !process.env.CALLBACK_SECRET,
      "Set E2E_WORKER_IMAGE (and CALLBACK_SECRET, as the app has it) to run the worker.",
    );
    test.setTimeout(900_000);
    const run = Date.now();
    const work = fs.mkdtempSync(path.join(os.tmpdir(), "e2e-"));

    /** Upload on a card, run the worker image on it as Batch would, wait for the store. */
    async function uploadAndProcess(
      card: "cells" | "molecules",
      input: string,
      name: string,
      tick: RegExp[] = [],
    ): Promise<string> {
      await page.goto("/", { waitUntil: "domcontentloaded" });
      await page
        .getByRole("button", {
          name: card === "cells" ? /single cell/i : /single molecule/i,
        })
        .first()
        .click();
      await page.getByRole("radio", { name: "Upload & process on server" }).click();
      await page
        .locator(card === "cells" ? "#sc-file-input-folder" : "#sm-file-input-file")
        .setInputFiles(input);
      for (const box of tick) await page.getByRole("checkbox", { name: box }).check();
      await page.getByLabel("Dataset name").fill(name);
      const initiate = page.waitForRequest("**/api/ingest/initiate");

      await page.getByRole("button", { name: "Start upload" }).click();
      const params = (await initiate).postDataJSON().processingParams;
      const type = card === "cells" ? "single_cell" : "single_molecule";

      await page.waitForURL("**/explore");
      await expect
        .poll(async () => (await mine(page, type, name))[0]?.status ?? "", {
          timeout: 120_000,
        })
        .toBe("QUEUED");
      const [{ id }] = await mine(page, type, name);

      execFileSync(
        "docker",
        [
          "run", "--rm", "--network", "host",
          "-e", `DATASET_ID=${id}`,
          "-e", `DATASET_KIND=${type}`,
          "-e", `PROCESSING_PARAMS=${JSON.stringify(params)}`,
          "-e", `AWS_S3_ENDPOINT=${process.env.AWS_S3_ENDPOINT}`,
          "-e", `AWS_S3_BUCKET=${process.env.AWS_S3_BUCKET}`,
          "-e", `AWS_ACCESS_KEY_ID=${process.env.AWS_ACCESS_KEY_ID}`,
          "-e", `AWS_SECRET_ACCESS_KEY=${process.env.AWS_SECRET_ACCESS_KEY}`,
          "-e", `CALLBACK_URL=${new URL(page.url()).origin}/api/ingest/${id}/callback`,
          "-e", `CALLBACK_SECRET=${process.env.CALLBACK_SECRET}`,
          WORKER_IMAGE!,
        ],
        { stdio: "inherit", timeout: 300_000 },
      );
      await expect
        .poll(async () => (await mine(page, "spatialdata", name)).length, {
          timeout: 120_000,
        })
        .toBe(1);

      return id;
    }

    await signInDev(page, `e2e-prepare-overlay-${run}@example.com`);

    // The molecules, uploaded alone
    const table = path.join(work, "spots.csv");

    fs.writeFileSync(
      table,
      [
        "feature_name,x_location,y_location",
        ...Array.from(
          { length: 2000 },
          (_, i) =>
            `gene_${i % 10},${(Math.random() * 500).toFixed(2)},${(Math.random() * 500).toFixed(2)}`,
        ),
      ].join("\n"),
    );
    const smName = `e2e-overlay-molecules-${run}`;
    const smId = await uploadAndProcess("molecules", table, smName);

    // The cells, uploaded later with an image and no molecules of their own
    const cells = path.join(work, "cells");

    fs.mkdirSync(cells);
    fs.copyFileSync(path.join(TINY, "h5ad", "tiny.h5ad"), path.join(cells, "cells.h5ad"));
    fs.cpSync(
      path.join(TINY, "zarr", "spatialdata_v3.zarr", "images", "dapi"),
      path.join(cells, "dapi.zarr"),
      { recursive: true },
    );
    // A run's own files, so it is never a duplicate of an earlier run's
    fs.writeFileSync(path.join(cells, "dapi.zarr", "e2e-run.txt"), String(run));
    const scId = await uploadAndProcess("cells", cells, `e2e-overlay-cells-${run}`, [
      /Also upload the images/,
    ]);

    fs.rmSync(work, { recursive: true });

    await test.step("the molecules store is offered and attached in the dialog", async () => {
      await page.goto(`/spatialdata-viewer/${scId}`);
      await page.waitForFunction(
        () => (window as any).__merfish?.stats?.pointCount === 300,
        undefined,
        { timeout: 120_000 },
      );
      // No molecules yet: the store has none of its own
      await expect(page.locator(SELECTORS.smLayerButton)).toHaveCount(0);
      await page.locator(SELECTORS.overlayButton).click();
      await page.getByPlaceholder("Search your molecule datasets…").fill(smName);
      await page.getByRole("button", { name: new RegExp(smName) }).click();
      // The viewer reopens with the overlay
      await page.waitForLoadState("load");
      await page.waitForFunction(
        () => (window as any).__merfish?.stats?.pointCount === 300,
        undefined,
        { timeout: 120_000 },
      );
      await expect(page.locator(SELECTORS.smLayerButton)).toBeVisible({
        timeout: 60_000,
      });
      // Drawn by the tiled points layer, which reports what it has drawn
      await expect(page.getByText(/[\d.]+ ?k? drawn/).first()).toBeVisible({
        timeout: 60_000,
      });
    });

    await test.step("mapping.json links it as a store and keeps the store's own entry", async () => {
      const mapping = await (
        await page.request.get(`/api/datasets/${scId}/mapping`)
      ).json();

      expect(mapping).toMatchObject({
        cells: `${scId}.zarr`,
        linkColumn: "__all__",
        links: { __all__: smId },
        source: "app",
        datasetId: smId,
        linkFormat: "store",
      });
    });

    await test.step("a cells store is not offered as molecules", async () => {
      const res = await page.request.post(`/api/datasets/${smId}/overlay`, {
        data: { smDatasetId: scId },
      });

      // The molecules store has no table, so it takes no overlay — and the
      // cells store has no molecules to give
      expect(res.status()).toBe(400);
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
