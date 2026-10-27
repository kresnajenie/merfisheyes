import crypto from "node:crypto";
import fs from "node:fs";
import os from "node:os";
import path from "node:path";

import { GetObjectCommand, S3Client } from "@aws-sdk/client-s3";

import { test, expect } from "../../fixtures/merfish-fixtures";
import { signInDev } from "../../helpers/auth";

/**
 * A raw file over the multipart threshold (100 MB) goes up in parts, several
 * at once (lib/upload/uploadRawToS3.ts `uploadParts`), and the object S3
 * assembles is the file, byte for byte. Tagged @upload.
 */

const UPLOAD_ENABLED = process.env.E2E_UPLOAD === "1";

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

async function sha256Of(stream: NodeJS.ReadableStream): Promise<string> {
  const hash = crypto.createHash("sha256");

  for await (const chunk of stream) hash.update(chunk as Buffer);

  return hash.digest("hex");
}

test.describe("@upload a big raw file in parallel parts", () => {
  test.skip(!UPLOAD_ENABLED, "Set E2E_UPLOAD=1 (and start docker infra) to run upload tests.");

  let file: string;

  test.beforeAll(() => {
    // ~200 MB of molecules: 7 parts of 32 MB, more than go up at once
    file = path.join(fs.mkdtempSync(path.join(os.tmpdir(), "e2e-")), "spots.csv");
    const fd = fs.openSync(file, "w");
    const rows = Array.from(
      { length: 50_000 },
      (_, i) => `Gene${i % 300},${(i * 0.37).toFixed(3)},${(i * 0.73).toFixed(3)}\n`,
    ).join("");

    fs.writeSync(fd, "gene,x,y\n");
    for (let written = 0; written < 200 * 1024 * 1024; written += rows.length) {
      fs.writeSync(fd, rows);
    }
    fs.closeSync(fd);
  });
  test.afterAll(() => fs.rmSync(path.dirname(file), { recursive: true }));

  test("arrives whole, every part in its place", async ({ page }) => {
    test.setTimeout(300_000);
    const name = `e2e-multipart-${Date.now()}`;

    await signInDev(page, `e2e-multipart-${Date.now()}@example.com`);
    await page.goto("/", { waitUntil: "domcontentloaded" });
    await page.getByRole("button", { name: /single molecule/i }).first().click();
    await page.getByRole("radio", { name: "Upload & process on server" }).click();
    const initiated = page.waitForResponse("**/api/ingest/initiate");

    await page.locator("#sm-file-input-file").setInputFiles(file);
    await page.getByLabel("Dataset name").fill(name);
    await page.getByRole("button", { name: "Start upload" }).click();

    const { datasetId, uploadUrls } = await (await initiated).json();
    const target = uploadUrls["spots.csv"];

    expect(target.mode).toBe("multipart");
    expect(target.parts.length).toBeGreaterThan(6);

    // Every part up and the dataset queued for processing
    await expect
      .poll(
        async () =>
          (await (await page.request.get(`/api/ingest/${datasetId}/status`)).json())
            .status,
        { timeout: 240_000 },
      )
      .toBe("QUEUED");

    const object = await s3().send(
      new GetObjectCommand({
        Bucket: process.env.AWS_S3_BUCKET,
        Key: `raw/${datasetId}/spots.csv`,
      }),
    );

    expect(object.ContentLength).toBe(fs.statSync(file).size);
    expect(await sha256Of(object.Body as NodeJS.ReadableStream)).toBe(
      await sha256Of(fs.createReadStream(file)),
    );
  });
});
