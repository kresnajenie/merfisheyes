import { execSync } from "node:child_process";

import {
  S3Client,
  HeadBucketCommand,
  CreateBucketCommand,
} from "@aws-sdk/client-s3";

/**
 * Global setup for the @upload suite. No-op unless E2E_UPLOAD=1. Prepares the
 * hermetic infrastructure started by docker-compose.test.yml:
 *   1. applies the Prisma schema to the test database
 *   2. ensures the S3 bucket exists (MinIO)
 *
 * Env is loaded from .env.test by playwright.config.ts before this runs.
 */
export default async function globalSetup(): Promise<void> {
  if (process.env.E2E_UPLOAD !== "1") return;

  console.log("[global-setup] preparing upload infra (Prisma migrate + S3 bucket)");

  execSync("npx prisma migrate deploy", { stdio: "inherit", env: process.env });

  const bucket = process.env.AWS_S3_BUCKET;
  if (!bucket) throw new Error("[global-setup] AWS_S3_BUCKET is required for upload tests");

  const s3 = new S3Client({
    region: process.env.AWS_REGION || "us-east-1",
    endpoint: process.env.AWS_S3_ENDPOINT,
    forcePathStyle: true,
    credentials: {
      accessKeyId: process.env.AWS_ACCESS_KEY_ID || "",
      secretAccessKey: process.env.AWS_SECRET_ACCESS_KEY || "",
    },
  });

  try {
    await s3.send(new HeadBucketCommand({ Bucket: bucket }));
    console.log(`[global-setup] bucket '${bucket}' present`);
  } catch {
    await s3.send(new CreateBucketCommand({ Bucket: bucket }));
    console.log(`[global-setup] created bucket '${bucket}'`);
  }
}
