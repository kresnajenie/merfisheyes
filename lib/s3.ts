// lib/s3.ts (or src/lib/s3.ts)
import {
  S3Client,
  GetObjectCommand,
  PutObjectCommand,
  CreateMultipartUploadCommand,
  UploadPartCommand,
  CompleteMultipartUploadCommand,
  AbortMultipartUploadCommand,
  ListMultipartUploadsCommand,
  ListObjectsV2Command,
  DeleteObjectsCommand,
} from "@aws-sdk/client-s3";
import { getSignedUrl } from "@aws-sdk/s3-request-presigner";
import { config } from "dotenv";
config({ path: ".env.local" });

// Optional S3-compatible endpoint (e.g. MinIO for local/E2E testing). When
// AWS_S3_ENDPOINT is set we switch to path-style addressing and pass explicit
// credentials; otherwise behaviour is unchanged (real AWS S3, default creds).
const S3_ENDPOINT = process.env.AWS_S3_ENDPOINT || process.env.S3_ENDPOINT || "";

// Initialize S3 client
export const s3Client = new S3Client({
  region: process.env.AWS_REGION || "us-east-1",
  // Credentials are automatically loaded from:
  // 1. Environment variables (AWS_ACCESS_KEY_ID, AWS_SECRET_ACCESS_KEY)
  // 2. IAM role (when running on EC2)
  // 3. AWS credentials file
  ...(S3_ENDPOINT
    ? {
        endpoint: S3_ENDPOINT,
        forcePathStyle: true,
        credentials: {
          accessKeyId: process.env.AWS_ACCESS_KEY_ID || "",
          secretAccessKey: process.env.AWS_SECRET_ACCESS_KEY || "",
        },
      }
    : {}),
});

export const S3_BUCKET =
  process.env.AWS_S3_BUCKET || process.env.S3_BUCKET || "";
export const AWS_REGION = process.env.AWS_REGION || "us-east-1";

// Helper to ensure bucket is configured
function ensureBucket() {
  if (!S3_BUCKET) {
    throw new Error("AWS_S3_BUCKET environment variable is required");
  }

  return S3_BUCKET;
}

/**
 * Generate a presigned URL for uploading a file to S3
 * @param key - S3 object key (path)
 * @param contentType - MIME type of the file
 * @param expiresIn - URL expiration time in seconds (default: 1 hour)
 */
export async function generatePresignedUploadUrl(
  key: string,
  contentType: string,
  expiresIn: number = 3600,
) {
  const command = new PutObjectCommand({
    Bucket: ensureBucket(),
    Key: key,
    ContentType: contentType,
    // Don't set ACL or other CORS-related parameters here
    // Let the bucket CORS policy handle it
  });

  const url = await getSignedUrl(s3Client, command, {
    expiresIn,
    // Don't sign additional headers that might conflict with CORS
    unhoistableHeaders: new Set(),
  });

  return {
    url,
    method: "PUT" as const,
    headers: {
      "Content-Type": contentType,
    },
  };
}

/**
 * Generate a presigned URL for downloading a file from S3
 * @param key - S3 object key (path)
 * @param expiresIn - URL expiration in seconds (default: 1 hour)
 */
export async function generatePresignedDownloadUrl(
  key: string,
  expiresIn: number = 3600,
): Promise<string> {
  const command = new GetObjectCommand({
    Bucket: ensureBucket(),
    Key: key,
  });

  const url = await getSignedUrl(s3Client, command, { expiresIn });

  return url;
}

/**
 * Generate presigned URLs for all files in a dataset
 * @param datasetId - The dataset ID (e.g., "ds_xyz123")
 * @param fileKeys - Array of file keys relative to dataset folder
 * @param expiresIn - URL expiration in seconds (default: 12 hours)
 */
export async function generateDatasetUrls(
  datasetId: string,
  fileKeys: string[],
  expiresIn: number = 12 * 3600, // 12 hours default
): Promise<Record<string, string>> {
  const urls: Record<string, string> = {};

  // Generate URLs in parallel for better performance
  const urlPromises = fileKeys.map(async (fileKey) => {
    const s3Key = `datasets/${datasetId}/${fileKey}`;
    const url = await generatePresignedDownloadUrl(s3Key, expiresIn);

    return { fileKey, url };
  });

  const results = await Promise.all(urlPromises);

  for (const { fileKey, url } of results) {
    urls[fileKey] = url;
  }

  return urls;
}

/**
 * Get the public URL for an S3 object (for public buckets only)
 * Note: This won't work if bucket is private - use generatePresignedDownloadUrl instead
 */
export function getPublicUrl(key: string): string {
  return `https://${ensureBucket()}.s3.${AWS_REGION}.amazonaws.com/${key}`;
}

/**
 * Generate presigned URL for manifest.json.gz
 * @param datasetId - The dataset ID
 * @param expiresIn - URL expiration in seconds (default: 12 hours)
 */
export async function generateManifestUrl(
  datasetId: string,
  expiresIn: number = 12 * 3600,
): Promise<string> {
  const manifestKey = `datasets/${datasetId}/manifest.json.gz`;

  return generatePresignedDownloadUrl(manifestKey, expiresIn);
}

// ─── Multipart upload (raw ingestion path) ──────────────────────
//
// Used by the server-side ingestion flow (app/api/ingest/*) for large raw
// uploads. Files above MULTIPART_THRESHOLD are uploaded as multiple presigned
// parts; smaller files keep the single-PUT path (generatePresignedUploadUrl).
// The S3 multipart UploadId is round-tripped through the client (no DB state);
// orphaned/incomplete uploads are cleaned up via listInProgressMultipartUploads
// + abortMultipartUpload (see /api/ingest/[id]/abort) and an S3 lifecycle rule.

/** Files larger than this use multipart upload. S3 single-PUT tops out ~5 GB. */
export const MULTIPART_THRESHOLD = 100 * 1024 * 1024; // 100 MB
/** Part size for multipart uploads. S3 requires every part except the last ≥ 5 MB. */
export const MULTIPART_PART_SIZE = 32 * 1024 * 1024; // 32 MB

/** Number of parts a file of `size` bytes splits into at MULTIPART_PART_SIZE. */
export function computePartCount(size: number): number {
  return Math.max(1, Math.ceil(size / MULTIPART_PART_SIZE));
}

/**
 * Start a multipart upload and return S3's UploadId.
 * @param key - full S3 object key (e.g. `raw/{datasetId}/{relativePath}`)
 * @param contentType - MIME type of the object
 */
export async function createMultipartUpload(
  key: string,
  contentType: string = "application/octet-stream",
): Promise<{ uploadId: string }> {
  const command = new CreateMultipartUploadCommand({
    Bucket: ensureBucket(),
    Key: key,
    ContentType: contentType,
  });

  const response = await s3Client.send(command);

  if (!response.UploadId) {
    throw new Error(`CreateMultipartUpload returned no UploadId for ${key}`);
  }

  return { uploadId: response.UploadId };
}

/**
 * Generate presigned PUT URLs for every part of a multipart upload.
 * @param key - full S3 object key
 * @param uploadId - S3 multipart UploadId from createMultipartUpload
 * @param partCount - number of parts (1-indexed part numbers)
 * @param expiresIn - URL expiration in seconds (default: 1 hour)
 */
export async function generatePresignedUploadPartUrls(
  key: string,
  uploadId: string,
  partCount: number,
  expiresIn: number = 3600,
): Promise<Array<{ partNumber: number; url: string }>> {
  const bucket = ensureBucket();

  const urlPromises = Array.from({ length: partCount }, async (_, i) => {
    const partNumber = i + 1;
    const command = new UploadPartCommand({
      Bucket: bucket,
      Key: key,
      UploadId: uploadId,
      PartNumber: partNumber,
    });
    const url = await getSignedUrl(s3Client, command, {
      expiresIn,
      unhoistableHeaders: new Set(),
    });

    return { partNumber, url };
  });

  return Promise.all(urlPromises);
}

/**
 * Finalize a multipart upload from the parts the client uploaded.
 * @param key - full S3 object key
 * @param uploadId - S3 multipart UploadId
 * @param parts - one entry per uploaded part with its S3-returned ETag
 */
export async function completeMultipartUpload(
  key: string,
  uploadId: string,
  parts: Array<{ partNumber: number; etag: string }>,
): Promise<void> {
  const command = new CompleteMultipartUploadCommand({
    Bucket: ensureBucket(),
    Key: key,
    UploadId: uploadId,
    MultipartUpload: {
      Parts: parts
        .slice()
        .sort((a, b) => a.partNumber - b.partNumber)
        .map((p) => ({ PartNumber: p.partNumber, ETag: p.etag })),
    },
  });

  await s3Client.send(command);
}

/**
 * Abort a single multipart upload, discarding any uploaded parts.
 */
export async function abortMultipartUpload(
  key: string,
  uploadId: string,
): Promise<void> {
  const command = new AbortMultipartUploadCommand({
    Bucket: ensureBucket(),
    Key: key,
    UploadId: uploadId,
  });

  await s3Client.send(command);
}

/**
 * List in-progress (incomplete) multipart uploads under a key prefix.
 * Used by /abort to find orphans to clean up without any client-held state.
 */
export async function listInProgressMultipartUploads(
  prefix: string,
): Promise<Array<{ key: string; uploadId: string }>> {
  const bucket = ensureBucket();
  const uploads: Array<{ key: string; uploadId: string }> = [];
  let keyMarker: string | undefined;
  let uploadIdMarker: string | undefined;

  do {
    const response = await s3Client.send(
      new ListMultipartUploadsCommand({
        Bucket: bucket,
        Prefix: prefix,
        KeyMarker: keyMarker,
        UploadIdMarker: uploadIdMarker,
      }),
    );

    for (const u of response.Uploads ?? []) {
      if (u.Key && u.UploadId) {
        uploads.push({ key: u.Key, uploadId: u.UploadId });
      }
    }

    if (response.IsTruncated) {
      keyMarker = response.NextKeyMarker;
      uploadIdMarker = response.NextUploadIdMarker;
    } else {
      keyMarker = undefined;
      uploadIdMarker = undefined;
    }
  } while (keyMarker || uploadIdMarker);

  return uploads;
}

/**
 * Delete every object under a key prefix (paginated list → batched delete).
 * @returns the number of objects deleted
 */
export async function deleteObjectsByPrefix(prefix: string): Promise<number> {
  const bucket = ensureBucket();
  let deleted = 0;
  let continuationToken: string | undefined;

  do {
    const listed = await s3Client.send(
      new ListObjectsV2Command({
        Bucket: bucket,
        Prefix: prefix,
        ContinuationToken: continuationToken,
      }),
    );

    const objects = (listed.Contents ?? [])
      .map((o) => o.Key)
      .filter((k): k is string => Boolean(k))
      .map((Key) => ({ Key }));

    if (objects.length > 0) {
      await s3Client.send(
        new DeleteObjectsCommand({
          Bucket: bucket,
          Delete: { Objects: objects, Quiet: true },
        }),
      );
      deleted += objects.length;
    }

    continuationToken = listed.IsTruncated
      ? listed.NextContinuationToken
      : undefined;
  } while (continuationToken);

  return deleted;
}
