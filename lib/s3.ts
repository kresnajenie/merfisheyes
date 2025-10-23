// lib/s3.ts (or src/lib/s3.ts)
import {
  S3Client,
  GetObjectCommand,
  PutObjectCommand,
} from "@aws-sdk/client-s3";
import { getSignedUrl } from "@aws-sdk/s3-request-presigner";
import { config } from "dotenv";
config({ path: ".env.local" });

// Initialize S3 client
export const s3Client = new S3Client({
  region: process.env.AWS_REGION || "us-east-1",
  // Credentials are automatically loaded from:
  // 1. Environment variables (AWS_ACCESS_KEY_ID, AWS_SECRET_ACCESS_KEY)
  // 2. IAM role (when running on EC2)
  // 3. AWS credentials file
});

export const S3_BUCKET = process.env.AWS_S3_BUCKET || process.env.S3_BUCKET || "";
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
