/**
 * S3 URL Parser Utility
 * Handles parsing and normalization of S3 URLs from various formats
 */

export interface ParsedS3Url {
  /** Base folder URL without trailing slash */
  baseUrl: string;
  /** AWS region extracted from URL */
  region: string;
  /** Bucket name */
  bucket: string;
  /** Path/prefix within bucket */
  prefix: string;
}

/**
 * Normalize S3 URL to base folder URL
 * Handles various input formats:
 * - https://bucket.s3.region.amazonaws.com/path/to/folder
 * - https://bucket.s3.region.amazonaws.com/path/to/folder/
 * - https://bucket.s3.region.amazonaws.com/path/to/folder/manifest.json.gz
 * - https://s3.region.amazonaws.com/bucket/path/to/folder
 *
 * @param url - S3 URL in any format
 * @returns Parsed S3 URL components with normalized base URL
 */
export function parseS3Url(url: string): ParsedS3Url {
  try {
    // Remove whitespace
    url = url.trim();

    // Parse URL
    const urlObj = new URL(url);
    const hostname = urlObj.hostname;
    const pathname = urlObj.pathname;

    let bucket: string;
    let region: string;
    let prefix: string;

    // Format 1: bucket.s3.region.amazonaws.com
    // Example: my-bucket.s3.us-east-1.amazonaws.com
    const virtualHostedMatch = hostname.match(/^(.+?)\.s3\.(.+?)\.amazonaws\.com$/);

    if (virtualHostedMatch) {
      bucket = virtualHostedMatch[1];
      region = virtualHostedMatch[2];
      prefix = pathname.substring(1); // Remove leading slash
    }
    // Format 2: s3.region.amazonaws.com/bucket/path
    // Example: s3.us-east-1.amazonaws.com/my-bucket/path
    else {
      const pathHostedMatch = hostname.match(/^s3\.(.+?)\.amazonaws\.com$/);

      if (pathHostedMatch) {
        region = pathHostedMatch[1];
        const pathParts = pathname.substring(1).split('/');
        bucket = pathParts[0];
        prefix = pathParts.slice(1).join('/');
      } else {
        throw new Error('Invalid S3 URL format. Expected format: https://bucket.s3.region.amazonaws.com/path or https://s3.region.amazonaws.com/bucket/path');
      }
    }

    // Remove trailing slash from prefix
    prefix = prefix.replace(/\/$/, '');

    // If URL ends with manifest.json.gz, strip it
    if (prefix.endsWith('/manifest.json.gz')) {
      prefix = prefix.slice(0, -'/manifest.json.gz'.length);
    } else if (prefix.endsWith('manifest.json.gz')) {
      // Handle case without leading slash
      prefix = prefix.slice(0, -'manifest.json.gz'.length).replace(/\/$/, '');
    }

    // Construct normalized base URL (virtual-hosted style)
    const baseUrl = `https://${bucket}.s3.${region}.amazonaws.com/${prefix}`;

    return {
      baseUrl,
      region,
      bucket,
      prefix,
    };
  } catch (error) {
    throw new Error(`Failed to parse S3 URL: ${error instanceof Error ? error.message : 'Invalid URL'}`);
  }
}

/**
 * Construct full S3 URL for a file within the dataset folder
 * @param baseUrl - Base folder URL from parseS3Url
 * @param filePath - Relative file path (e.g., "manifest.json.gz", "coords/spatial.bin.gz")
 * @returns Full S3 URL
 */
export function getS3FileUrl(baseUrl: string, filePath: string): string {
  // Remove leading slash from filePath if present
  filePath = filePath.replace(/^\//, '');

  // Ensure baseUrl doesn't have trailing slash
  baseUrl = baseUrl.replace(/\/$/, '');

  return `${baseUrl}/${filePath}`;
}

/**
 * Validate S3 URL format
 * @param url - URL to validate
 * @returns true if valid S3 URL
 */
export function isValidS3Url(url: string): boolean {
  try {
    parseS3Url(url);
    return true;
  } catch {
    return false;
  }
}

/**
 * Test if manifest exists at the given S3 URL
 * @param baseUrl - Base folder URL
 * @returns Promise<boolean> - true if manifest exists and is accessible
 */
export async function testManifestAccess(baseUrl: string): Promise<boolean> {
  try {
    const manifestUrl = getS3FileUrl(baseUrl, 'manifest.json.gz');
    const response = await fetch(manifestUrl, { method: 'HEAD' });
    return response.ok;
  } catch {
    return false;
  }
}
