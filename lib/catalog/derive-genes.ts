import { gunzipSync } from "zlib";

const MANIFEST_MAX_BYTES = 5 * 1024 * 1024; // 5 MB cap on any single fetch
const FETCH_TIMEOUT_MS = 10_000;

/** SSRF guard: only https, and reject internal / private / metadata hosts. */
function isSafePublicUrl(u: URL): boolean {
  if (u.protocol !== "https:") return false;
  const host = u.hostname.toLowerCase();

  if (host === "localhost" || host.endsWith(".local")) return false;
  // Literal private / link-local / loopback IPv4 ranges + IPv6 loopback.
  if (
    /^127\./.test(host) ||
    /^10\./.test(host) ||
    /^192\.168\./.test(host) ||
    /^169\.254\./.test(host) ||
    /^172\.(1[6-9]|2\d|3[01])\./.test(host) ||
    host === "0.0.0.0" ||
    host === "::1"
  ) {
    return false;
  }

  return true;
}

/**
 * Fetch + parse a JSON file from a public S3 URL. Tries the URL as given, then
 * the gzipped variant (`<url>.gz`) — some datasets store manifests/indexes only
 * gzipped (plain `.json` 403s). SSRF-guarded, https-only, no-redirect, timed,
 * and size-capped. Returns `null` on any failure (never throws).
 */
export async function fetchJsonFromS3(url: string): Promise<any | null> {
  for (const target of [url, `${url}.gz`]) {
    try {
      const u = new URL(target);

      if (!isSafePublicUrl(u)) return null;

      const controller = new AbortController();
      const timeout = setTimeout(() => controller.abort(), FETCH_TIMEOUT_MS);
      const res = await fetch(u.toString(), {
        signal: controller.signal,
        redirect: "error", // don't follow redirects (SSRF)
      });

      clearTimeout(timeout);
      if (!res.ok) continue;

      const buf = Buffer.from(await res.arrayBuffer());

      if (buf.byteLength > MANIFEST_MAX_BYTES) return null;

      const text = target.endsWith(".gz")
        ? gunzipSync(buf).toString("utf-8")
        : buf.toString("utf-8");

      return JSON.parse(text);
    } catch {
      // try the next variant
    }
  }

  return null;
}

/** Strip trailing slashes so we can append a `/file.json` suffix cleanly. */
function normalizeBase(s3BaseUrl: string): string {
  return s3BaseUrl.replace(/\/+$/, "");
}

/**
 * Resolve the full gene list for a single catalog entry by reading its dataset
 * manifest from S3. Gene names live in different files per dataset type:
 *   - single_molecule → `manifest.json(.gz)` → `genes.unique_gene_names`
 *   - single_cell     → `expr/index.json(.gz)` → `genes[].name`
 *
 * Throws when the manifest/index can't be fetched or holds no usable gene list,
 * so callers can fail the import item rather than silently write empty genes.
 */
export async function deriveGenesForEntry({
  s3BaseUrl,
  datasetType,
}: {
  s3BaseUrl: string;
  datasetType: string;
}): Promise<string[]> {
  const base = normalizeBase(s3BaseUrl);

  if (datasetType === "single_molecule") {
    const manifest = await fetchJsonFromS3(`${base}/manifest.json`);
    const names = manifest?.genes?.unique_gene_names;

    if (!Array.isArray(names) || names.length === 0) {
      throw new Error(
        `single-molecule manifest at ${base} has no genes.unique_gene_names`,
      );
    }

    return names.map(String);
  }

  // single_cell (default)
  const index = await fetchJsonFromS3(`${base}/expr/index.json`);
  const genes = index?.genes;

  if (!Array.isArray(genes) || genes.length === 0) {
    throw new Error(`single-cell index at ${base}/expr/index.json has no genes`);
  }

  const names = genes
    .map((g: { name?: unknown }) => g?.name)
    .filter((n: unknown): n is string => typeof n === "string" && n.length > 0);

  if (names.length === 0) {
    throw new Error(`single-cell index at ${base}/expr/index.json has no gene names`);
  }

  return names;
}
