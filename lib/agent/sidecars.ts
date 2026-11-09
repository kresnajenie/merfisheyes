import { LimitedFetchError, limitedFetch } from "@/lib/net/limited-fetch";

/**
 * The JSON files written beside a store — catalog.json, stats.json, qc.json
 * (docs/specs/schemas/, checked by scripts/validate_sidecars.py) — read over
 * public HTTPS and merged into one description of the dataset
 * (docs/specs/MCP-PROTOCOL.md §4). Three small GETs, no Zarr opened. A
 * missing or unreadable file is normal: legacy and non-BIL folders have none.
 */

const NAMES = ["catalog", "stats", "qc"] as const;

type Name = (typeof NAMES)[number];
type Doc = Record<string, unknown>;

export interface Sidecars {
  /** The folder the files were looked for in. */
  base: string;
  catalog: Doc | null;
  stats: Doc | null;
  qc: Doc | null;
}

/** Only S3 is asked: an s3BaseUrl can be registered by any signed-in user, and
 * the server must not be pointed at an arbitrary host. */
export const isS3Host = (host: string) =>
  /(^|\.)s3([.-][a-z0-9-]+)*\.amazonaws\.com$/i.test(host);

/** The folder holding a store: the parent of `…/<name>.zarr`, or the URL
 * itself when it is already a folder. */
export function sidecarBase(storeUrl: string): string | null {
  let u: URL;

  try {
    u = new URL(storeUrl);
  } catch {
    return null;
  }
  const path = u.pathname.replace(/\/+$/, "");
  const folder = /\.zarr$/i.test(path)
    ? path.slice(0, path.lastIndexOf("/"))
    : path;

  return `${u.protocol}//${u.host}${folder}/`;
}

const cache = new Map<string, { at: number; value: Doc | null }>();
const TTL_MS = 60 * 60 * 1000;
const MAX_ENTRIES = 500;

async function readOne(base: string, name: Name): Promise<Doc | null> {
  const url = `${base}${name}.json`;
  const hit = cache.get(url);

  if (hit && Date.now() - hit.at < TTL_MS) return hit.value;
  let value: Doc | null = null;

  try {
    const res = await limitedFetch(url, {
      maxBytes: 256 * 1024,
      timeoutMs: 5000,
      allowHost: isS3Host,
    });

    if (res.status === 200) {
      const parsed = JSON.parse(res.body);

      value =
        parsed && typeof parsed === "object" && !Array.isArray(parsed)
          ? parsed
          : null;
    }
  } catch (e) {
    if (
      !(e instanceof LimitedFetchError) &&
      !(e instanceof SyntaxError) &&
      !(e instanceof Error && e.name === "TimeoutError")
    ) {
      throw e;
    }
  }
  if (cache.size >= MAX_ENTRIES) cache.delete(cache.keys().next().value!);
  cache.set(url, { at: Date.now(), value });

  return value;
}

export async function readSidecars(storeUrl: string): Promise<Sidecars | null> {
  const base = sidecarBase(storeUrl);

  if (!base || !isS3Host(new URL(base).hostname)) return null;
  const [catalog, stats, qc] = await Promise.all(
    NAMES.map((n) => readOne(base, n)),
  );

  return { base, catalog, stats, qc };
}

const num = (v: unknown) =>
  typeof v === "number" && Number.isFinite(v) ? v : null;
const str = (v: unknown) => (typeof v === "string" && v ? v : null);

/** One flat description an agent can read at a glance. Absent values are
 * null, never guessed. */
export function summarise(s: Sidecars) {
  const { stats, qc, catalog } = s;
  const levels = (str(qc?.mmc_levels) ?? "").split(",").filter(Boolean);

  return {
    cells: num(stats?.cells) ?? num(qc?.n_cells),
    genes: num(stats?.genes),
    molecules: num(stats?.molecules),
    moleculeGenes: num(stats?.molecule_genes),
    sections: num(qc?.n_samples),
    dims: num(stats?.dims),
    images: Array.isArray(stats?.images) ? (stats.images as string[]) : [],
    deColumns: Array.isArray(stats?.de_stats)
      ? (stats.de_stats as string[])
      : [],
    input: str(stats?.input),
    pipeline: str(stats?.pipeline),
    qc: qc
      ? {
          medianCountsPerCell: num(qc.median_counts_per_cell),
          meanCountsPerCell: num(qc.mean_counts_per_cell),
          genesDetected: num(qc.n_genes_detected),
          maskPercentile: num(qc.mask_percentile),
          pctArtifact: num(qc.pct_artifact),
        }
      : null,
    taxonomy: levels.length
      ? {
          method: str(qc?.mmc_methods),
          levels: levels.map((level) => ({
            level,
            groups: num(qc?.[`n_${level}`]),
            meanBootstrap: num(qc?.[`mean_boot_${level}`]),
          })),
        }
      : null,
    investigator: str((catalog?.metadata as Doc | undefined)?.investigator),
    license: str((catalog?.metadata as Doc | undefined)?.license),
  };
}
