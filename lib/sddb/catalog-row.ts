/**
 * How one row of a SpatialData-DB catalog release becomes one Explore card:
 * pure functions, so the sync script and its tests share them exactly.
 * Spec: docs/specs/sddb-catalog.md.
 */
import type { CatalogRow } from "./catalog-schema";

import { CONVENTIONS } from "@/lib/spatialdata/conventions";

/**
 * Where each public bucket lives (`s3://` URLs carry no region), and the
 * spellings Explore filters on: docs/specs/conventions.json.
 */
const BUCKET_REGION = CONVENTIONS.catalog.bucketRegions;
const SPECIES = CONVENTIONS.catalog.species;
const PLATFORM = CONVENTIONS.catalog.platforms;

/**
 * `s3://bucket/key.zarr` → the virtual-hosted HTTPS URL the browser reads.
 * Throws on a bucket whose region is unknown rather than guessing one.
 */
export function storeUrlFromZarrUrl(zarrUrl: string): string {
  const m = /^s3:\/\/([^/]+)\/(.+?)\/*$/.exec(zarrUrl);

  if (!m) throw new Error(`not an s3:// URL: ${zarrUrl}`);
  const [, bucket, key] = m;
  const region =
    BUCKET_REGION[bucket] ?? /-([a-z]{2}-[a-z]+-\d)$/.exec(bucket)?.[1];

  if (!region) throw new Error(`unknown region for bucket ${bucket}`);

  return `https://${bucket}.s3.${region}.amazonaws.com/${key}`;
}

const str = (v: unknown): string | null => {
  if (v === null || v === undefined) return null;
  const s = String(v).trim();

  return s === "" || s === "None" ? null : s;
};
const num = (v: unknown): number | null => {
  if (v === null || v === undefined) return null;
  const n = Number(v);

  return Number.isFinite(n) ? n : null;
};
const bool = (v: unknown): boolean | null =>
  typeof v === "boolean" ? v : null;
const capitalize = (s: string) => s.charAt(0).toUpperCase() + s.slice(1);
const cap = (s: string, n: number) => (s.length > n ? s.slice(0, n) : s);

export function formatBytes(n: number): string {
  if (n >= 1e12) return `${(n / 1e12).toFixed(1)} TB`;
  if (n >= 1e9) return `${(n / 1e9).toFixed(1)} GB`;
  if (n >= 1e6) return `${(n / 1e6).toFixed(1)} MB`;

  return `${Math.round(n / 1e3)} KB`;
}

export interface MappedCatalogRow {
  uid: string;
  /** Upsert payload for `CatalogDataset` (minus the audit columns). */
  data: {
    title: string;
    description: string;
    species: string | null;
    disease: string | null;
    tissue: string | null;
    platform: string;
    tags: string[];
    thumbnailUrl: string | null;
    metadata: Record<string, unknown>;
    externalLink: string | null;
    publicationLink: string | null;
    numCells: number | null;
    numGenes: number | null;
    source: "sddb";
    sddbUid: string;
    license: string | null;
    licenseNoncommercial: boolean | null;
    sourceRelease: string;
    sourceGeneratedAt: string | null;
    isPublished: true;
    isBil: false;
    isCommunity: false;
    isInternal: false;
  };
  /** The row's one entry: the store, opened by its public URL. */
  entry: {
    label: string;
    datasetType: "spatialdata";
    s3BaseUrl: string;
  };
}

/**
 * One catalog row → one card. `release` is the GitHub release tag the file
 * came from and `generatedAt` the file's `sddb_generated_at` metadata; both
 * are recorded on the row so a later run can tell what it is overwriting.
 */
export function mapCatalogRow(
  row: CatalogRow,
  { release, generatedAt }: { release: string; generatedAt: string | null },
): MappedCatalogRow {
  const uid = str(row.uid);

  if (!uid) throw new Error("catalog row without uid");
  const zarrUrl = str(row.zarr_url);

  if (!zarrUrl) throw new Error(`catalog row ${uid} without zarr_url`);

  const technology = str(row.technology) ?? "unknown";
  const platform = PLATFORM[technology] ?? technology;
  const organism = str(row.organism);
  const species = organism
    ? (SPECIES[organism.toLowerCase()] ?? capitalize(organism))
    : null;
  const tissue = str(row.tissue);
  const disease = str(row.disease);
  const assay = str(row.assay) ?? platform;
  const nObs = num(row.n_obs);
  const nFeatures = num(row.n_features);
  const doi = str(row.doi);
  const license = str(row.license);
  const sizeBytes = num(row.size_bytes);
  const nTranscripts = num(row.n_transcripts);
  const collections = Array.isArray(row.collections)
    ? (row.collections as unknown[]).map(String).filter(Boolean)
    : [];

  const parts = [assay, organism ? `${organism} ${tissue ?? ""}`.trim() : null]
    .filter(Boolean)
    .join(" of ");
  const description = cap(
    [
      parts + (disease ? `, ${disease}` : "") + ".",
      nObs != null && nFeatures != null
        ? `${nObs.toLocaleString("en-US")} observations × ${nFeatures.toLocaleString("en-US")} features.`
        : null,
      "From SpatialData-DB" +
        (str(row.tier) ? ` (${str(row.tier)} tier)` : "") +
        ".",
    ]
      .filter(Boolean)
      .join(" "),
    2000,
  );

  const metadata: Record<string, unknown> = {
    cardUrl: str(row.card_url),
    tier: str(row.tier),
    assay: str(row.assay),
    product: str(row.product),
    chemistryVersion: str(row.chemistry_version),
    sampleId: str(row.sample_id),
    studyId: str(row.study_id),
    canonicalKey: str(row.canonical_key),
    developmentStage: str(row.development_stage),
    captureArea: str(row.capture_area),
    defaultTable: str(row.default_table),
    nTranscripts:
      nTranscripts != null ? nTranscripts.toLocaleString("en-US") : null,
    storeSize: sizeBytes != null ? formatBytes(sizeBytes) : null,
    publicationYear: str(row.publication_year),
    licenseUrl: str(row.license_url),
    collections: collections.length ? collections : null,
  };

  for (const k of Object.keys(metadata)) {
    if (metadata[k] === null) delete metadata[k];
  }

  return {
    uid,
    data: {
      title: cap(str(row.title) ?? str(row.sample_id) ?? uid, 500),
      description,
      species: species ? cap(species, 200) : null,
      disease: disease ? cap(capitalize(disease), 200) : null,
      tissue: tissue ? cap(capitalize(tissue), 200) : null,
      platform: cap(platform, 200),
      tags: [
        "spatialdata-db",
        technology.toLowerCase(),
        organism?.toLowerCase(),
        tissue?.toLowerCase(),
      ]
        .filter((t): t is string => Boolean(t))
        .filter((t, i, a) => a.indexOf(t) === i),
      thumbnailUrl: str(row.thumbnail_url),
      metadata,
      externalLink: str(row.dataset_url) ?? str(row.card_url),
      publicationLink: doi ? `https://doi.org/${doi}` : null,
      numCells: nObs,
      numGenes: nFeatures,
      source: "sddb",
      sddbUid: uid,
      license,
      licenseNoncommercial: bool(row.license_noncommercial),
      sourceRelease: release,
      sourceGeneratedAt: generatedAt,
      isPublished: true,
      isBil: false,
      isCommunity: false,
      isInternal: false,
    },
    entry: {
      label: `${platform} store`,
      datasetType: "spatialdata",
      s3BaseUrl: storeUrlFromZarrUrl(zarrUrl),
    },
  };
}

/**
 * `n` rows spread over the catalog's technologies, round-robin in file
 * order, so a partial sync (a trial on the dev database) shows every kind of
 * store rather than the first `n` Visium rows.
 */
export function sampleRows<T extends { technology?: unknown }>(
  rows: T[],
  n: number,
): T[] {
  const groups = new Map<string, T[]>();

  for (const row of rows) {
    const key = String(row.technology ?? "");

    (groups.get(key) ?? groups.set(key, []).get(key)!).push(row);
  }
  const out: T[] = [];

  for (let i = 0; out.length < Math.min(n, rows.length); i++) {
    for (const group of groups.values()) {
      if (i < group.length && out.length < n) out.push(group[i]);
    }
  }

  return out;
}

/**
 * Which previously synced uids are no longer in the release: those rows are
 * unpublished (never deleted, so their views and any curation survive).
 */
export function staleUids(
  synced: Iterable<string>,
  inRelease: Iterable<string>,
): string[] {
  const keep = new Set(inRelease);

  return [...synced].filter((uid) => !keep.has(uid));
}
