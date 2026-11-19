import { readFile } from "node:fs/promises";

import {
  type AsyncBuffer,
  type FileMetaData,
  asyncBufferFromUrl,
  cachedAsyncBuffer,
  parquetMetadataAsync,
  parquetReadObjects,
} from "hyparquet";
import { compressors } from "hyparquet-compressors";

import { INDEX_URL } from "../config.ts";

/**
 * The published query index (MCP.md pillar B; `docs/specs/metadata.md` §7 on
 * `feature/metadata`): `catalog.parquet` — one row per live dataset, SDDB's
 * schema-3 columns then the `meyes_*` ones — with `ontology_rollup.json`,
 * `technology.json` and `release_info.json` beside it, and the gene index.
 *
 * One reader for both servers (D23): the package on a person's computer and
 * the site's `/mcp` load the same files the same way. Only the columns a
 * search reads are decoded (all 131 cost ~230 MB of objects); the files are
 * fetched again only when `release_info.json` names a new release.
 */

/** The catalog columns search reads, filters on or returns. */
export const CATALOG_COLUMNS = [
  "uid",
  "canonical_key",
  "study_id",
  "collection_name",
  "organism",
  "organism_id",
  "tissue",
  "tissue_id",
  "tissue_system",
  "tissue_system_id",
  "tissue_organ",
  "tissue_organ_id",
  "disease",
  "disease_id",
  "disease_class",
  "disease_group",
  "development_stage",
  "technology",
  "vendor",
  "assay",
  "panel_size",
  "validation_status",
  "tier",
  "license",
  "n_obs",
  "n_features",
  "n_transcripts",
  "published_at",
  "publication_date",
  "title",
  "sample_title",
  "sample_id",
  "donor_id",
  "doi",
  "pmid",
  "modality",
  "obs_unit",
  "timepoint",
  "timepoint_kind",
  "timepoint_days",
  "treatment",
  "genotype",
  "condition",
  "same_section_group",
  "same_section_relation",
  "same_section_release",
  "meyes_id",
  "meyes_source",
  "meyes_kind",
  "meyes_page_url",
  "meyes_viewer_url",
  "meyes_description",
  "meyes_sex",
  "meyes_n_sections",
  "meyes_dims",
  "meyes_format",
  "meyes_metadata_from",
  "meyes_bil_code",
  "meyes_investigator",
  "meyes_institution",
  "meyes_cell_type_levels",
] as const;

type Column = (typeof CATALOG_COLUMNS)[number];

/** One catalog row: strings, numbers (Int64 converted), string lists, or
 * null. */
export type Row = { [K in Column]: string | number | string[] | null } & {
  /** Position in `catalog.parquet` — the bit a gene bitset sets. */
  row: number;
};

export interface RollupTerm {
  ancestors_present: string[];
  name: string;
}

export interface Technology {
  id: string;
  label: string;
  vendor_id?: string;
  legacy_assay?: string[];
}

export interface GeneStat {
  symbol: string;
  feature_id: string;
  total_counts: number | null;
  fraction_obs_detected: number | null;
  meyes_id: string;
}

export interface LibraryIndex {
  /** `release_info.json` `release`, e.g. "meyes-index-2026-11-18". */
  release: string;
  builtAt: string | null;
  source: string;
  rows: Row[];
  byMeyesId: Map<string, Row>;
  rollup: Record<string, RollupTerm>;
  technologies: Technology[];
  /** Loaded on the first expression filter. */
  geneIndex(): Promise<GeneIndex>;
  /** The datasets that measure a gene, over every live row: one range read
   * of `gene_lookup.parquet`. */
  measured(symbol: string): Promise<Measured>;
}

export interface Measured {
  /** The gene entities whose symbols include it (one per organism and
   * annotation: ENSMUSG…, ENSG…, a bare symbol). */
  entities: { entity_id: string; feature_id: string; n_datasets: number }[];
  /** Catalog row positions (`Row.row`). */
  rows: Set<number>;
}

export interface GeneIndex {
  /** Lower-case symbol or feature id → the stats of every dataset of ours
   * that measures it. */
  bySymbol: Map<string, GeneStat[]>;
  /** The meyes_ids the gene index covers (stores and Explore-only rows). */
  covered: Set<string>;
}

/** A file of the index: over HTTP, or from a folder (tests, a local
 * build). */
async function fetchFile(
  base: string,
  name: string,
  signal?: AbortSignal,
): Promise<ArrayBuffer> {
  if (!/^https?:\/\//.test(base)) {
    const b = await readFile(`${base.replace(/^file:\/\//, "")}/${name}`);

    return b.buffer.slice(b.byteOffset, b.byteOffset + b.byteLength);
  }
  const res = await fetch(`${base}/${name}`, { signal, cache: "no-store" });

  if (!res.ok)
    throw new Error(
      `The library index could not be read (${name}: HTTP ${res.status}). Try again later.`,
    );

  return res.arrayBuffer();
}

async function fetchJson<T>(
  base: string,
  name: string,
  signal?: AbortSignal,
): Promise<T> {
  return JSON.parse(
    new TextDecoder().decode(await fetchFile(base, name, signal)),
  ) as T;
}

/** BigInt (parquet Int64) → number; numpy object lists stay arrays. */
function plain(v: unknown): string | number | string[] | null {
  if (v === undefined || v === null) return null;
  if (typeof v === "bigint") return Number(v);
  if (Array.isArray(v)) return v.map(String);
  if (typeof v === "number" && Number.isNaN(v)) return null;

  return v as string | number;
}

async function readParquet(
  buf: ArrayBuffer,
  columns: readonly string[],
): Promise<Record<string, unknown>[]> {
  return (await parquetReadObjects({
    file: buf,
    columns: [...columns],
    compressors,
  })) as Record<string, unknown>[];
}

async function loadIndex(
  base: string,
  release: { release: string; built_at?: string },
  signal?: AbortSignal,
): Promise<LibraryIndex> {
  const [catalog, rollup, tech] = await Promise.all([
    fetchFile(base, "catalog.parquet", signal),
    fetchJson<Record<string, RollupTerm>>(base, "ontology_rollup.json", signal),
    fetchJson<{ technologies: Technology[] }>(base, "technology.json", signal),
  ]);
  const raw = await readParquet(catalog, CATALOG_COLUMNS);
  const rows = raw.map((r, i) => {
    const out = { row: i } as Row;

    for (const c of CATALOG_COLUMNS) out[c] = plain(r[c]);

    return out;
  });
  let genes: Promise<GeneIndex> | null = null;
  let lookup: Promise<{ file: AsyncBuffer; metadata: FileMetaData }> | null =
    null;
  const measuredCache = new Map<string, Promise<Measured>>();

  const openLookup = async () => {
    // The bitsets are over catalog.parquet's row order: refuse them if the
    // files on either side were rebuilt apart.
    const meta = await fetchJson<{ meyes_ids: string[] }>(
      base,
      "gene_search_meta.json",
    );

    if (
      meta.meyes_ids.length !== rows.length ||
      meta.meyes_ids.some((id, i) => id !== rows[i].meyes_id)
    )
      throw new Error(
        "The gene index and the catalog are from different builds; try again in a few minutes.",
      );
    const file = /^https?:\/\//.test(base)
      ? cachedAsyncBuffer(
          await asyncBufferFromUrl({ url: `${base}/gene_lookup.parquet` }),
        )
      : // A folder (tests, a local build): read whole.
        await fetchFile(base, "gene_lookup.parquet");

    return { file, metadata: await parquetMetadataAsync(file) };
  };

  return {
    release: release.release,
    builtAt: release.built_at ?? null,
    source: base,
    rows,
    byMeyesId: new Map(rows.map((r) => [String(r.meyes_id), r])),
    rollup,
    technologies: tech.technologies ?? [],
    geneIndex() {
      genes ??= (async () => {
        const buf = await fetchFile(base, "gene_index.parquet");
        const recs = await readParquet(buf, [
          "symbol",
          "feature_id",
          "total_counts",
          "fraction_obs_detected",
          "meyes_id",
        ]);
        const bySymbol = new Map<string, GeneStat[]>();
        const covered = new Set<string>();

        for (const r of recs) {
          const s: GeneStat = {
            symbol: String(r.symbol),
            feature_id: String(r.feature_id),
            total_counts: plain(r.total_counts) as number | null,
            fraction_obs_detected: plain(r.fraction_obs_detected) as
              | number
              | null,
            meyes_id: String(r.meyes_id),
          };

          covered.add(s.meyes_id);
          for (const k of new Set([
            s.symbol.toLowerCase(),
            s.feature_id.toLowerCase(),
          ])) {
            const list = bySymbol.get(k);

            if (list) list.push(s);
            else bySymbol.set(k, [s]);
          }
        }

        return { bySymbol, covered };
      })();
      genes.catch(() => (genes = null));

      return genes;
    },
    measured(symbol: string) {
      const q = symbol.toLowerCase().trim();
      let p = measuredCache.get(q);

      if (p) return p;
      p = (async () => {
        lookup ??= openLookup();
        lookup.catch(() => (lookup = null));
        const { file, metadata } = await lookup;
        const out: Measured = { entities: [], rows: new Set() };
        let start = 0;

        // Sorted by symbol: the groups whose min/max hold it (usually one).
        for (const g of metadata.row_groups) {
          const n = Number(g.num_rows);
          const st = g.columns[0].meta_data?.statistics;
          const lo = st?.min_value as string | undefined;
          const hi = st?.max_value as string | undefined;

          if (lo === undefined || hi === undefined || (lo <= q && q <= hi)) {
            const recs = (await parquetReadObjects({
              file,
              metadata,
              compressors,
              utf8: false,
              rowStart: start,
              rowEnd: start + n,
            })) as {
              symbol: string;
              entity_id: string;
              feature_id: string;
              n_datasets: number;
              bitset: Uint8Array;
            }[];

            for (const r of recs) {
              if (r.symbol !== q) continue;
              out.entities.push({
                entity_id: r.entity_id,
                feature_id: r.feature_id,
                n_datasets: Number(r.n_datasets),
              });
              r.bitset.forEach((byte, i) => {
                for (let b = 0; byte && b < 8; b++, byte >>= 1)
                  if (byte & 1) out.rows.add(i * 8 + b);
              });
            }
          }
          start += n;
        }

        return out;
      })();
      p.catch(() => measuredCache.delete(q));
      if (measuredCache.size > 500) measuredCache.clear();
      measuredCache.set(q, p);

      return p;
    },
  };
}

/** How often a long-running server asks whether a new release is out. */
const RECHECK_MS = 10 * 60 * 1000;

let cached: { base: string; index: LibraryIndex; checkedAt: number } | null =
  null;
let loading: Promise<LibraryIndex> | null = null;

/**
 * The index, loaded once per release and shared by every tool call. A
 * failed recheck keeps the copy in hand; a failed first load throws.
 */
export function libraryIndex(
  base: string = INDEX_URL,
  signal?: AbortSignal,
): Promise<LibraryIndex> {
  const now = Date.now();

  if (cached && cached.base === base && now - cached.checkedAt < RECHECK_MS)
    return Promise.resolve(cached.index);
  loading ??= (async () => {
    try {
      const info = await fetchJson<{ release: string; built_at?: string }>(
        base,
        "release_info.json",
        signal,
      ).catch((e) => {
        if (cached?.base === base) return null;
        throw e;
      });

      if (!info) {
        cached!.checkedAt = now;

        return cached!.index;
      }
      if (cached?.base === base && cached.index.release === info.release) {
        cached.checkedAt = now;

        return cached.index;
      }
      const index = await loadIndex(base, info, signal);

      cached = { base, index, checkedAt: now };

      return index;
    } finally {
      loading = null;
    }
  })();

  return loading;
}

/** For tests: forget the loaded index. */
export function resetLibraryIndex() {
  cached = null;
  loading = null;
}
