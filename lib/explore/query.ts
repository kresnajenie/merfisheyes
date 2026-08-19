import type { Prisma } from "@prisma/client";

import { prisma } from "@/lib/prisma";

/**
 * Shared query pieces for the Explore catalog listing — used by both the SSR
 * seed (app/explore/page.tsx) and /api/explore so their payload shapes can't
 * drift.
 *
 * The core rule: list payloads NEVER include the `genes` array. A catalog row
 * carries hundreds of gene symbols (550+ for a typical panel); with 20 cards a
 * page plus the featured carousel, shipping them made every explore response
 * ~700 KB before this select existed. Gene matching happens in SQL (via the
 * pg_trgm index on catalog_genes_text) and only the per-card matches are
 * returned, capped, as `matchedGenes`.
 */

export const CARD_ENTRY_SELECT = {
  id: true,
  catalogId: true,
  label: true,
  datasetType: true,
  s3BaseUrl: true,
  datasetId: true,
  thumbnailUrl: true,
  sortOrder: true,
} satisfies Prisma.CatalogDatasetEntrySelect;

export const CARD_SELECT = {
  id: true,
  title: true,
  description: true,
  species: true,
  disease: true,
  institute: true,
  tissue: true,
  platform: true,
  tags: true,
  thumbnailUrl: true,
  bilCode: true,
  metadata: true,
  externalLink: true,
  publicationLink: true,
  isPublished: true,
  isFeatured: true,
  isBil: true,
  isInternal: true,
  isCommunity: true,
  sortOrder: true,
  numCells: true,
  numGenes: true,
  entries: {
    select: CARD_ENTRY_SELECT,
    orderBy: { sortOrder: "asc" as const },
  },
} satisfies Prisma.CatalogDatasetSelect;

/** Split a raw gene query into search tokens: commas and whitespace both
 * separate terms, so "ntrk, bdnf" and "ntrk bdnf" mean the same thing. */
export function parseGeneTokens(raw: string): string[] {
  return raw
    .split(/[,\s]+/)
    .map((t) => t.trim())
    .filter(Boolean)
    .slice(0, 10);
}

const escapeIlike = (s: string) => s.replace(/[\\%_]/g, (m) => `\\${m}`);

/**
 * Tables with a `genes text[]` column + a trgm index on
 * catalog_genes_text(genes). The table name is interpolated into SQL, so it
 * must come from this literal union — never from user input.
 */
export type GenesTable = "catalog_datasets" | "datasets";

/**
 * IDs of rows whose gene list matches EVERY token as a case-insensitive
 * substring. Each predicate hits the pg_trgm GIN index on
 * catalog_genes_text(genes) (migrations 20260803120000 / 20260819100000).
 */
export async function findGeneMatchIds(
  tokens: string[],
  table: GenesTable = "catalog_datasets",
): Promise<string[]> {
  if (tokens.length === 0) return [];
  const predicates = tokens
    .map((_, i) => `catalog_genes_text(genes) ILIKE $${i + 1} ESCAPE '\\'`)
    .join(" AND ");
  const rows = await prisma.$queryRawUnsafe<{ id: string }[]>(
    `SELECT id FROM ${table} WHERE ${predicates}`,
    ...tokens.map((t) => `%${escapeIlike(t)}%`),
  );

  return rows.map((r) => r.id);
}

/**
 * For the current page of items only, pull the actual gene symbols that
 * matched the tokens (any token, per gene), capped per dataset — this is what
 * the card renders as "matching genes" chips. Kept separate from the listing
 * query so the full genes array never leaves the database.
 */
export async function findMatchedGenes(
  ids: string[],
  tokens: string[],
  cap = 12,
  table: GenesTable = "catalog_datasets",
): Promise<Record<string, string[]>> {
  if (ids.length === 0 || tokens.length === 0) return {};

  // $1 = ids, $2..$(n+1) = tokens
  const tokenPreds = tokens
    .map((_, i) => `g ILIKE $${i + 2} ESCAPE '\\'`)
    .join(" OR ");
  const rows = await prisma.$queryRawUnsafe<
    { id: string; matched: string[] }[]
  >(
    `SELECT id,
            ARRAY(SELECT g FROM unnest(genes) g
                  WHERE ${tokenPreds}
                  ORDER BY g
                  LIMIT ${cap}) AS matched
       FROM ${table}
      WHERE id = ANY($1)`,
    ids,
    ...tokens.map((t) => `%${escapeIlike(t)}%`),
  );

  return Object.fromEntries(rows.map((r) => [r.id, r.matched]));
}
