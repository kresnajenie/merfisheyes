import type { Prisma } from "@prisma/client";

import { prisma } from "@/lib/prisma";
import {
  CARD_ORDER_BY,
  CARD_SELECT,
  findGeneMatchIds,
  findMatchedGenes,
  parseGeneTokens,
} from "@/lib/explore/query";
import {
  type ExploreCategory,
  catalogFacetWhere,
  catalogFeaturedWhere,
  catalogListWhere,
} from "@/lib/explore/visibility";

/**
 * The Explore search — one implementation for /api/explore and the agent API
 * (/api/agent/v1/catalog/search). Visibility is lib/explore/visibility.ts.
 */

export interface CatalogSearch {
  search?: string;
  species?: string;
  tissue?: string;
  platform?: string;
  /** Substring tokens, comma or space separated: every token must match. */
  genes?: string;
  /** Comma-separated symbols, each matched exactly (case-insensitive). */
  genesExact?: string;
  datasetType?: string;
  category?: ExploreCategory;
  skip: number;
  take: number;
}

const escapeIlike = (s: string) => s.replace(/[\\%_]/g, (m) => `\\${m}`);

export async function searchCatalog(q: CatalogSearch, isAdmin: boolean) {
  const conditions: Prisma.CatalogDatasetWhereInput[] = [
    catalogListWhere(q.category ?? "all", isAdmin),
  ];

  const search = q.search?.trim() ?? "";

  if (search) {
    // Tokenize on whitespace: each token must match SOME searched field
    // (title / description / bil_code / a tag element / metadata.investigator),
    // tokens AND'd together. Lets "mouse brain" match a row where one field
    // contains "mouse" and another contains "brain".
    const tokens = search.split(/\s+/).filter(Boolean).slice(0, 10);

    if (tokens.length > 0) {
      const args = tokens.map((t) => `%${escapeIlike(t)}%`);
      const perTokenSql = tokens
        .map((_, i) => {
          const p = `$${i + 1}`;

          return `(
            title ILIKE ${p} ESCAPE '\\' OR
            description ILIKE ${p} ESCAPE '\\' OR
            bil_code ILIKE ${p} ESCAPE '\\' OR
            metadata->>'investigator' ILIKE ${p} ESCAPE '\\' OR
            EXISTS (SELECT 1 FROM unnest(tags) tg WHERE tg ILIKE ${p} ESCAPE '\\')
          )`;
        })
        .join(" AND ");
      const matchingIds = await prisma.$queryRawUnsafe<{ id: string }[]>(
        `SELECT id FROM catalog_datasets WHERE ${perTokenSql}`,
        ...args,
      );

      conditions.push({ id: { in: matchingIds.map((r) => r.id) } });
    }
  }
  if (q.species)
    conditions.push({ species: { equals: q.species, mode: "insensitive" } });
  if (q.tissue)
    conditions.push({ tissue: { equals: q.tissue, mode: "insensitive" } });
  if (q.platform)
    conditions.push({ platform: { equals: q.platform, mode: "insensitive" } });
  if (q.datasetType)
    conditions.push({ entries: { some: { datasetType: q.datasetType } } });

  // Gene search: live as-you-type, case-insensitive substring per token.
  // "ntrk, bdnf" (or "ntrk bdnf") = datasets whose gene list matches BOTH
  // "ntrk" and "bdnf" somewhere. Each token's ILIKE hits the pg_trgm GIN
  // index on catalog_genes_text(genes).
  const geneTokens = parseGeneTokens(q.genes?.trim() ?? "");

  if (geneTokens.length > 0) {
    const ids = await findGeneMatchIds(geneTokens);

    conditions.push({ id: { in: ids } });
  }

  // Gene exact match: case-insensitive, each chip must match a gene exactly
  const geneList = (q.genesExact ?? "")
    .split(",")
    .map((g) => g.trim())
    .filter(Boolean);

  if (geneList.length > 0) {
    const matchingIds = await prisma.$queryRawUnsafe<{ id: string }[]>(
      `SELECT id FROM catalog_datasets WHERE ${geneList
        .map(
          (_, i) =>
            `EXISTS (SELECT 1 FROM unnest(genes) g WHERE g ILIKE $${i + 1})`,
        )
        .join(" AND ")}`,
      ...geneList,
    );

    conditions.push({ id: { in: matchingIds.map((r) => r.id) } });
  }

  const where: Prisma.CatalogDatasetWhereInput = { AND: conditions };

  const [items, total] = await Promise.all([
    prisma.catalogDataset.findMany({
      where,
      select: CARD_SELECT,
      orderBy: CARD_ORDER_BY,
      skip: q.skip,
      take: q.take,
    }),
    prisma.catalogDataset.count({ where }),
  ]);

  // Which genes actually matched, per card, for the current page only. The
  // exact chips are included so their hits show up too.
  const highlightTokens = [...geneTokens, ...geneList];
  const matchedByItem =
    highlightTokens.length > 0
      ? await findMatchedGenes(
          items.map((i) => i.id),
          highlightTokens,
        )
      : {};
  const itemsOut = items.map((i) =>
    highlightTokens.length > 0
      ? { ...i, matchedGenes: matchedByItem[i.id] ?? [] }
      : i,
  );

  return { items: itemsOut, total };
}

/** The featured carousel. */
export function featuredCatalog(isAdmin: boolean) {
  return prisma.catalogDataset.findMany({
    where: catalogFeaturedWhere(isAdmin),
    select: CARD_SELECT,
    orderBy: { sortOrder: "asc" },
  });
}

/** Distinct species / tissue / platform values over the rows "all" shows. */
export async function catalogFacets(isAdmin: boolean) {
  const where = catalogFacetWhere(isAdmin);
  const [speciesRaw, tissueRaw, platformRaw] = await Promise.all([
    prisma.catalogDataset.findMany({
      where: { AND: [where, { species: { not: null } }] },
      select: { species: true },
      distinct: ["species"],
      orderBy: { species: "asc" },
    }),
    prisma.catalogDataset.findMany({
      where: { AND: [where, { tissue: { not: null } }] },
      select: { tissue: true },
      distinct: ["tissue"],
      orderBy: { tissue: "asc" },
    }),
    prisma.catalogDataset.findMany({
      where: { AND: [where, { platform: { not: null } }] },
      select: { platform: true },
      distinct: ["platform"],
      orderBy: { platform: "asc" },
    }),
  ]);

  return {
    species: speciesRaw.map((r) => r.species).filter((v): v is string => !!v),
    tissues: tissueRaw.map((r) => r.tissue).filter((v): v is string => !!v),
    platforms: platformRaw
      .map((r) => r.platform)
      .filter((v): v is string => !!v),
  };
}
