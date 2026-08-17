import { NextRequest, NextResponse } from "next/server";
import { Prisma } from "@prisma/client";

import { prisma } from "@/lib/prisma";
import { auth } from "@/lib/auth";

const includeEntries = { entries: { orderBy: { sortOrder: "asc" as const } } };

// GET /api/explore — public endpoint for browsing published datasets
export async function GET(req: NextRequest) {
  const session = await auth();
  const isAdmin = session?.user?.role === "ADMIN" || session?.user?.role === "SUPER_ADMIN";

  const url = new URL(req.url);
  const search = url.searchParams.get("search")?.trim() ?? "";
  const species = url.searchParams.get("species") ?? "";
  const tissue = url.searchParams.get("tissue") ?? "";
  const platform = url.searchParams.get("platform") ?? "";
  const genesParam = url.searchParams.get("genes")?.trim() ?? "";
  const genesExactParam = url.searchParams.get("genesExact")?.trim() ?? "";
  const datasetType = url.searchParams.get("datasetType") ?? "";
  const tab = url.searchParams.get("tab") ?? "";
  const page = Math.max(1, Number(url.searchParams.get("page") ?? "1"));
  const limit = Math.min(100, Math.max(1, Number(url.searchParams.get("limit") ?? "50")));
  const skip = (page - 1) * limit;

  // Category filter (was tabs): "all" | "featured" | "bil" | "community" |
  // "internal". "all" now includes approved community submissions alongside
  // curated content.
  const category = tab || "all";
  const includeCommunity = category === "all" || category === "community";

  // App-uploaded (community) datasets are viewable by datasetId (no public
  // s3BaseUrl); curated content requires s3BaseUrl. When community is in scope,
  // accept either.
  const hasViewableEntry: Prisma.CatalogDatasetWhereInput = includeCommunity
    ? {
        entries: {
          some: {
            OR: [{ s3BaseUrl: { not: null } }, { datasetId: { not: null } }],
          },
        },
      }
    : { entries: { some: { s3BaseUrl: { not: null } } } };

  // Build where clause — only datasets with at least one viewable entry
  const conditions: Prisma.CatalogDatasetWhereInput[] = [
    { isPublished: true },
    hasViewableEntry,
  ];

  if (category === "internal" && isAdmin) {
    conditions.push({ isInternal: true }, { isCommunity: false });
  } else {
    // Every non-internal category excludes internal datasets.
    conditions.push({ isInternal: false });
  }

  if (category === "community") {
    // User submissions, only once an admin has approved them.
    conditions.push({ isCommunity: true }, { reviewStatus: "APPROVED" });
  } else if (category === "all") {
    // Curated content plus approved community submissions.
    conditions.push({
      OR: [
        { isCommunity: false },
        { AND: [{ isCommunity: true }, { reviewStatus: "APPROVED" }] },
      ],
    });
  } else if (category !== "internal") {
    // featured / bil — curated only (internal already excluded community).
    conditions.push({ isCommunity: false });
  }

  if (category === "featured") conditions.push({ isFeatured: true });
  if (category === "bil") conditions.push({ isBil: true });

  if (search) {
    // Tokenize on whitespace: each token must match SOME searched field
    // (title / description / bil_code / a tag element / metadata.investigator),
    // tokens AND'd together. Lets "mouse brain" match a row where one field
    // contains "mouse" and another contains "brain".
    const tokens = search.split(/\s+/).filter(Boolean).slice(0, 10);
    if (tokens.length > 0) {
      const escapeIlike = (s: string) => s.replace(/[\\%_]/g, (m) => `\\${m}`);
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
  if (species) conditions.push({ species: { equals: species, mode: "insensitive" } });
  if (tissue) conditions.push({ tissue: { equals: tissue, mode: "insensitive" } });
  if (platform) conditions.push({ platform: { equals: platform, mode: "insensitive" } });
  if (datasetType) conditions.push({ entries: { some: { datasetType } } });

  // Gene search: case-insensitive substring match (ILIKE %term%) — live as-you-type.
  // Matches against the genes joined into one string so the pg_trgm GIN index on
  // catalog_genes_text(genes) is used (equivalent to "some gene contains term"
  // for single-token gene symbols, which never contain spaces).
  if (genesParam) {
    const term = genesParam.trim().replace(/\s+/g, "");
    if (term) {
      const matchingIds = await prisma.$queryRawUnsafe<{ id: string }[]>(
        `SELECT id FROM catalog_datasets WHERE catalog_genes_text(genes) ILIKE $1`,
        `%${term}%`,
      );
      const ids = matchingIds.map((r) => r.id);
      conditions.push({ id: { in: ids } });
    }
  }

  // Gene exact match: case-insensitive, each chip must match a gene exactly
  if (genesExactParam) {
    const geneList = genesExactParam.split(",").map((g) => g.trim()).filter(Boolean);
    if (geneList.length > 0) {
      const matchingIds = await prisma.$queryRawUnsafe<{ id: string }[]>(
        `SELECT id FROM catalog_datasets WHERE ${geneList
          .map((_, i) => `EXISTS (SELECT 1 FROM unnest(genes) g WHERE g ILIKE $${i + 1})`)
          .join(" AND ")}`,
        ...geneList,
      );
      const ids = matchingIds.map((r) => r.id);
      conditions.push({ id: { in: ids } });
    }
  }

  const where: Prisma.CatalogDatasetWhereInput = { AND: conditions };

  // Base filter for bil: curated only, exclude internal for non-admins.
  const hasEntry = { entries: { some: { s3BaseUrl: { not: null } } } };
  const publicBase: Prisma.CatalogDatasetWhereInput = isAdmin
    ? { isPublished: true, isCommunity: false, ...hasEntry }
    : { isPublished: true, isCommunity: false, isInternal: false, ...hasEntry };

  // Featured carousel: curated PLUS approved community submissions (a community
  // dataset can be flagged featured). Viewable by s3BaseUrl or datasetId.
  const hasEntryAny = {
    entries: {
      some: { OR: [{ s3BaseUrl: { not: null } }, { datasetId: { not: null } }] },
    },
  };
  const communityOr: Prisma.CatalogDatasetWhereInput = {
    OR: [
      { isCommunity: false },
      { AND: [{ isCommunity: true }, { reviewStatus: "APPROVED" }] },
    ],
  };
  const featuredBase: Prisma.CatalogDatasetWhereInput = isAdmin
    ? { isPublished: true, ...communityOr, ...hasEntryAny }
    : { isPublished: true, isInternal: false, ...communityOr, ...hasEntryAny };

  // Fetch results, featured separately
  const [items, total, featured, bil, filters] = await Promise.all([
    prisma.catalogDataset.findMany({
      where,
      include: includeEntries,
      orderBy: [{ sortOrder: "asc" }, { createdAt: "desc" }],
      skip,
      take: limit,
    }),
    prisma.catalogDataset.count({ where }),
    prisma.catalogDataset.findMany({
      where: { ...featuredBase, isFeatured: true },
      include: includeEntries,
      orderBy: { sortOrder: "asc" },
    }),
    prisma.catalogDataset.findMany({
      where: { ...publicBase, isBil: true },
      include: includeEntries,
      orderBy: { sortOrder: "asc" },
    }),
    // Get distinct filter values
    getDistinctFilters(isAdmin),
  ]);

  return NextResponse.json({ items, total, page, limit, featured, bil, filters });
}

async function getDistinctFilters(isAdmin: boolean) {
  // Filter options reflect everything shown in "All": curated content plus
  // approved community submissions (internal excluded for non-admins).
  const communityOr: Prisma.CatalogDatasetWhereInput = {
    OR: [
      { isCommunity: false },
      { AND: [{ isCommunity: true }, { reviewStatus: "APPROVED" }] },
    ],
  };
  const published: Prisma.CatalogDatasetWhereInput = isAdmin
    ? { isPublished: true, ...communityOr }
    : { isPublished: true, isInternal: false, ...communityOr };

  const [speciesRaw, tissueRaw, platformRaw] = await Promise.all([
    prisma.catalogDataset.findMany({
      where: { ...published, species: { not: null } },
      select: { species: true },
      distinct: ["species"],
      orderBy: { species: "asc" },
    }),
    prisma.catalogDataset.findMany({
      where: { ...published, tissue: { not: null } },
      select: { tissue: true },
      distinct: ["tissue"],
      orderBy: { tissue: "asc" },
    }),
    prisma.catalogDataset.findMany({
      where: { ...published, platform: { not: null } },
      select: { platform: true },
      distinct: ["platform"],
      orderBy: { platform: "asc" },
    }),
  ]);

  return {
    species: speciesRaw.map((r) => r.species).filter(Boolean) as string[],
    tissues: tissueRaw.map((r) => r.tissue).filter(Boolean) as string[],
    platforms: platformRaw.map((r) => r.platform).filter(Boolean) as string[],
  };
}
