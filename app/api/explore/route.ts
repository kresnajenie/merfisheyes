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
  const datasetType = url.searchParams.get("datasetType") ?? "";
  const tab = url.searchParams.get("tab") ?? "";
  const page = Math.max(1, Number(url.searchParams.get("page") ?? "1"));
  const limit = Math.min(100, Math.max(1, Number(url.searchParams.get("limit") ?? "50")));
  const skip = (page - 1) * limit;

  // Build where clause
  const conditions: Prisma.CatalogDatasetWhereInput[] = [{ isPublished: true }];

  // Tab-specific filters
  if (tab === "internal" && isAdmin) {
    conditions.push({ isInternal: true });
  } else {
    // Non-internal tabs: always exclude internal datasets
    conditions.push({ isInternal: false });
  }
  if (tab === "featured") conditions.push({ isFeatured: true });
  if (tab === "bil") conditions.push({ isBil: true });

  if (search) {
    conditions.push({
      OR: [
        { title: { contains: search, mode: "insensitive" } },
        { description: { contains: search, mode: "insensitive" } },
        { tags: { hasSome: [search] } },
      ],
    });
  }
  if (species) conditions.push({ species: { equals: species, mode: "insensitive" } });
  if (tissue) conditions.push({ tissue: { equals: tissue, mode: "insensitive" } });
  if (platform) conditions.push({ platform: { equals: platform, mode: "insensitive" } });
  if (datasetType) conditions.push({ entries: { some: { datasetType } } });

  const where: Prisma.CatalogDatasetWhereInput = { AND: conditions };

  // Base filter for featured/bil: exclude internal for non-admins
  const publicBase: Prisma.CatalogDatasetWhereInput = isAdmin
    ? { isPublished: true }
    : { isPublished: true, isInternal: false };

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
      where: { ...publicBase, isFeatured: true },
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
  const published: Prisma.CatalogDatasetWhereInput = isAdmin
    ? { isPublished: true }
    : { isPublished: true, isInternal: false };

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
