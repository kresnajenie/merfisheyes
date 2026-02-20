import { NextRequest, NextResponse } from "next/server";
import { Prisma } from "@prisma/client";

import { prisma } from "@/lib/prisma";

const includeEntries = { entries: { orderBy: { sortOrder: "asc" as const } } };

// GET /api/explore — public endpoint for browsing published datasets
export async function GET(req: NextRequest) {
  const url = new URL(req.url);
  const search = url.searchParams.get("search")?.trim() ?? "";
  const species = url.searchParams.get("species") ?? "";
  const tissue = url.searchParams.get("tissue") ?? "";
  const platform = url.searchParams.get("platform") ?? "";
  const datasetType = url.searchParams.get("datasetType") ?? "";
  const page = Math.max(1, Number(url.searchParams.get("page") ?? "1"));
  const limit = Math.min(50, Math.max(1, Number(url.searchParams.get("limit") ?? "20")));
  const skip = (page - 1) * limit;

  // Build where clause
  const conditions: Prisma.CatalogDatasetWhereInput[] = [{ isPublished: true }];

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
      where: { isPublished: true, isFeatured: true },
      include: includeEntries,
      orderBy: { sortOrder: "asc" },
    }),
    prisma.catalogDataset.findMany({
      where: { isPublished: true, isBil: true },
      include: includeEntries,
      orderBy: { sortOrder: "asc" },
    }),
    // Get distinct filter values
    getDistinctFilters(),
  ]);

  return NextResponse.json({ items, total, page, limit, featured, bil, filters });
}

async function getDistinctFilters() {
  const published = { isPublished: true };

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
