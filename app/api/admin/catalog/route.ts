import { NextRequest, NextResponse } from "next/server";

import { prisma } from "@/lib/prisma";
import { requireAdmin } from "@/lib/admin-auth";

// GET /api/admin/catalog — list all catalog datasets
export async function GET(req: NextRequest) {
  const { error } = await requireAdmin();
  if (error) return error;

  const url = new URL(req.url);
  const search = url.searchParams.get("search") ?? "";
  const page = Math.max(1, Number(url.searchParams.get("page") ?? "1"));
  const limit = Math.min(100, Math.max(1, Number(url.searchParams.get("limit") ?? "50")));
  const skip = (page - 1) * limit;

  const where = search
    ? {
        OR: [
          { title: { contains: search, mode: "insensitive" as const } },
          { description: { contains: search, mode: "insensitive" as const } },
        ],
      }
    : {};

  const [items, total] = await Promise.all([
    prisma.catalogDataset.findMany({
      where,
      orderBy: [{ sortOrder: "asc" }, { createdAt: "desc" }],
      skip,
      take: limit,
    }),
    prisma.catalogDataset.count({ where }),
  ]);

  return NextResponse.json({ items, total, page, limit });
}

// POST /api/admin/catalog — create a new catalog dataset
export async function POST(req: NextRequest) {
  const { error, session } = await requireAdmin();
  if (error) return error;

  const body = await req.json();

  const item = await prisma.catalogDataset.create({
    data: {
      title: body.title,
      description: body.description ?? null,
      species: body.species ?? null,
      disease: body.disease ?? null,
      institute: body.institute ?? null,
      tissue: body.tissue ?? null,
      platform: body.platform ?? null,
      tags: body.tags ?? [],
      thumbnailUrl: body.thumbnailUrl ?? null,
      datasetType: body.datasetType,
      s3BaseUrl: body.s3BaseUrl ?? null,
      datasetId: body.datasetId ?? null,
      externalLink: body.externalLink ?? null,
      isPublished: body.isPublished ?? false,
      isFeatured: body.isFeatured ?? false,
      sortOrder: body.sortOrder ?? 0,
      numCells: body.numCells ?? null,
      numGenes: body.numGenes ?? null,
      createdBy: session!.user.id,
    },
  });

  return NextResponse.json(item, { status: 201 });
}
