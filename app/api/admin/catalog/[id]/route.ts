import { NextRequest, NextResponse } from "next/server";

import { prisma } from "@/lib/prisma";
import { requireAdmin } from "@/lib/admin-auth";

const includeEntries = { entries: { orderBy: { sortOrder: "asc" as const } } };

// GET /api/admin/catalog/[id]
export async function GET(
  _req: NextRequest,
  { params }: { params: Promise<{ id: string }> },
) {
  const { error } = await requireAdmin();
  if (error) return error;

  const { id } = await params;
  const item = await prisma.catalogDataset.findUnique({
    where: { id },
    include: includeEntries,
  });

  if (!item) return NextResponse.json({ error: "Not found" }, { status: 404 });

  return NextResponse.json(item);
}

// PATCH /api/admin/catalog/[id]
export async function PATCH(
  req: NextRequest,
  { params }: { params: Promise<{ id: string }> },
) {
  const { error } = await requireAdmin();
  if (error) return error;

  const { id } = await params;
  const body = await req.json();

  // Only allow updating known fields
  const allowedFields = [
    "title",
    "description",
    "species",
    "disease",
    "institute",
    "tissue",
    "platform",
    "tags",
    "thumbnailUrl",
    "externalLink",
    "isPublished",
    "isFeatured",
    "isBil",
    "sortOrder",
    "numCells",
    "numGenes",
  ] as const;

  const data: Record<string, unknown> = {};
  for (const key of allowedFields) {
    if (key in body) {
      data[key] = body[key];
    }
  }

  // Handle entries with delete-and-recreate
  if ("entries" in body && Array.isArray(body.entries)) {
    await prisma.catalogDatasetEntry.deleteMany({ where: { catalogId: id } });
    data.entries = {
      create: body.entries.map((e: { label: string; datasetType: string; s3BaseUrl?: string; datasetId?: string; sortOrder?: number }, i: number) => ({
        label: e.label,
        datasetType: e.datasetType,
        s3BaseUrl: e.s3BaseUrl || null,
        datasetId: e.datasetId || null,
        sortOrder: e.sortOrder ?? i,
      })),
    };
  }

  const item = await prisma.catalogDataset.update({
    where: { id },
    data,
    include: includeEntries,
  });

  return NextResponse.json(item);
}

// DELETE /api/admin/catalog/[id]
export async function DELETE(
  _req: NextRequest,
  { params }: { params: Promise<{ id: string }> },
) {
  const { error } = await requireAdmin();
  if (error) return error;

  const { id } = await params;

  await prisma.catalogDataset.delete({ where: { id } });

  return NextResponse.json({ success: true });
}
