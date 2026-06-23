import { NextResponse } from "next/server";

import { prisma } from "@/lib/prisma";
import { requireAdmin } from "@/lib/admin-auth";

/**
 * GET /api/admin/catalog/export
 *
 * Streams every catalog dataset (including unpublished + internal) as a
 * downloadable JSON array. Round-trippable: the response shape matches the
 * `items` payload accepted by `/api/admin/catalog/import`. Includes `id`
 * and `bilCode` so the importer can match-and-update existing rows.
 */
export async function GET() {
  const { error } = await requireAdmin();

  if (error) return error;

  const rows = await prisma.catalogDataset.findMany({
    include: { entries: { orderBy: { sortOrder: "asc" } } },
    orderBy: [{ sortOrder: "asc" }, { createdAt: "desc" }],
  });

  // Strip audit / FK columns; keep everything an admin would want when
  // editing the JSON by hand and re-uploading.
  const items = rows.map((r) => ({
    id: r.id,
    title: r.title,
    description: r.description,
    species: r.species,
    disease: r.disease,
    institute: r.institute,
    tissue: r.tissue,
    platform: r.platform,
    tags: r.tags,
    genes: r.genes,
    thumbnailUrl: r.thumbnailUrl,
    bilCode: r.bilCode,
    metadata: r.metadata,
    externalLink: r.externalLink,
    publicationLink: r.publicationLink,
    isPublished: r.isPublished,
    isFeatured: r.isFeatured,
    isBil: r.isBil,
    isInternal: r.isInternal,
    sortOrder: r.sortOrder,
    numCells: r.numCells,
    numGenes: r.numGenes,
    entries: r.entries.map((e) => ({
      label: e.label,
      datasetType: e.datasetType,
      s3BaseUrl: e.s3BaseUrl,
      datasetId: e.datasetId,
      thumbnailUrl: e.thumbnailUrl,
      sortOrder: e.sortOrder,
    })),
  }));

  const stamp = new Date().toISOString().slice(0, 10);

  return new NextResponse(JSON.stringify(items, null, 2), {
    status: 200,
    headers: {
      "Content-Type": "application/json",
      "Content-Disposition": `attachment; filename="catalog-${stamp}.json"`,
    },
  });
}
