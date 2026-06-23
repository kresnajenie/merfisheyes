import { NextRequest, NextResponse } from "next/server";
import { Prisma } from "@prisma/client";

import { prisma } from "@/lib/prisma";
import { auth } from "@/lib/auth";
import { requireAdmin } from "@/lib/admin-auth";

interface EntryPayload {
  label?: string;
  datasetType?: string;
  s3BaseUrl?: string | null;
  datasetId?: string | null;
  thumbnailUrl?: string | null;
  sortOrder?: number;
}

interface ItemPayload {
  id?: string;
  title?: string;
  description?: string | null;
  species?: string | null;
  disease?: string | null;
  institute?: string | null;
  tissue?: string | null;
  platform?: string | null;
  tags?: string[];
  genes?: string[];
  thumbnailUrl?: string | null;
  bilCode?: string | null;
  metadata?: Record<string, unknown> | null;
  externalLink?: string | null;
  publicationLink?: string | null;
  isPublished?: boolean;
  isFeatured?: boolean;
  isBil?: boolean;
  isInternal?: boolean;
  sortOrder?: number;
  numCells?: number | null;
  numGenes?: number | null;
  entries?: EntryPayload[];
}

interface ImportBody {
  mode?: "add" | "overwrite" | "replace_all";
  items?: ItemPayload | ItemPayload[];
}

/**
 * If `link` looks like a merfisheyes viewer URL, extract dataset type +
 * either a datasetId or an inner s3BaseUrl. Returns null when the link
 * isn't a known viewer path (e.g. a raw S3 URL or unrelated host).
 *
 * Supported shapes (host-agnostic — only the path matters):
 *   /viewer/{id}
 *   /sm-viewer/{id}
 *   /viewer/from-s3?url={encodedS3Url}
 *   /sm-viewer/from-s3?url={encodedS3Url}
 */
function parseViewerLink(link: string): {
  datasetType: string;
  datasetId: string | null;
  s3BaseUrl: string | null;
} | null {
  try {
    const url = new URL(link);
    const parts = url.pathname.split("/").filter(Boolean);

    if (parts.length < 1) return null;

    const isSm = parts[0] === "sm-viewer";
    const isSc = parts[0] === "viewer";

    if (!isSm && !isSc) return null;

    const datasetType = isSm ? "single_molecule" : "single_cell";

    if (parts[1] === "from-s3") {
      const inner = url.searchParams.get("url");

      if (!inner) return null;

      return {
        datasetType,
        datasetId: null,
        s3BaseUrl: decodeURIComponent(inner),
      };
    }
    if (parts.length >= 2 && parts[1] !== "from-local") {
      return { datasetType, datasetId: parts[1], s3BaseUrl: null };
    }

    return null;
  } catch {
    return null;
  }
}

/**
 * Resolve an entry payload into its final create-shape: trims labels,
 * auto-derives type/id/s3BaseUrl from a viewer link pasted into the
 * s3BaseUrl slot, and applies a sensible default label.
 */
function normalizeEntry(e: EntryPayload, i: number) {
  let datasetType = e.datasetType ?? "single_cell";
  let datasetId = e.datasetId || null;
  let s3BaseUrl = e.s3BaseUrl || null;

  // If s3BaseUrl is actually a merfisheyes link, parse it and let it drive
  // the three fields. Raw S3 URLs (anything parseViewerLink can't read)
  // pass through unchanged, in which case the caller's `datasetType` wins.
  if (s3BaseUrl) {
    const parsed = parseViewerLink(s3BaseUrl);

    if (parsed) {
      datasetType = parsed.datasetType;
      datasetId = parsed.datasetId;
      s3BaseUrl = parsed.s3BaseUrl;
    }
  }

  return {
    label:
      e.label?.trim() ||
      (datasetType === "single_molecule" ? "Single Molecule" : "Single Cell"),
    datasetType,
    s3BaseUrl,
    datasetId,
    thumbnailUrl: e.thumbnailUrl || null,
    sortOrder: e.sortOrder ?? i,
  };
}

const buildCreateData = (item: ItemPayload, createdBy: string) => ({
  title: item.title!,
  description: item.description ?? null,
  species: item.species ?? null,
  disease: item.disease ?? null,
  institute: item.institute ?? null,
  tissue: item.tissue ?? null,
  platform: item.platform ?? null,
  tags: item.tags ?? [],
  genes: item.genes ?? [],
  thumbnailUrl: item.thumbnailUrl ?? null,
  bilCode: item.bilCode ?? null,
  metadata: (item.metadata as Prisma.InputJsonValue) ?? {},
  externalLink: item.externalLink ?? null,
  publicationLink: item.publicationLink ?? null,
  isPublished: item.isPublished ?? false,
  isFeatured: item.isFeatured ?? false,
  isBil: item.isBil ?? false,
  isInternal: item.isInternal ?? false,
  sortOrder: item.sortOrder ?? 0,
  numCells: item.numCells ?? null,
  numGenes: item.numGenes ?? null,
  createdBy,
  entries: {
    create: (item.entries ?? []).map(normalizeEntry),
  },
});

const buildUpdateData = (item: ItemPayload) => {
  const data: Prisma.CatalogDatasetUpdateInput = {};

  if ("title" in item) data.title = item.title!;
  if ("description" in item) data.description = item.description ?? null;
  if ("species" in item) data.species = item.species ?? null;
  if ("disease" in item) data.disease = item.disease ?? null;
  if ("institute" in item) data.institute = item.institute ?? null;
  if ("tissue" in item) data.tissue = item.tissue ?? null;
  if ("platform" in item) data.platform = item.platform ?? null;
  if ("tags" in item) data.tags = item.tags ?? [];
  if ("genes" in item) data.genes = item.genes ?? [];
  if ("thumbnailUrl" in item) data.thumbnailUrl = item.thumbnailUrl ?? null;
  if ("bilCode" in item) data.bilCode = item.bilCode ?? null;
  if ("metadata" in item)
    data.metadata = (item.metadata as Prisma.InputJsonValue) ?? {};
  if ("externalLink" in item) data.externalLink = item.externalLink ?? null;
  if ("publicationLink" in item)
    data.publicationLink = item.publicationLink ?? null;
  if ("isPublished" in item) data.isPublished = item.isPublished ?? false;
  if ("isFeatured" in item) data.isFeatured = item.isFeatured ?? false;
  if ("isBil" in item) data.isBil = item.isBil ?? false;
  if ("isInternal" in item) data.isInternal = item.isInternal ?? false;
  if ("sortOrder" in item) data.sortOrder = item.sortOrder ?? 0;
  if ("numCells" in item) data.numCells = item.numCells ?? null;
  if ("numGenes" in item) data.numGenes = item.numGenes ?? null;

  return data;
};

/**
 * POST /api/admin/catalog/import
 *
 * Body: { mode: "add" | "overwrite", items: ItemPayload | ItemPayload[] }
 *
 *   add        — match by id (then bilCode); matched rows are skipped, the
 *                rest are created
 *   overwrite  — same matching; matched rows are updated in place (entries
 *                replaced via delete-and-recreate), unmatched rows are created
 *
 * Returns { mode, created, updated, skipped, errors } where each list
 * holds { index, id?, title?, reason? } records the caller can surface.
 */
export async function POST(req: NextRequest) {
  const { error, session } = await requireAdmin();

  if (error) return error;

  let body: ImportBody;

  try {
    body = await req.json();
  } catch {
    return NextResponse.json({ error: "Invalid JSON" }, { status: 400 });
  }

  const mode: "add" | "overwrite" | "replace_all" =
    body.mode === "overwrite"
      ? "overwrite"
      : body.mode === "replace_all"
        ? "replace_all"
        : "add";

  // replace_all is destructive at catalog scope (can delete rows) — gate
  // it to SUPER_ADMIN regardless of who can reach the rest of the endpoint.
  if (mode === "replace_all") {
    const session = await auth();

    if (session?.user.role !== "SUPER_ADMIN") {
      return NextResponse.json(
        { error: "replace_all requires SUPER_ADMIN" },
        { status: 403 },
      );
    }
  }

  const rawItems = body.items;

  if (!rawItems) {
    return NextResponse.json({ error: "`items` is required" }, { status: 400 });
  }
  const items = Array.isArray(rawItems) ? rawItems : [rawItems];

  const created: { index: number; id: string; title: string }[] = [];
  const updated: { index: number; id: string; title: string }[] = [];
  const skipped: { index: number; title: string; reason: string }[] = [];
  const errors: { index: number; title?: string; error: string }[] = [];
  const deleted: { id: string; title: string }[] = [];
  // Track every row touched by the import so replace_all can delete the rest.
  const touchedIds = new Set<string>();

  for (let i = 0; i < items.length; i++) {
    const item = items[i];

    if (!item || typeof item !== "object") {
      errors.push({ index: i, error: "item is not an object" });
      continue;
    }
    if (!item.title || typeof item.title !== "string" || !item.title.trim()) {
      errors.push({ index: i, title: item.title, error: "title is required" });
      continue;
    }

    // Locate an existing row: id first, then bilCode.
    let existing: { id: string } | null = null;

    if (item.id) {
      existing = await prisma.catalogDataset.findUnique({
        where: { id: item.id },
        select: { id: true },
      });
    }
    if (!existing && item.bilCode) {
      existing = await prisma.catalogDataset.findUnique({
        where: { bilCode: item.bilCode },
        select: { id: true },
      });
    }

    try {
      if (existing) {
        if (mode === "add") {
          skipped.push({
            index: i,
            title: item.title,
            reason: item.id
              ? `id ${item.id} already exists`
              : `bilCode ${item.bilCode} already exists`,
          });
          continue;
        }
        // overwrite / replace_all: entries are delete-and-recreated in a
        // transaction so a failure doesn't leave the row half-updated.
        await prisma.$transaction(async (tx) => {
          await tx.catalogDatasetEntry.deleteMany({
            where: { catalogId: existing!.id },
          });

          const updateData = buildUpdateData(item);

          if ("entries" in item) {
            updateData.entries = {
              create: (item.entries ?? []).map(normalizeEntry),
            };
          }

          await tx.catalogDataset.update({
            where: { id: existing!.id },
            data: updateData,
          });
        });

        updated.push({ index: i, id: existing.id, title: item.title });
        touchedIds.add(existing.id);
      } else {
        const row = await prisma.catalogDataset.create({
          data: buildCreateData(item, session!.user.id),
        });

        created.push({ index: i, id: row.id, title: item.title });
        touchedIds.add(row.id);
      }
    } catch (e) {
      errors.push({
        index: i,
        title: item.title,
        error: (e as Error).message,
      });
    }
  }

  // replace_all: anything in the DB that the JSON didn't touch gets removed
  // so the catalog ends up as a faithful mirror of the upload. Entries
  // cascade-delete via the schema. Only runs when there were no upstream
  // errors — leaves the DB intact when something already failed.
  if (mode === "replace_all" && errors.length === 0) {
    const stale = await prisma.catalogDataset.findMany({
      where: { id: { notIn: Array.from(touchedIds) } },
      select: { id: true, title: true },
    });

    if (stale.length > 0) {
      await prisma.catalogDataset.deleteMany({
        where: { id: { in: stale.map((r) => r.id) } },
      });
      deleted.push(...stale);
    }
  }

  return NextResponse.json({
    mode,
    created,
    updated,
    skipped,
    errors,
    deleted,
  });
}
