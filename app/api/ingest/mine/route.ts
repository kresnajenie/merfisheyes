import { NextRequest, NextResponse } from "next/server";
import { Prisma, DatasetStatus } from "@prisma/client";

import { requireUser } from "@/lib/admin-auth";
import { prisma } from "@/lib/prisma";
import { isStale, reconcileWithBatch } from "@/lib/ingest/reconcile";

/**
 * List the signed-in user's own datasets for the "Your uploads" strip on
 * /explore. Mirrors the single-dataset status route, but returns the whole set
 * so the page can show in-progress uploads next to finished ones.
 *
 * Self-heals: any in-progress row that has gone quiet (a lost callback) is
 * reconciled against Batch here too, so the strip never shows a spinner that
 * can't resolve.
 */
const select = {
  id: true,
  title: true,
  status: true,
  datasetType: true,
  numCells: true,
  numGenes: true,
  viewCount: true,
  s3BaseUrl: true,
  ingestSource: true,
  adminOwned: true,
  thumbnailUrl: true,
  processingProgress: true,
  errorMessage: true,
  batchJobId: true,
  createdAt: true,
  completedAt: true,
  // Catalog-style metadata (owner-editable, carried into a submission).
  description: true,
  species: true,
  disease: true,
  tissue: true,
  platform: true,
  institute: true,
  tags: true,
  externalLink: true,
  publicationLink: true,
  metadata: true,
} as const;

const SORT_FIELDS: Record<string, "createdAt" | "numCells" | "title"> = {
  date: "createdAt",
  molecules: "numCells",
  name: "title",
};

export async function GET(request: NextRequest) {
  const { error, session } = await requireUser();

  if (error) return error;

  const params = request.nextUrl.searchParams;
  const sortParam = params.get("sort") ?? "date";
  const dirParam = params.get("dir") === "asc" ? "asc" : "desc";
  const orderBy = {
    [SORT_FIELDS[sortParam] ?? "createdAt"]: dirParam,
  } as const;

  // Filters (all optional). Mirrors the Explore pill filters plus owner-only
  // Type / Status.
  const search = params.get("search")?.trim() ?? "";
  const species = params.get("species") ?? "";
  const tissue = params.get("tissue") ?? "";
  const platform = params.get("platform") ?? "";
  const type = params.get("type") ?? ""; // single_cell | single_molecule
  const statusFilter = params.get("status") ?? "";

  // Pagination is opt-in: only when a page/limit is supplied (the Explore
  // "Your uploads" strip calls with no params and expects the full list).
  const paginate = params.has("page") || params.has("limit");
  const page = Math.max(1, Number(params.get("page") ?? "1"));
  const limit = paginate
    ? Math.min(100, Math.max(1, Number(params.get("limit") ?? "20")))
    : 200;
  const skip = (page - 1) * limit;

  try {
    // Exclude UPLOADING entirely. The dataset currently transferring is shown
    // in the top upload banner (which owns the live progress), so a card here
    // would just duplicate it; and any *other* UPLOADING row is an abandoned
    // attempt (tab closed / interrupted before it finished) that would read as
    // a stuck "Uploading…" spinner. A dataset appears here once it reaches
    // QUEUED — i.e. as soon as the bytes are up and the job is submitted.
    // Personal datasets, plus admin-owned (shared) ones for admins.
    const isAdmin =
      session.user.role === "ADMIN" || session.user.role === "SUPER_ADMIN";
    // Ownership scope, reused for the where clause and the facet queries.
    const ownerScope: Prisma.DatasetWhereInput = {
      status: { not: "UPLOADING" as const },
      ...(isAdmin
        ? { OR: [{ ownerId: session.user.id }, { adminOwned: true }] }
        : { ownerId: session.user.id }),
    };

    const insensitive = Prisma.QueryMode.insensitive;
    const filterAnd: Prisma.DatasetWhereInput[] = [];

    if (type) filterAnd.push({ datasetType: type });
    if (statusFilter && statusFilter in DatasetStatus)
      filterAnd.push({ status: statusFilter as DatasetStatus });
    if (species)
      filterAnd.push({ species: { equals: species, mode: insensitive } });
    if (tissue)
      filterAnd.push({ tissue: { equals: tissue, mode: insensitive } });
    if (platform)
      filterAnd.push({ platform: { equals: platform, mode: insensitive } });
    if (search) {
      filterAnd.push({
        OR: [
          { title: { contains: search, mode: insensitive } },
          { description: { contains: search, mode: insensitive } },
          { species: { contains: search, mode: insensitive } },
          { tissue: { contains: search, mode: insensitive } },
          { platform: { contains: search, mode: insensitive } },
          { institute: { contains: search, mode: insensitive } },
          { tags: { has: search } },
        ],
      });
    }

    const where: Prisma.DatasetWhereInput = filterAnd.length
      ? { ...ownerScope, AND: filterAnd }
      : ownerScope;

    let [datasets, total, facets] = await Promise.all([
      prisma.dataset.findMany({ where, select, orderBy, skip, take: limit }),
      prisma.dataset.count({ where }),
      getFacets(ownerScope),
    ]);

    const stale = datasets.filter((d) => isStale(d));

    if (stale.length > 0) {
      await Promise.all(
        stale.map((d) =>
          reconcileWithBatch(d).catch((e) =>
            console.error(`mine: reconcile ${d.id} failed:`, e?.message),
          ),
        ),
      );
      datasets = await prisma.dataset.findMany({
        where,
        select,
        orderBy,
        skip,
        take: limit,
      });
    }

    // Attach the community-submission status (if any) for each dataset, so the
    // account cards can show "Pending review" / "Published" / "Rejected".
    const submissions = await prisma.catalogDataset.findMany({
      where: {
        isCommunity: true,
        sourceDatasetId: { in: datasets.map((d) => d.id) },
      },
      select: {
        id: true,
        sourceDatasetId: true,
        reviewStatus: true,
        isPublished: true,
        reviewNote: true,
      },
    });
    const submissionByDataset = new Map(
      submissions.map((s) => [s.sourceDatasetId, s]),
    );

    // Project memberships (the caller's own projects) for the card badge.
    const memberships = await prisma.projectDataset.findMany({
      where: {
        datasetId: { in: datasets.map((d) => d.id) },
        project: { ownerId: session.user.id },
      },
      select: { datasetId: true, project: { select: { title: true } } },
    });
    const projectNamesByDataset = new Map<string, string[]>();

    for (const m of memberships) {
      const list = projectNamesByDataset.get(m.datasetId) ?? [];

      list.push(m.project.title);
      projectNamesByDataset.set(m.datasetId, list);
    }

    return NextResponse.json({
      total,
      page,
      limit,
      filters: facets,
      datasets: datasets.map((d) => {
        const viewerPath =
          d.datasetType === "single_molecule" ? "sm-viewer" : "viewer";
        // S3-registered rows load via the from-s3 URL flow, not by id.
        const viewerUrl =
          d.status !== "COMPLETE"
            ? null
            : d.ingestSource === "s3_registered" && d.s3BaseUrl
              ? `/${viewerPath}/from-s3?url=${encodeURIComponent(d.s3BaseUrl)}`
              : `/${viewerPath}/${d.id}`;
        const sub = submissionByDataset.get(d.id);

        return {
          id: d.id,
          title: d.title,
          status: d.status,
          datasetType: d.datasetType,
          progress: d.processingProgress,
          errorMessage: d.errorMessage,
          numCells: d.numCells,
          numGenes: d.numGenes,
          viewCount: d.viewCount,
          ingestSource: d.ingestSource,
          adminOwned: d.adminOwned,
          thumbnailUrl: d.thumbnailUrl,
          createdAt: d.createdAt,
          completedAt: d.completedAt,
          viewerUrl,
          // Owner-editable metadata.
          description: d.description,
          species: d.species,
          disease: d.disease,
          tissue: d.tissue,
          platform: d.platform,
          institute: d.institute,
          tags: d.tags,
          externalLink: d.externalLink,
          publicationLink: d.publicationLink,
          metadata: d.metadata,
          projectNames: projectNamesByDataset.get(d.id) ?? [],
          submission: sub
            ? {
                catalogId: sub.id,
                reviewStatus: sub.reviewStatus,
                isPublished: sub.isPublished,
                reviewNote: sub.reviewNote,
              }
            : null,
        };
      }),
    });
  } catch (err: any) {
    console.error("Ingest mine error:", err);

    return NextResponse.json(
      { error: "Internal server error", message: err.message },
      { status: 500 },
    );
  }
}

// Distinct species / tissue / platform across the owner's datasets (ignoring
// the active filters, so the pills always offer every value). Mirrors the
// Explore facets shape.
async function getFacets(ownerScope: Prisma.DatasetWhereInput) {
  const [speciesRaw, tissueRaw, platformRaw] = await Promise.all([
    prisma.dataset.findMany({
      where: { ...ownerScope, species: { not: null } },
      select: { species: true },
      distinct: ["species"],
      orderBy: { species: "asc" },
    }),
    prisma.dataset.findMany({
      where: { ...ownerScope, tissue: { not: null } },
      select: { tissue: true },
      distinct: ["tissue"],
      orderBy: { tissue: "asc" },
    }),
    prisma.dataset.findMany({
      where: { ...ownerScope, platform: { not: null } },
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
