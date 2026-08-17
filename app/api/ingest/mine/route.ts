import { NextRequest, NextResponse } from "next/server";

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

  const sortParam = request.nextUrl.searchParams.get("sort") ?? "date";
  const dirParam =
    request.nextUrl.searchParams.get("dir") === "asc" ? "asc" : "desc";
  const orderBy = {
    [SORT_FIELDS[sortParam] ?? "createdAt"]: dirParam,
  } as const;

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
    const where = {
      status: { not: "UPLOADING" as const },
      ...(isAdmin
        ? { OR: [{ ownerId: session.user.id }, { adminOwned: true }] }
        : { ownerId: session.user.id }),
    };

    let datasets = await prisma.dataset.findMany({
      where,
      select,
      orderBy,
      take: 200,
    });

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
        take: 200,
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

    return NextResponse.json({
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
