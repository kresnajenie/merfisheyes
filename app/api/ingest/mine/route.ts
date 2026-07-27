import { NextResponse } from "next/server";

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
  processingProgress: true,
  errorMessage: true,
  batchJobId: true,
  createdAt: true,
  completedAt: true,
} as const;

export async function GET() {
  const { error, session } = await requireUser();

  if (error) return error;

  try {
    // Exclude UPLOADING entirely. The dataset currently transferring is shown
    // in the top upload banner (which owns the live progress), so a card here
    // would just duplicate it; and any *other* UPLOADING row is an abandoned
    // attempt (tab closed / interrupted before it finished) that would read as
    // a stuck "Uploading…" spinner. A dataset appears here once it reaches
    // QUEUED — i.e. as soon as the bytes are up and the job is submitted.
    const where = {
      ownerId: session.user.id,
      status: { not: "UPLOADING" as const },
    };

    let datasets = await prisma.dataset.findMany({
      where,
      select,
      orderBy: { createdAt: "desc" },
      take: 30,
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
        orderBy: { createdAt: "desc" },
        take: 30,
      });
    }

    return NextResponse.json({
      datasets: datasets.map((d) => {
        const viewerPath =
          d.datasetType === "single_molecule" ? "sm-viewer" : "viewer";

        return {
          id: d.id,
          title: d.title,
          status: d.status,
          datasetType: d.datasetType,
          progress: d.processingProgress,
          errorMessage: d.errorMessage,
          numCells: d.numCells,
          numGenes: d.numGenes,
          createdAt: d.createdAt,
          completedAt: d.completedAt,
          viewerUrl: d.status === "COMPLETE" ? `/${viewerPath}/${d.id}` : null,
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
