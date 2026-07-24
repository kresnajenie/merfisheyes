import { NextRequest, NextResponse } from "next/server";

import { requireUser } from "@/lib/admin-auth";
import { prisma } from "@/lib/prisma";
import { isStale, reconcileWithBatch } from "@/lib/ingest/reconcile";

/**
 * Polling endpoint for ingestion progress (design §5.5). The browser polls this
 * while a dataset moves QUEUED → PROCESSING → COMPLETE/FAILED. (Supabase
 * Realtime can replace the polling later; this stays as the fallback.)
 */
export async function GET(
  request: NextRequest,
  { params }: { params: Promise<{ id: string }> },
) {
  const { error, session } = await requireUser();

  if (error) return error;

  try {
    const { id: datasetId } = await params;

    const select = {
      id: true,
      title: true,
      status: true,
      ownerId: true,
      datasetType: true,
      numCells: true,
      numGenes: true,
      processingProgress: true,
      errorMessage: true,
      batchJobId: true,
      createdAt: true,
      completedAt: true,
    };

    let dataset = await prisma.dataset.findUnique({
      where: { id: datasetId },
      select,
    });

    if (!dataset) {
      return NextResponse.json({ error: "Dataset not found" }, { status: 404 });
    }

    if (dataset.ownerId !== session.user.id) {
      return NextResponse.json({ error: "Forbidden" }, { status: 403 });
    }

    // A lost callback would otherwise leave this dataset at QUEUED/PROCESSING
    // forever, with the client polling indefinitely. Batch still knows the
    // job's fate, so ask it — but only once the row has gone quiet, so healthy
    // jobs cost no extra API calls. Re-read afterwards to return the truth in
    // the same response rather than on the next poll.
    if (isStale(dataset)) {
      const outcome = await reconcileWithBatch(dataset);

      if (outcome.changed) {
        dataset =
          (await prisma.dataset.findUnique({
            where: { id: datasetId },
            select,
          })) ?? dataset;
      }
    }

    const viewerPath =
      dataset.datasetType === "single_molecule" ? "sm-viewer" : "viewer";

    return NextResponse.json({
      id: dataset.id,
      title: dataset.title,
      status: dataset.status,
      progress: dataset.processingProgress,
      errorMessage: dataset.errorMessage,
      batchJobId: dataset.batchJobId,
      numCells: dataset.numCells,
      numGenes: dataset.numGenes,
      createdAt: dataset.createdAt,
      completedAt: dataset.completedAt,
      viewerUrl:
        dataset.status === "COMPLETE" ? `/${viewerPath}/${dataset.id}` : null,
    });
  } catch (err: any) {
    console.error("Ingest status error:", err);

    return NextResponse.json(
      { error: "Internal server error", message: err.message },
      { status: 500 },
    );
  }
}
