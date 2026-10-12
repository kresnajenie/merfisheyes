import { NextRequest, NextResponse } from "next/server";

import { ownerOrAdminError, requireUser } from "@/lib/admin-auth";
import {
  deletableDatasetSelect,
  deleteDatasetCompletely,
  exploreCardFor,
  overlayOf,
} from "@/lib/ingest/delete-dataset";
import { prisma } from "@/lib/prisma";

/**
 * DELETE /api/ingest/[id] — permanently delete one of your datasets, in any
 * state (a job still processing it is stopped).
 *
 * `?withOverlay=1` also deletes the single-molecule dataset overlaid on a
 * single-cell one, when you own that too. Without it the molecule dataset is
 * kept; deleting a molecule dataset always unlinks it from the cell datasets
 * that overlay it.
 *
 * 409 `on_explore` when the dataset (or the overlay going with it) is shown
 * on Explore or awaiting review: it has to be withdrawn first.
 */
export async function DELETE(
  request: NextRequest,
  { params }: { params: Promise<{ id: string }> },
) {
  const { error, session } = await requireUser();

  if (error) return error;

  try {
    const { id } = await params;
    const dataset = await prisma.dataset.findUnique({
      where: { id },
      select: deletableDatasetSelect,
    });

    if (!dataset) {
      return NextResponse.json({ error: "Dataset not found" }, { status: 404 });
    }

    const forbidden = ownerOrAdminError(dataset, session);

    if (forbidden) return forbidden;

    const targets = [dataset];

    if (request.nextUrl.searchParams.get("withOverlay") === "1") {
      const smId = await overlayOf(dataset);
      const sm = smId
        ? await prisma.dataset.findUnique({
            where: { id: smId },
            select: deletableDatasetSelect,
          })
        : null;

      if (sm && !ownerOrAdminError(sm, session)) targets.push(sm);
    }

    // Check everything before deleting anything
    for (const target of targets) {
      const card = await exploreCardFor(target.id);

      if (card) {
        return NextResponse.json(
          {
            error: "on_explore",
            message: `"${target.title ?? target.id}" is on Explore (or awaiting review) as "${card.title}". Withdraw it from Explore first.`,
          },
          { status: 409 },
        );
      }
    }

    const unlinkedFrom: string[] = [];

    for (const target of targets) {
      unlinkedFrom.push(
        ...(await deleteDatasetCompletely(target)).unlinkedFrom,
      );
    }

    return NextResponse.json({
      success: true,
      deleted: targets.map((t) => t.id),
      unlinkedFrom,
    });
  } catch (err: any) {
    console.error("Dataset delete error:", err);

    return NextResponse.json(
      { error: "Internal server error", message: err.message },
      { status: 500 },
    );
  }
}
