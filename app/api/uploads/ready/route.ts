import { NextRequest, NextResponse } from "next/server";

import { requireUser } from "@/lib/admin-auth";
import { prisma } from "@/lib/prisma";
import { sendDatasetReadyEmail } from "@/lib/ses";

/**
 * A save that the browser uploads as it is — an image, a SpatialData store,
 * a processed folder — has finished: the signed-in user gets one "ready"
 * email for it, as a save processed on the server already does from
 * `finalize.ts` (docs/plans/PREVIEW-SAVE.md, D6). Such a save never reaches
 * the worker, and its completion routes are shared with saves that email
 * already, so the browser says when the whole save is done.
 *
 * `datasetIds` are the save's datasets, the one to open first first (a
 * folder's cells, which show its molecules and image over them). Only the
 * signed-in user's own finished datasets count, and the email only ever goes
 * to the signed-in user, whatever the request says.
 */
export async function POST(request: NextRequest) {
  const { error, session } = await requireUser();

  if (error) return error;

  let body: { datasetIds?: unknown };

  try {
    body = await request.json();
  } catch {
    return NextResponse.json({ error: "Invalid JSON" }, { status: 400 });
  }
  const ids = Array.isArray(body.datasetIds)
    ? body.datasetIds.filter((id): id is string => typeof id === "string")
    : [];

  if (!ids.length || ids.length > 10) {
    return NextResponse.json(
      { error: "datasetIds: 1 to 10 dataset ids" },
      { status: 400 },
    );
  }

  const found = await prisma.dataset.findMany({
    where: { id: { in: ids }, ownerId: session.user.id, status: "COMPLETE" },
    select: {
      id: true,
      title: true,
      datasetType: true,
      numCells: true,
      numGenes: true,
      platform: true,
    },
  });
  // In the order the save gave them
  const datasets = ids.flatMap((id) => found.filter((d) => d.id === id));

  // How many of the save's datasets counted, whatever became of the email
  const counted = datasets.length;

  if (!counted || !session.user.email) {
    return NextResponse.json({ success: true, emailed: false, counted });
  }
  const [first, ...rest] = datasets;

  try {
    await sendDatasetReadyEmail({
      email: session.user.email,
      datasetId: first.id,
      datasetType: first.datasetType,
      molecules: first.datasetType === "single_molecule",
      datasetName: first.title ?? undefined,
      metadata: {
        // An image counts neither cells nor genes
        numCells: first.numCells || undefined,
        numGenes: first.numGenes || undefined,
        platform: first.platform ?? first.datasetType ?? undefined,
      },
      alsoSaved: rest.map((d) => ({
        datasetId: d.id,
        datasetType: d.datasetType,
        datasetName: d.title ?? undefined,
      })),
    });
  } catch (e: any) {
    // The save is done either way; the email is a courtesy
    console.error("Ready email failed:", e?.message);

    return NextResponse.json({ success: true, emailed: false, counted });
  }

  return NextResponse.json({ success: true, emailed: true, counted });
}
