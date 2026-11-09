import { NextRequest, NextResponse } from "next/server";

import { auth } from "@/lib/auth";
import {
  hasMolecules,
  moleculeLinkFor,
  ownedDataset,
} from "@/lib/ingest/overlay";

/**
 * How the viewer opens one of the user's molecule datasets as an overlay:
 * the link a saved dataset's overlay would write for it (moleculeLinkFor).
 * A local preview, which has no dataset to write to, draws it from this to
 * check the two line up before saving.
 *
 *   GET → { linkColumn: "__all__", links: { __all__ }, source, linkFormat?, datasetId }
 */
export async function GET(
  _req: NextRequest,
  { params }: { params: Promise<{ datasetId: string }> },
) {
  const session = await auth();

  if (!session?.user?.id) {
    return NextResponse.json({ error: "Not signed in" }, { status: 401 });
  }
  const isAdmin =
    session.user.role === "ADMIN" || session.user.role === "SUPER_ADMIN";
  const { datasetId } = await params;
  const sm = await ownedDataset(datasetId, session.user.id, isAdmin);

  if (!sm || !(await hasMolecules(sm))) {
    return NextResponse.json(
      { error: "Pick a dataset of yours with molecules" },
      { status: 404 },
    );
  }

  return NextResponse.json(await moleculeLinkFor(datasetId));
}
