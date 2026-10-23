import { NextRequest, NextResponse } from "next/server";

import { requireUser } from "@/lib/admin-auth";
import { prisma } from "@/lib/prisma";
import { sendDatasetFailedEmail } from "@/lib/ses";
import { classifyFailure } from "@/lib/ingest/error-classification";
import { markDatasetFailed } from "@/lib/ingest/mark-failed";

/**
 * An upload that failed in the browser (Upload & Save, or the transfer of a
 * server upload): the signed-in user gets the reason by email, as a failed
 * server processing run already does (kresnajenie/merfisheyes#166).
 *
 * When the upload had already made its dataset — theirs, and not finished —
 * the dataset is marked FAILED with the reason (`markDatasetFailed`, which
 * also sends the email), so it no longer sits at UPLOADING and the admin
 * failures page shows it. The email only ever goes to the signed-in user.
 */
export async function POST(request: NextRequest) {
  const { error, session } = await requireUser();

  if (error) return error;

  let body: { datasetId?: unknown; title?: unknown; error?: unknown };

  try {
    body = await request.json();
  } catch {
    return NextResponse.json({ error: "Invalid JSON" }, { status: 400 });
  }
  const message =
    typeof body.error === "string" && body.error.trim()
      ? body.error.slice(0, 2000)
      : "Upload failed";
  const title = typeof body.title === "string" ? body.title.slice(0, 200) : "";

  if (typeof body.datasetId === "string" && body.datasetId) {
    const dataset = await prisma.dataset.findUnique({
      where: { id: body.datasetId },
      select: { ownerId: true, status: true },
    });

    if (
      dataset?.ownerId === session.user.id &&
      dataset.status !== "COMPLETE" &&
      dataset.status !== "FAILED"
    ) {
      await markDatasetFailed(body.datasetId, message);

      return NextResponse.json({ success: true, marked: true });
    }
  }

  // No dataset of theirs to mark (the upload failed before it was made):
  // the email alone
  if (session.user.email) {
    try {
      await sendDatasetFailedEmail({
        email: session.user.email,
        datasetId: "none — the upload stopped before the dataset was created",
        datasetName: title || undefined,
        errorMessage: message,
        userHint: classifyFailure(message).userHint,
      });
    } catch (e: any) {
      console.error("Upload failure email failed:", e?.message);
    }
  }

  return NextResponse.json({ success: true, marked: false });
}
