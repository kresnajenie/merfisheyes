import { NextRequest } from "next/server";

import { sendDatasetReadyEmail } from "@/lib/ses";
import { prisma } from "@/lib/prisma";
import { requireUploadOwner } from "@/lib/admin-auth";

// POST /api/send-email — "your dataset is ready", sent by the browser after a
// save. Only the dataset's signed-in owner may ask, and the mail goes to
// their own account address: an `email` in the body is ignored (it used to
// be the recipient, which let anyone send our mail to anyone).
export async function POST(req: NextRequest) {
  try {
    const { datasetId, datasetName, metadata } = await req.json();

    if (!datasetId) {
      return Response.json({ error: "datasetId is required" }, { status: 400 });
    }
    const owner = await requireUploadOwner(String(datasetId));

    if (owner.error) return owner.error;
    const email = owner.session.user.email;

    if (!email) {
      return Response.json(
        { error: "Your account has no email address" },
        { status: 400 },
      );
    }

    // The link opens the viewer of what the dataset is
    const dataset = await prisma.dataset.findUnique({
      where: { id: String(datasetId) },
      select: { datasetType: true },
    });

    await sendDatasetReadyEmail({
      email,
      datasetId,
      datasetType: dataset?.datasetType,
      datasetName,
      metadata,
    });

    return Response.json({ success: true });
  } catch (error) {
    console.error("SES Error:", error);

    return Response.json(
      {
        error: error instanceof Error ? error.message : "Failed to send email",
      },
      { status: 500 },
    );
  }
}
