import { NextRequest, NextResponse } from "next/server";

import { requireUser } from "@/lib/admin-auth";
import { prisma } from "@/lib/prisma";

interface CompleteIngestRequest {
  uploadId: string;
}

export async function POST(
  request: NextRequest,
  { params }: { params: Promise<{ id: string }> },
) {
  const { error, session } = await requireUser();

  if (error) return error;

  try {
    const { id: datasetId } = await params;
    const { uploadId }: CompleteIngestRequest = await request.json();

    if (!uploadId) {
      return NextResponse.json(
        { error: "uploadId is required" },
        { status: 400 },
      );
    }

    const uploadSession = await prisma.uploadSession.findUnique({
      where: { id: uploadId },
      include: { dataset: true },
    });

    if (!uploadSession || uploadSession.datasetId !== datasetId) {
      return NextResponse.json(
        { error: "Upload session not found" },
        { status: 404 },
      );
    }

    if (uploadSession.dataset.ownerId !== session.user.id) {
      return NextResponse.json({ error: "Forbidden" }, { status: 403 });
    }

    if (uploadSession.completedFiles !== uploadSession.totalFiles) {
      return NextResponse.json(
        {
          error: "Upload incomplete",
          message: `Only ${uploadSession.completedFiles} of ${uploadSession.totalFiles} files uploaded`,
        },
        { status: 400 },
      );
    }

    // Raw bytes are in place. Hand off to the queue.
    // (SubmitJob to AWS Batch is Phase 3 — this phase stops at QUEUED.)
    await prisma.dataset.update({
      where: { id: datasetId },
      data: { status: "QUEUED" },
    });

    return NextResponse.json({
      success: true,
      message: "Raw upload complete; dataset queued for processing",
      datasetId,
    });
  } catch (err: any) {
    console.error("Ingest complete error:", err);

    return NextResponse.json(
      { error: "Internal server error", message: err.message },
      { status: 500 },
    );
  }
}
