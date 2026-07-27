import { NextRequest, NextResponse } from "next/server";

import { requireUser } from "@/lib/admin-auth";
import { isBatchConfigured, submitIngestJob } from "@/lib/batch";
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
      include: { dataset: { include: { owner: true } } },
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

    // Raw bytes are in place — mark QUEUED before submitting so the state is
    // correct even if SubmitJob fails.
    await prisma.dataset.update({
      where: { id: datasetId },
      data: { status: "QUEUED" },
    });

    // Fire-and-forget the Batch job: SubmitJob returns immediately and every
    // status update arrives via /api/ingest/[id]/callback. Never block this
    // route on processing (Vercel functions are short-lived).
    let batchJobId: string | null = null;
    let submitError: string | null = null;

    if (isBatchConfigured()) {
      try {
        batchJobId = await submitIngestJob({
          datasetId,
          kind:
            uploadSession.dataset.datasetType === "single_molecule"
              ? "single_molecule"
              : "single_cell",
          processingParams: uploadSession.dataset.processingParams,
          computeTier: uploadSession.dataset.owner?.computeTier,
        });
        await prisma.dataset.update({
          where: { id: datasetId },
          data: { batchJobId },
        });
      } catch (e: any) {
        // The upload itself succeeded; leave the dataset QUEUED so the job can
        // be resubmitted rather than losing the raw data.
        submitError = e?.message || "SubmitJob failed";
        console.error(`Ingest complete: SubmitJob failed for ${datasetId}:`, e);
      }
    } else {
      submitError = "Batch not configured (CALLBACK_SECRET / callback base URL)";
      console.warn(`Ingest complete: ${submitError}; ${datasetId} stays QUEUED`);
    }

    return NextResponse.json({
      success: true,
      message: batchJobId
        ? "Raw upload complete; processing job submitted"
        : "Raw upload complete; dataset queued (job not submitted)",
      datasetId,
      batchJobId,
      submitError,
    });
  } catch (err: any) {
    console.error("Ingest complete error:", err);

    return NextResponse.json(
      { error: "Internal server error", message: err.message },
      { status: 500 },
    );
  }
}
