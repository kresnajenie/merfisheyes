import { NextRequest, NextResponse } from "next/server";

import { requireUser } from "@/lib/admin-auth";
import { prisma } from "@/lib/prisma";
import { completeMultipartUpload } from "@/lib/s3";

interface CompleteFileRequest {
  uploadId: string; // our UploadSession id (up_...)
  // Present only for multipart uploads: S3's UploadId + the ETag of every part.
  s3UploadId?: string;
  parts?: Array<{ partNumber: number; etag: string }>;
}

export async function POST(
  request: NextRequest,
  { params }: { params: Promise<{ id: string; key: string }> },
) {
  const { error, session } = await requireUser();

  if (error) return error;

  try {
    const { id: datasetId, key: encodedFileKey } = await params;
    const fileKey = decodeURIComponent(encodedFileKey);
    const body: CompleteFileRequest = await request.json();
    const { uploadId, s3UploadId, parts } = body;

    if (!uploadId) {
      return NextResponse.json(
        { error: "uploadId is required" },
        { status: 400 },
      );
    }

    // Verify the session belongs to this dataset and this user owns it.
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

    // Finalize a multipart upload before recording the file as complete.
    if (s3UploadId) {
      if (!parts || parts.length === 0) {
        return NextResponse.json(
          { error: "parts[] is required to complete a multipart upload" },
          { status: 400 },
        );
      }
      await completeMultipartUpload(
        `raw/${datasetId}/${fileKey}`,
        s3UploadId,
        parts,
      );
    }

    await prisma.$transaction(async (tx) => {
      const updated = await tx.uploadFile.updateMany({
        where: { uploadSessionId: uploadId, fileKey, status: { not: "COMPLETE" } },
        data: { status: "COMPLETE", uploadedAt: new Date() },
      });

      // Only bump the counter if this call actually flipped a pending row,
      // so a retried request can't over-count completedFiles.
      if (updated.count > 0) {
        await tx.uploadSession.update({
          where: { id: uploadId },
          data: { completedFiles: { increment: updated.count } },
        });
      }
    });

    return NextResponse.json({ success: true, message: "File marked as complete" });
  } catch (err: any) {
    console.error("Ingest mark file complete error:", err);

    return NextResponse.json(
      { error: "Internal server error", message: err.message },
      { status: 500 },
    );
  }
}
