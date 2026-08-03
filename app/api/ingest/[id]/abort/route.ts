import { NextRequest, NextResponse } from "next/server";

import { requireUser } from "@/lib/admin-auth";
import { prisma } from "@/lib/prisma";
import {
  listInProgressMultipartUploads,
  abortMultipartUpload,
  deleteObjectsByPrefix,
} from "@/lib/s3";

export async function POST(
  request: NextRequest,
  { params }: { params: Promise<{ id: string }> },
) {
  const { error, session } = await requireUser();

  if (error) return error;

  try {
    const { id: datasetId } = await params;

    const dataset = await prisma.dataset.findUnique({
      where: { id: datasetId },
    });

    if (!dataset) {
      return NextResponse.json({ error: "Dataset not found" }, { status: 404 });
    }

    if (dataset.ownerId !== session.user.id) {
      return NextResponse.json({ error: "Forbidden" }, { status: 403 });
    }

    // Only an in-progress (pre-processing) upload can be aborted.
    if (dataset.status !== "UPLOADING" && dataset.status !== "QUEUED") {
      return NextResponse.json(
        {
          error: "Cannot abort",
          message: `Dataset is ${dataset.status}; only UPLOADING or QUEUED uploads can be aborted`,
        },
        { status: 409 },
      );
    }

    const prefix = `raw/${datasetId}/`;

    // 1. Abort any incomplete multipart uploads (releases their staged parts).
    //    Found server-side via ListMultipartUploads — no client state needed.
    //    Degrade gracefully: if enumeration fails (e.g. missing
    //    s3:ListBucketMultipartUploads permission), the S3 lifecycle rule that
    //    expires incomplete multipart uploads is the documented backstop
    //    (design §8.6/§10), so we still clean up objects + the DB row.
    let abortedMultipartUploads = 0;

    try {
      const inProgress = await listInProgressMultipartUploads(prefix);

      await Promise.all(
        inProgress.map((u) => abortMultipartUpload(u.key, u.uploadId)),
      );
      abortedMultipartUploads = inProgress.length;
    } catch (cleanupErr: any) {
      console.warn(
        `Ingest abort: could not enumerate/abort multipart uploads for ${prefix} (falling back to lifecycle rule):`,
        cleanupErr.message,
      );
    }

    // 2. Delete any fully-uploaded raw objects.
    const deletedObjects = await deleteObjectsByPrefix(prefix);

    // 3. Remove the dataset record (cascades to UploadSession + UploadFile).
    await prisma.dataset.delete({ where: { id: datasetId } });

    return NextResponse.json({
      success: true,
      message: "Upload aborted and cleaned up",
      datasetId,
      abortedMultipartUploads,
      deletedObjects,
    });
  } catch (err: any) {
    console.error("Ingest abort error:", err);

    return NextResponse.json(
      { error: "Internal server error", message: err.message },
      { status: 500 },
    );
  }
}
