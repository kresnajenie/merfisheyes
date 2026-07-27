import { NextRequest, NextResponse } from "next/server";
import { nanoid } from "nanoid";

import { requireUser } from "@/lib/admin-auth";
import { prisma } from "@/lib/prisma";
import {
  generatePresignedUploadUrl,
  createMultipartUpload,
  generatePresignedUploadPartUrls,
  computePartCount,
  MULTIPART_THRESHOLD,
  MULTIPART_PART_SIZE,
} from "@/lib/s3";

type DatasetKind = "single_cell" | "single_molecule";

interface InitiateIngestRequest {
  metadata?: {
    title?: string;
  };
  // The §5.8 processing pipeline spec: { kind, stages: {...} }. Stored verbatim
  // on Dataset.processingParams and consumed by the worker (Phase 3).
  processingParams: {
    kind: DatasetKind;
    stages?: Record<string, unknown>;
    [key: string]: unknown;
  };
  files: Array<{
    key: string; // relative path (folders preserve their path)
    size: number; // bytes
    contentType?: string;
  }>;
}

type UploadTarget =
  | { mode: "single"; url: string }
  | {
      mode: "multipart";
      s3UploadId: string;
      partSize: number;
      parts: Array<{ partNumber: number; url: string }>;
    };

export async function POST(request: NextRequest) {
  // Logged-in users only (design §2).
  const { error, session } = await requireUser();

  if (error) return error;

  try {
    const body: InitiateIngestRequest = await request.json();
    const { metadata, processingParams, files } = body;

    const kind = processingParams?.kind;

    if (
      !processingParams ||
      (kind !== "single_cell" && kind !== "single_molecule") ||
      !Array.isArray(files) ||
      files.length === 0
    ) {
      return NextResponse.json(
        {
          error:
            "Missing/invalid fields: processingParams.kind ('single_cell'|'single_molecule') and a non-empty files[] are required",
        },
        { status: 400 },
      );
    }

    if (files.some((f) => !f.key || typeof f.size !== "number" || f.size < 0)) {
      return NextResponse.json(
        { error: "Each file needs a non-empty key and a non-negative size" },
        { status: 400 },
      );
    }

    // Generate IDs. datasetType drives viewer routing later; the id prefix just
    // mirrors the existing sm_/ds_ convention.
    const datasetId =
      kind === "single_molecule" ? `sm_${nanoid(10)}` : `ds_${nanoid(10)}`;
    const uploadId = `up_${nanoid(10)}`;
    const expiresAt = new Date(Date.now() + 3600000); // 1 hour

    // The canonical fingerprint is computed by the worker post-processing
    // (README §6); we don't have processed data yet. Use a unique placeholder
    // (fingerprint is @unique and non-nullable) that the worker overwrites.
    const placeholderFingerprint = `raw_${datasetId}`;

    // DB writes in a transaction (fast); S3 calls happen after (mirrors the
    // datasets/initiate pattern — never hold a txn open across network calls).
    await prisma.$transaction(
      async (tx) => {
        await tx.dataset.create({
          data: {
            id: datasetId,
            fingerprint: placeholderFingerprint,
            title: metadata?.title || "Untitled Dataset",
            numCells: 0, // filled in by the worker post-processing
            numGenes: 0, // filled in by the worker post-processing
            datasetType: kind,
            status: "UPLOADING",
            ownerId: session.user.id,
            ingestSource: "server",
            processingParams: processingParams as object,
          },
        });

        await tx.uploadSession.create({
          data: {
            id: uploadId,
            datasetId,
            totalFiles: files.length,
            completedFiles: 0,
            expiresAt,
          },
        });

        await tx.uploadFile.createMany({
          data: files.map((file) => ({
            uploadSessionId: uploadId,
            fileKey: file.key,
            fileSize: BigInt(file.size),
            status: "PENDING" as const,
          })),
        });
      },
      { maxWait: 10000, timeout: 20000 },
    );

    // Presigned upload targets (outside the txn). Large files → multipart.
    const uploadUrls: Record<string, UploadTarget> = {};

    const batchSize = 25;

    for (let i = 0; i < files.length; i += batchSize) {
      const batch = files.slice(i, i + batchSize);
      const targets = await Promise.all(
        batch.map(async (file): Promise<[string, UploadTarget]> => {
          const s3Key = `raw/${datasetId}/${file.key}`;
          const contentType = file.contentType || "application/octet-stream";

          if (file.size > MULTIPART_THRESHOLD) {
            const { uploadId: s3UploadId } = await createMultipartUpload(
              s3Key,
              contentType,
            );
            const parts = await generatePresignedUploadPartUrls(
              s3Key,
              s3UploadId,
              computePartCount(file.size),
            );

            return [
              file.key,
              {
                mode: "multipart",
                s3UploadId,
                partSize: MULTIPART_PART_SIZE,
                parts,
              },
            ];
          }

          const presigned = await generatePresignedUploadUrl(
            s3Key,
            contentType,
          );

          return [file.key, { mode: "single", url: presigned.url }];
        }),
      );

      for (const [key, target] of targets) {
        uploadUrls[key] = target;
      }
    }

    return NextResponse.json({
      success: true,
      datasetId,
      uploadId,
      uploadUrls,
      expiresIn: 3600,
      expiresAt: expiresAt.toISOString(),
    });
  } catch (err: any) {
    console.error("Ingest initiate error:", err);

    return NextResponse.json(
      { error: "Internal server error", message: err.message },
      { status: 500 },
    );
  }
}
