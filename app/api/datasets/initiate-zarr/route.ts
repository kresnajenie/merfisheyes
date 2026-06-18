// app/api/datasets/initiate-zarr/route.ts
// Zarr-folder upload variant of /api/datasets/initiate.
// Returns ONE presigned POST policy that authorizes uploads to
// `datasets/${datasetId}/...` for an hour. The client then POSTs each file
// in the zarr tree against that single policy — no per-file signing.
import { createPresignedPost } from "@aws-sdk/s3-presigned-post";
import { NextRequest, NextResponse } from "next/server";
import { nanoid } from "nanoid";

import { prisma } from "@/lib/prisma";
import { S3_BUCKET, s3Client } from "@/lib/s3";

const corsHeaders = {
  "Access-Control-Allow-Origin": process.env.CORS_ORIGIN || "*",
  "Access-Control-Allow-Methods": "POST, OPTIONS",
  "Access-Control-Allow-Headers": "Content-Type",
};

export async function OPTIONS() {
  return NextResponse.json({}, { headers: corsHeaders });
}

interface InitiateZarrRequest {
  fingerprint: string;
  metadata: {
    title?: string;
    numCells: number;
    numGenes: number;
    platform?: string;
  };
  totalFiles: number;
  // Optional: caps the size of any single uploaded file. Defaults to 200 MB
  // which is well above any sensible single zarr chunk.
  maxFileSizeBytes?: number;
}

export async function POST(request: NextRequest) {
  try {
    const body: InitiateZarrRequest = await request.json();
    const { fingerprint, metadata, totalFiles } = body;

    if (!fingerprint || !metadata || !totalFiles) {
      return NextResponse.json(
        {
          error:
            "Missing required fields: fingerprint, metadata, or totalFiles",
        },
        { status: 400, headers: corsHeaders },
      );
    }
    if (!metadata.numCells || !metadata.numGenes) {
      return NextResponse.json(
        { error: "metadata.numCells and metadata.numGenes are required" },
        { status: 400, headers: corsHeaders },
      );
    }

    if (!S3_BUCKET) {
      return NextResponse.json(
        { error: "AWS_S3_BUCKET is not configured" },
        { status: 500, headers: corsHeaders },
      );
    }

    const datasetId = `ds_${nanoid(10)}`;
    const uploadId = `up_${nanoid(10)}`;
    const expiresAt = new Date(Date.now() + 3600_000);
    const maxFileSize = body.maxFileSizeBytes ?? 200 * 1024 * 1024;
    const keyPrefix = `datasets/${datasetId}/`;

    // Sign a single POST policy that authorizes uploads anywhere under the
    // dataset's S3 prefix until expiresAt.
    const { url, fields } = await createPresignedPost(s3Client, {
      Bucket: S3_BUCKET,
      // Placeholder key; client overrides per upload via the form's `key` field.
      Key: `${keyPrefix}\${filename}`,
      Conditions: [
        ["starts-with", "$key", keyPrefix],
        ["content-length-range", 0, maxFileSize],
      ],
      Expires: 3600,
    });

    // Create dataset row + upload session in a single transaction.
    // No per-file UploadFile rows for zarr — we only track aggregate progress.
    await prisma.$transaction(async (tx) => {
      await tx.dataset.create({
        data: {
          id: datasetId,
          fingerprint,
          title: metadata.title || "Untitled Zarr Dataset",
          numCells: metadata.numCells,
          numGenes: metadata.numGenes,
          datasetType: metadata.platform || "single_cell",
          formatVersion: "zarr",
          status: "UPLOADING",
        },
      });

      await tx.uploadSession.create({
        data: {
          id: uploadId,
          datasetId,
          totalFiles,
          completedFiles: 0,
          expiresAt,
        },
      });
    });

    return NextResponse.json(
      {
        success: true,
        datasetId,
        uploadId,
        keyPrefix,
        url,
        fields,
        expiresIn: 3600,
        expiresAt: expiresAt.toISOString(),
      },
      { headers: corsHeaders },
    );
  } catch (error: any) {
    console.error("Initiate zarr upload error:", error);

    if (error.code === "P2002" && error.meta?.target?.includes("fingerprint")) {
      return NextResponse.json(
        {
          error: "Duplicate dataset",
          message: "A dataset with this fingerprint already exists",
        },
        { status: 409, headers: corsHeaders },
      );
    }

    return NextResponse.json(
      { error: "Internal server error", message: error.message },
      { status: 500, headers: corsHeaders },
    );
  }
}
