// app/api/datasets/initiate/route.ts
// OR src/app/api/datasets/initiate/route.ts
import { NextRequest, NextResponse } from "next/server";
import { nanoid } from "nanoid";

import { prisma } from "@/lib/prisma";
import { generatePresignedUploadUrl } from "@/lib/s3";

const corsHeaders = {
  "Access-Control-Allow-Origin": process.env.CORS_ORIGIN || "*",
  "Access-Control-Allow-Methods": "POST, OPTIONS",
  "Access-Control-Allow-Headers": "Content-Type",
};

export async function OPTIONS() {
  return NextResponse.json({}, { headers: corsHeaders });
}

interface InitiateUploadRequest {
  fingerprint: string;
  metadata: {
    title?: string;
    numCells: number;
    numGenes: number;
    platform?: string;
    description?: string;
  };
  files: Array<{
    key: string;
    size: number;
    contentType?: string;
  }>;
}

export async function POST(request: NextRequest) {
  try {
    const body: InitiateUploadRequest = await request.json();
    const { fingerprint, metadata, files } = body;

    // Validate required fields
    if (!fingerprint || !metadata || !files || files.length === 0) {
      return NextResponse.json(
        { error: "Missing required fields: fingerprint, metadata, or files" },
        { status: 400, headers: corsHeaders },
      );
    }

    if (!metadata.numCells || !metadata.numGenes) {
      return NextResponse.json(
        { error: "metadata.numCells and metadata.numGenes are required" },
        { status: 400, headers: corsHeaders },
      );
    }

    // Generate IDs
    const datasetId = `ds_${nanoid(10)}`;
    const uploadId = `up_${nanoid(10)}`;
    const expiresAt = new Date(Date.now() + 3600000); // 1 hour from now

    // Start transaction
    const result = await prisma.$transaction(async (tx) => {
      // 1. Create dataset record
      const dataset = await tx.dataset.create({
        data: {
          id: datasetId,
          fingerprint,
          title: metadata.title || "Untitled Dataset",
          numCells: metadata.numCells,
          numGenes: metadata.numGenes,
          datasetType: "single_cell",
          status: "UPLOADING",
        },
      });

      // 2. Create upload session
      const uploadSession = await tx.uploadSession.create({
        data: {
          id: uploadId,
          datasetId,
          totalFiles: files.length,
          completedFiles: 0,
          expiresAt,
        },
      });

      // 3. Create upload file records and generate presigned URLs
      const uploadUrls: Record<string, any> = {};

      for (const file of files) {
        // Create file record in database
        await tx.uploadFile.create({
          data: {
            uploadSessionId: uploadId,
            fileKey: file.key,
            fileSize: BigInt(file.size),
            status: "PENDING",
          },
        });

        // Generate S3 key (path in bucket)
        const s3Key = `datasets/${datasetId}/${file.key}`;

        // Generate presigned URL
        const presignedUrl = await generatePresignedUploadUrl(
          s3Key,
          file.contentType || "application/octet-stream",
          3600, // 1 hour expiration
        );

        uploadUrls[file.key] = presignedUrl.url;
      }

      return {
        dataset,
        uploadSession,
        uploadUrls,
      };
    });

    // Return success response
    return NextResponse.json(
      {
        success: true,
        datasetId,
        uploadId,
        uploadUrls: result.uploadUrls,
        expiresIn: 3600,
        expiresAt: expiresAt.toISOString(),
      },
      { headers: corsHeaders },
    );
  } catch (error: any) {
    console.error("Initiate upload error:", error);

    // Handle unique constraint violation (duplicate fingerprint)
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
      {
        error: "Internal server error",
        message: error.message,
      },
      { status: 500, headers: corsHeaders },
    );
  }
}
