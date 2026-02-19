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

interface InitiateSingleMoleculeUploadRequest {
  fingerprint: string;
  metadata: {
    title?: string;
    numMolecules: number;
    numGenes: number;
    platform?: string;
    description?: string;
  };
  manifest: any; // The manifest JSON to store
  files: Array<{
    key: string;
    size: number;
    contentType?: string;
  }>;
}

export async function POST(request: NextRequest) {
  try {
    const body: InitiateSingleMoleculeUploadRequest = await request.json();
    const { fingerprint, metadata, manifest, files } = body;

    // Validate required fields
    if (
      !fingerprint ||
      !metadata ||
      !manifest ||
      !files ||
      files.length === 0
    ) {
      return NextResponse.json(
        {
          error:
            "Missing required fields: fingerprint, metadata, manifest, or files",
        },
        { status: 400, headers: corsHeaders },
      );
    }

    // Debug: log what the API actually receives
    console.log("[API initiate] metadata:", JSON.stringify(metadata));
    console.log("[API initiate] manifest stats:", JSON.stringify(manifest?.statistics));
    console.log("[API initiate] manifest genes count:", manifest?.genes?.unique_gene_names?.length);

    // Resolve molecule/gene counts from metadata or manifest fallback
    const numMolecules =
      metadata.numMolecules ||
      manifest?.statistics?.total_molecules ||
      0;
    const numGenes =
      metadata.numGenes ||
      manifest?.statistics?.unique_genes ||
      manifest?.genes?.unique_gene_names?.length ||
      0;

    if (!numMolecules || !numGenes) {
      return NextResponse.json(
        { error: "Could not determine numMolecules or numGenes from metadata or manifest" },
        { status: 400, headers: corsHeaders },
      );
    }

    // Generate IDs
    const datasetId = `sm_${nanoid(10)}`;
    const uploadId = `up_${nanoid(10)}`;
    const expiresAt = new Date(Date.now() + 3600000); // 1 hour from now

    // Start transaction with extended timeout for many files
    const result = await prisma.$transaction(
      async (tx) => {
        // 1. Create dataset record
        const dataset = await tx.dataset.create({
          data: {
            id: datasetId,
            fingerprint,
            title: metadata.title || "Untitled Single Molecule Dataset",
            numCells: numMolecules, // Store molecule count in numCells field
            numGenes: numGenes,
            datasetType: "single_molecule",
            manifestJson: manifest,
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
      },
      {
        timeout: 60000, // 60 second timeout for large datasets with many genes
      },
    );

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
    console.error("Initiate single molecule upload error:", error);

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
