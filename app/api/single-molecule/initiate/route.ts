import { NextRequest, NextResponse } from "next/server";
import { nanoid } from "nanoid";

import { prisma } from "@/lib/prisma";
import { auth } from "@/lib/auth";
import { generatePresignedUploadUrl } from "@/lib/s3";
import { resolveDatasetOwner } from "@/lib/datasets/owner";

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
  /** Admins only: true → collective "admin" ownership instead of personal. */
  asAdmin?: boolean;
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

    // Uploads require sign-in so the dataset is always owned. Admins may opt
    // into collective "admin" ownership via asAdmin.
    const session = await auth();

    if (!session?.user) {
      return NextResponse.json(
        { error: "Sign in to upload" },
        { status: 401, headers: corsHeaders },
      );
    }
    const isAdmin =
      session.user.role === "ADMIN" || session.user.role === "SUPER_ADMIN";
    const { ownerId, adminOwned } = resolveDatasetOwner({
      isAdmin,
      asAdmin: !!body.asAdmin,
      userId: session.user.id,
    });

    // Generate IDs
    const datasetId = `sm_${nanoid(10)}`;
    const uploadId = `up_${nanoid(10)}`;
    const expiresAt = new Date(Date.now() + 3600000); // 1 hour from now

    // Start transaction (optimized for large file counts — one dataset with
    // hundreds of genes means hundreds of file rows, so batch-insert them and
    // keep presigned-URL generation out of the transaction entirely; the old
    // per-file create+presign loop inside the tx took ~33s for 530 files,
    // past Vercel's function timeout)
    await prisma.$transaction(
      async (tx) => {
        // 1. Create dataset record
        await tx.dataset.create({
          data: {
            id: datasetId,
            fingerprint,
            title: metadata.title || "Untitled Single Molecule Dataset",
            numCells: numMolecules, // Store molecule count in numCells field
            numGenes: numGenes,
            datasetType: "single_molecule",
            manifestJson: manifest,
            status: "UPLOADING",
            ownerId,
            adminOwned,
          },
        });

        // 2. Create upload session
        await tx.uploadSession.create({
          data: {
            id: uploadId,
            datasetId,
            totalFiles: files.length,
            completedFiles: 0,
            expiresAt,
          },
        });

        // 3. Batch create upload file records (much faster than individual creates)
        await tx.uploadFile.createMany({
          data: files.map((file) => ({
            uploadSessionId: uploadId,
            fileKey: file.key,
            fileSize: BigInt(file.size),
            status: "PENDING",
          })),
        });
      },
      {
        maxWait: 10000, // Maximum time to wait for a transaction slot (10s)
        timeout: 20000, // Maximum time for the transaction to complete (20s)
      },
    );

    // 4. Generate presigned URLs outside transaction (parallel processing)
    const uploadUrls: Record<string, any> = {};

    // Process in batches of 50 to avoid overwhelming S3
    const batchSize = 50;

    for (let i = 0; i < files.length; i += batchSize) {
      const batch = files.slice(i, i + batchSize);
      const batchResults = await Promise.all(
        batch.map(async (file) => {
          const s3Key = `datasets/${datasetId}/${file.key}`;
          const presignedUrl = await generatePresignedUploadUrl(
            s3Key,
            file.contentType || "application/octet-stream",
            3600, // 1 hour expiration
          );

          return { key: file.key, url: presignedUrl.url };
        }),
      );

      for (const { key, url } of batchResults) {
        uploadUrls[key] = url;
      }
    }

    // Return success response
    return NextResponse.json(
      {
        success: true,
        datasetId,
        uploadId,
        uploadUrls,
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
