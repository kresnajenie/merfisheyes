// app/api/datasets/[datasetId]/route.ts
import { NextRequest, NextResponse } from "next/server";

import { prisma } from "@/lib/prisma";
import { generateDatasetUrls } from "@/lib/s3";

const corsHeaders = {
  "Access-Control-Allow-Origin": process.env.CORS_ORIGIN || "*",
  "Access-Control-Allow-Methods": "GET, OPTIONS",
  "Access-Control-Allow-Headers": "Content-Type",
};

export async function OPTIONS() {
  return NextResponse.json({}, { headers: corsHeaders });
}

export async function GET(
  request: NextRequest,
  { params }: { params: Promise<{ datasetId: string }> },
) {
  try {
    const { datasetId } = await params;

    // Validate dataset ID format
    if (!datasetId || !datasetId.startsWith("ds_")) {
      return NextResponse.json(
        { error: "Invalid dataset ID format" },
        { status: 400, headers: corsHeaders },
      );
    }

    // Fetch dataset from database
    const dataset = await prisma.dataset.findUnique({
      where: { id: datasetId },
      include: {
        uploadSessions: {
          include: {
            files: true,
          },
          orderBy: { createdAt: "desc" },
          take: 1, // Get most recent upload session
        },
      },
    });

    if (!dataset) {
      return NextResponse.json(
        { error: "Dataset not found" },
        { status: 404, headers: corsHeaders },
      );
    }

    // Check if dataset upload is complete
    if (dataset.status !== "COMPLETE") {
      const statusMessages = {
        UPLOADING: "Dataset is still being uploaded",
        PROCESSING: "Dataset is being processed",
        FAILED: "Dataset upload failed",
      };

      return NextResponse.json(
        {
          error: "Dataset not ready",
          status: dataset.status,
          message:
            statusMessages[dataset.status as keyof typeof statusMessages] ||
            "Dataset is not available",
        },
        { status: 202, headers: corsHeaders }, // 202 Accepted - processing not complete
      );
    }

    // Get file list from most recent upload session
    const uploadSession = dataset.uploadSessions[0];
    const fileKeys =
      uploadSession?.files
        .filter((f) => f.status === "COMPLETE")
        .map((f) => f.fileKey) || [];

    // Generate presigned URLs for all files (valid for 12 hours)
    const fileUrls = await generateDatasetUrls(datasetId, fileKeys, 12 * 3600);

    // Calculate expiration time
    const expiresIn = 12 * 3600; // 12 hours in seconds
    const expiresAt = new Date(Date.now() + expiresIn * 1000);

    return NextResponse.json(
      {
        id: dataset.id,
        title: dataset.title,
        fingerprint: dataset.fingerprint,
        numCells: dataset.numCells,
        numGenes: dataset.numGenes,
        status: dataset.status,
        manifestUrl: dataset.manifestUrl,
        createdAt: dataset.createdAt,
        completedAt: dataset.completedAt,
        files: fileUrls,
        expiresIn,
        expiresAt: expiresAt.toISOString(),
      },
      { headers: corsHeaders },
    );
  } catch (error) {
    console.error("Get dataset error:", error);

    return NextResponse.json(
      { error: "Internal server error" },
      { status: 500, headers: corsHeaders },
    );
  }
}
