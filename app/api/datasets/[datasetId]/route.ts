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
  { params }: { params: Promise<{ datasetId: string }> }
) {
  try {
    const { datasetId } = await params;

    // Validate dataset ID format
    if (!datasetId || !datasetId.startsWith("ds_")) {
      return NextResponse.json(
        { error: "Invalid dataset ID format" },
        { status: 400, headers: corsHeaders }
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
        { status: 404, headers: corsHeaders }
      );
    }

    // Check if dataset upload is complete
    if (dataset.status !== "COMPLETE") {
      return NextResponse.json(
        {
          error: "Dataset upload not complete",
          status: dataset.status,
          message:
            dataset.status === "UPLOADING"
              ? "Dataset is still being uploaded"
              : dataset.status === "PROCESSING"
                ? "Dataset is being processed"
                : "Dataset upload failed",
        },
        { status: 409, headers: corsHeaders }
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
      { headers: corsHeaders }
    );
  } catch (error) {
    console.error("Get dataset error:", error);
    return NextResponse.json(
      { error: "Internal server error" },
      { status: 500, headers: corsHeaders }
    );
  }
}
