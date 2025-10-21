import { NextRequest, NextResponse } from "next/server";
import { prisma } from "@/lib/prisma";
import { generateManifestUrl } from "@/lib/s3";

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
  { params }: { params: Promise<{ id: string }> }
) {
  try {
    const { id: datasetId } = await params;

    // Validate dataset ID format for single molecule
    if (!datasetId || !datasetId.startsWith("sm_")) {
      return NextResponse.json(
        { error: "Invalid single molecule dataset ID format" },
        { status: 400, headers: corsHeaders }
      );
    }

    // Fetch dataset from database
    const dataset = await prisma.dataset.findUnique({
      where: { id: datasetId },
      select: {
        id: true,
        title: true,
        fingerprint: true,
        numCells: true, // stores molecule count for single molecule
        numGenes: true,
        datasetType: true,
        manifestJson: true,
        status: true,
        createdAt: true,
        completedAt: true,
      },
    });

    if (!dataset) {
      return NextResponse.json(
        { error: "Dataset not found" },
        { status: 404, headers: corsHeaders }
      );
    }

    // Verify this is a single molecule dataset
    if (dataset.datasetType !== "single_molecule") {
      return NextResponse.json(
        { error: "Not a single molecule dataset" },
        { status: 400, headers: corsHeaders }
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
        { status: 202, headers: corsHeaders } // 202 Accepted - processing not complete
      );
    }

    // Generate presigned URL for manifest (valid for 12 hours)
    const manifestUrl = await generateManifestUrl(datasetId, 12 * 3600);

    // Calculate expiration time
    const expiresIn = 12 * 3600; // 12 hours in seconds
    const expiresAt = new Date(Date.now() + expiresIn * 1000);

    return NextResponse.json(
      {
        id: dataset.id,
        title: dataset.title,
        fingerprint: dataset.fingerprint,
        numMolecules: dataset.numCells, // Return as numMolecules for clarity
        numGenes: dataset.numGenes,
        status: dataset.status,
        datasetType: dataset.datasetType,
        manifest: dataset.manifestJson, // Include manifest JSON from database
        manifestUrl, // Presigned URL to download manifest.json.gz from S3
        createdAt: dataset.createdAt,
        completedAt: dataset.completedAt,
        expiresIn,
        expiresAt: expiresAt.toISOString(),
      },
      { headers: corsHeaders }
    );
  } catch (error) {
    console.error("Get single molecule dataset error:", error);
    return NextResponse.json(
      { error: "Internal server error" },
      { status: 500, headers: corsHeaders }
    );
  }
}
