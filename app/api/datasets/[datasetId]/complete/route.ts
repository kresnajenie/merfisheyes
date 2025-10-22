import { NextRequest, NextResponse } from "next/server";

import { prisma } from "@/lib/prisma";
import { generateManifestUrl } from "@/lib/s3";

const corsHeaders = {
  "Access-Control-Allow-Origin": process.env.CORS_ORIGIN || "*",
  "Access-Control-Allow-Methods": "POST, OPTIONS",
  "Access-Control-Allow-Headers": "Content-Type",
};

export async function OPTIONS() {
  return NextResponse.json({}, { headers: corsHeaders });
}

interface CompleteUploadRequest {
  uploadId: string;
  manifestKey?: string; // Optional: if manifest was uploaded separately
}

export async function POST(
  request: NextRequest,
  { params }: { params: Promise<{ datasetId: string }> },
) {
  try {
    const { datasetId } = await params;
    const body: CompleteUploadRequest = await request.json();
    const { uploadId, manifestKey } = body;

    // Validate required fields
    if (!uploadId) {
      return NextResponse.json(
        { error: "uploadId is required" },
        { status: 400, headers: corsHeaders },
      );
    }

    // Verify dataset exists and belongs to this upload session
    const dataset = await prisma.dataset.findUnique({
      where: { id: datasetId },
      include: {
        uploadSessions: {
          where: { id: uploadId },
          include: {
            files: true,
          },
        },
      },
    });

    if (!dataset) {
      return NextResponse.json(
        { error: "Dataset not found" },
        { status: 404, headers: corsHeaders },
      );
    }

    const uploadSession = dataset.uploadSessions[0];

    if (!uploadSession) {
      return NextResponse.json(
        { error: "Upload session not found" },
        { status: 404, headers: corsHeaders },
      );
    }

    // Check if all files were uploaded
    const completedFiles = uploadSession.files.filter(
      (f) => f.status === "COMPLETE",
    ).length;

    if (completedFiles !== uploadSession.totalFiles) {
      return NextResponse.json(
        {
          error: "Not all files have been uploaded",
          uploaded: completedFiles,
          total: uploadSession.totalFiles,
        },
        { status: 400, headers: corsHeaders },
      );
    }

    // Generate manifest URL if manifest was uploaded
    let manifestUrl: string | null = null;

    if (manifestKey) {
      manifestUrl = await generateManifestUrl(datasetId);
    }

    // Mark dataset as complete
    const updatedDataset = await prisma.dataset.update({
      where: { id: datasetId },
      data: {
        status: "COMPLETE",
        completedAt: new Date(),
        manifestUrl: manifestUrl,
      },
    });

    // Get the base URL for the dataset
    const baseUrl =
      process.env.NEXT_PUBLIC_APP_URL || process.env.VERCEL_URL
        ? `https://${process.env.VERCEL_URL}`
        : "http://localhost:3000";

    return NextResponse.json(
      {
        success: true,
        datasetId: updatedDataset.id,
        status: updatedDataset.status,
        completedAt: updatedDataset.completedAt,
        viewUrl: `${baseUrl}/view/${datasetId}`,
        shareUrl: `${baseUrl}/view/${datasetId}`,
      },
      { headers: corsHeaders },
    );
  } catch (error: any) {
    console.error("Complete upload error:", error);

    return NextResponse.json(
      {
        error: "Internal server error",
        message: error.message,
      },
      { status: 500, headers: corsHeaders },
    );
  }
}
