// app/api/datasets/check-duplicate/[fingerprint]/route.ts
// OR src/app/api/datasets/check-duplicate/[fingerprint]/route.ts
import { NextRequest, NextResponse } from "next/server";
import { prisma } from "@/lib/prisma";

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
  { params }: { params: Promise<{ fingerprint: string }> }
) {
  try {
    const { fingerprint } = await params;

    // Validate fingerprint format
    if (!fingerprint || fingerprint.length < 10 || fingerprint.length > 64) {
      return NextResponse.json(
        { error: "Invalid fingerprint format" },
        { status: 400, headers: corsHeaders }
      );
    }

    // Check if dataset with this fingerprint exists
    const existingDataset = await prisma.dataset.findUnique({
      where: { fingerprint },
      select: {
        id: true,
        title: true,
        status: true,
        numCells: true,
        numGenes: true,
        createdAt: true,
        completedAt: true,
        manifestUrl: true,
      },
    });

    if (existingDataset) {
      return NextResponse.json(
        {
          exists: true,
          dataset: {
            id: existingDataset.id,
            title: existingDataset.title,
            status: existingDataset.status,
            numCells: existingDataset.numCells,
            numGenes: existingDataset.numGenes,
            uploadedAt: existingDataset.createdAt,
            completedAt: existingDataset.completedAt,
            manifestUrl: existingDataset.manifestUrl,
          },
        },
        { headers: corsHeaders }
      );
    }

    return NextResponse.json(
      {
        exists: false,
        message: "No duplicate found. Safe to upload.",
      },
      { headers: corsHeaders }
    );
  } catch (error) {
    console.error("Check duplicate error:", error);
    return NextResponse.json(
      { error: "Internal server error" },
      { status: 500, headers: corsHeaders }
    );
  }
}
