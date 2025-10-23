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
  { params }: { params: Promise<{ fingerprint: string }> },
) {
  try {
    const { fingerprint } = await params;

    if (!fingerprint || fingerprint.length < 10 || fingerprint.length > 64) {
      return NextResponse.json(
        { error: "Invalid fingerprint format" },
        { status: 400, headers: corsHeaders },
      );
    }

    // Find dataset with this fingerprint (only single molecule datasets)
    const existingDataset = await prisma.dataset.findFirst({
      where: {
        fingerprint,
        datasetType: "single_molecule",
        status: {
          in: ["COMPLETE", "PROCESSING"], // Only check completed/processing datasets
        },
      },
      select: {
        id: true,
        title: true,
        createdAt: true,
        numCells: true, // This stores molecule count for single molecule
        numGenes: true,
      },
    });

    if (!existingDataset) {
      return NextResponse.json(
        {
          exists: false,
          message: "No duplicate found. Safe to upload.",
        },
        { headers: corsHeaders },
      );
    }

    return NextResponse.json(
      {
        exists: true,
        dataset: {
          id: existingDataset.id,
          title: existingDataset.title,
          createdAt: existingDataset.createdAt.toISOString(),
          numMolecules: existingDataset.numCells,
          numGenes: existingDataset.numGenes,
        },
      },
      { headers: corsHeaders },
    );
  } catch (error: any) {
    console.error("Check duplicate error:", error);

    return NextResponse.json(
      {
        error: "Internal server error",
        message: error.message,
      },
      { status: 500, headers: corsHeaders },
    );
  }
}
