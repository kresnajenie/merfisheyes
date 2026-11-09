// app/api/datasets/check-duplicate/[fingerprint]/route.ts
// OR src/app/api/datasets/check-duplicate/[fingerprint]/route.ts
import { NextRequest, NextResponse } from "next/server";

import { prisma } from "@/lib/prisma";
import { isAdminRole, requireUser } from "@/lib/admin-auth";
import { ownedDatasetsWhere } from "@/lib/datasets/owner";

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
    // Only among the caller's own datasets: answering for anyone's told a
    // stranger the id and title of a dataset whose file they could hash.
    const { error, session } = await requireUser();

    if (error) return error;
    const { fingerprint } = await params;

    // Validate fingerprint format
    if (!fingerprint || fingerprint.length < 10 || fingerprint.length > 64) {
      return NextResponse.json(
        { error: "Invalid fingerprint format" },
        { status: 400, headers: corsHeaders },
      );
    }

    // Check if dataset with this fingerprint exists
    // findFirst, not findUnique: fingerprint is no longer a unique key.
    const existingDataset = await prisma.dataset.findFirst({
      where: {
        fingerprint,
        ...ownedDatasetsWhere(session.user.id, isAdminRole(session.user.role)),
      },
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
        { headers: corsHeaders },
      );
    }

    return NextResponse.json(
      {
        exists: false,
        message: "No duplicate found. Safe to upload.",
      },
      { headers: corsHeaders },
    );
  } catch (error) {
    console.error("Check duplicate error:", error);

    return NextResponse.json(
      { error: "Internal server error" },
      { status: 500, headers: corsHeaders },
    );
  }
}
