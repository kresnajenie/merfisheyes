// app/api/datasets/[datasetId]/upload-progress/route.ts
// Bumps UploadSession.completedFiles by `delta`. Used by the zarr uploader
// in lieu of per-file UploadFile rows.
import { NextRequest, NextResponse } from "next/server";

import { prisma } from "@/lib/prisma";

const corsHeaders = {
  "Access-Control-Allow-Origin": process.env.CORS_ORIGIN || "*",
  "Access-Control-Allow-Methods": "POST, OPTIONS",
  "Access-Control-Allow-Headers": "Content-Type",
};

export async function OPTIONS() {
  return NextResponse.json({}, { headers: corsHeaders });
}

interface ProgressRequest {
  uploadId: string;
  delta: number;
}

export async function POST(
  request: NextRequest,
  { params }: { params: Promise<{ datasetId: string }> },
) {
  try {
    const { datasetId } = await params;
    const body: ProgressRequest = await request.json();
    const { uploadId, delta } = body;

    if (!uploadId || typeof delta !== "number" || delta <= 0) {
      return NextResponse.json(
        { error: "uploadId and positive delta are required" },
        { status: 400, headers: corsHeaders },
      );
    }

    const session = await prisma.uploadSession.findUnique({
      where: { id: uploadId },
      select: { datasetId: true, totalFiles: true },
    });

    if (!session || session.datasetId !== datasetId) {
      return NextResponse.json(
        { error: "Upload session not found" },
        { status: 404, headers: corsHeaders },
      );
    }

    const updated = await prisma.uploadSession.update({
      where: { id: uploadId },
      data: { completedFiles: { increment: delta } },
      select: { completedFiles: true, totalFiles: true },
    });

    return NextResponse.json(
      { success: true, ...updated },
      { headers: corsHeaders },
    );
  } catch (error: any) {
    console.error("upload-progress error:", error);

    return NextResponse.json(
      { error: "Internal server error", message: error.message },
      { status: 500, headers: corsHeaders },
    );
  }
}
