// app/api/datasets/[datasetId]/files/complete/route.ts
// Batched form of files/[fileKey]/complete: marks many UploadFile rows
// COMPLETE in one request, so a many-thousand-file upload doesn't cost one
// round trip per file.
import { NextRequest, NextResponse } from "next/server";

import { prisma } from "@/lib/prisma";

const corsHeaders = {
  "Access-Control-Allow-Origin": process.env.CORS_ORIGIN || "*",
  "Access-Control-Allow-Methods": "POST, OPTIONS",
  "Access-Control-Allow-Headers": "Content-Type",
};

const MAX_KEYS_PER_REQUEST = 1000;

export async function OPTIONS() {
  return NextResponse.json({}, { headers: corsHeaders });
}

export async function POST(
  request: NextRequest,
  { params }: { params: Promise<{ datasetId: string }> },
) {
  try {
    const { datasetId } = await params;
    const { uploadId, fileKeys } = await request.json();

    if (
      !uploadId ||
      !Array.isArray(fileKeys) ||
      fileKeys.length === 0 ||
      fileKeys.length > MAX_KEYS_PER_REQUEST
    ) {
      return NextResponse.json(
        {
          error: `uploadId and 1-${MAX_KEYS_PER_REQUEST} fileKeys are required`,
        },
        { status: 400, headers: corsHeaders },
      );
    }

    const uploadSession = await prisma.uploadSession.findUnique({
      where: { id: uploadId },
      select: { datasetId: true },
    });

    if (!uploadSession || uploadSession.datasetId !== datasetId) {
      return NextResponse.json(
        { error: "Upload session not found" },
        { status: 404, headers: corsHeaders },
      );
    }

    const updated = await prisma.uploadFile.updateMany({
      where: { uploadSessionId: uploadId, fileKey: { in: fileKeys } },
      data: { status: "COMPLETE" },
    });

    return NextResponse.json(
      { success: true, updated: updated.count },
      { headers: corsHeaders },
    );
  } catch (error: any) {
    console.error("Mark files complete error:", error);

    return NextResponse.json(
      { error: "Internal server error", message: error.message },
      { status: 500, headers: corsHeaders },
    );
  }
}
