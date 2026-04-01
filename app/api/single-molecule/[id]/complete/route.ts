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

export async function POST(
  request: NextRequest,
  { params }: { params: Promise<{ id: string }> },
) {
  try {
    const { id: datasetId } = await params;
    const { uploadId } = await request.json();

    if (!uploadId) {
      return NextResponse.json(
        { error: "uploadId is required" },
        { status: 400, headers: corsHeaders },
      );
    }

    // Verify all files are complete
    const uploadSession = await prisma.uploadSession.findUnique({
      where: { id: uploadId },
      include: {
        files: true,
        dataset: true,
      },
    });

    if (!uploadSession) {
      return NextResponse.json(
        { error: "Upload session not found" },
        { status: 404, headers: corsHeaders },
      );
    }

    if (uploadSession.completedFiles !== uploadSession.totalFiles) {
      return NextResponse.json(
        {
          error: "Upload incomplete",
          message: `Only ${uploadSession.completedFiles} of ${uploadSession.totalFiles} files uploaded`,
        },
        { status: 400, headers: corsHeaders },
      );
    }

    // Mark dataset as complete
    await prisma.dataset.update({
      where: { id: datasetId },
      data: {
        status: "COMPLETE",
        completedAt: new Date(),
      },
    });

    // Email notification is sent from the client side to avoid
    // Vercel deployment protection blocking server-to-server requests

    return NextResponse.json(
      {
        success: true,
        message: "Upload completed successfully",
        datasetId,
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
