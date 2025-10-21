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
  { params }: { params: Promise<{ id: string; key: string }> }
) {
  try {
    const { id: datasetId, key: encodedFileKey } = await params;
    const fileKey = decodeURIComponent(encodedFileKey);
    const { uploadId } = await request.json();

    if (!uploadId) {
      return NextResponse.json(
        { error: "uploadId is required" },
        { status: 400, headers: corsHeaders }
      );
    }

    // Update file status and increment completed files count
    await prisma.$transaction(async (tx) => {
      // 1. Mark file as complete
      await tx.uploadFile.updateMany({
        where: {
          uploadSessionId: uploadId,
          fileKey: fileKey,
        },
        data: {
          status: "COMPLETE",
          uploadedAt: new Date(),
        },
      });

      // 2. Increment completed files count
      await tx.uploadSession.update({
        where: { id: uploadId },
        data: {
          completedFiles: {
            increment: 1,
          },
        },
      });
    });

    return NextResponse.json(
      { success: true, message: "File marked as complete" },
      { headers: corsHeaders }
    );
  } catch (error: any) {
    console.error("Mark file complete error:", error);
    return NextResponse.json(
      {
        error: "Internal server error",
        message: error.message,
      },
      { status: 500, headers: corsHeaders }
    );
  }
}
