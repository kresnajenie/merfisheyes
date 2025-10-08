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
  { params }: { params: Promise<{ datasetId: string; fileKey: string }> }
) {
  try {
    const { datasetId, fileKey } = await params;
    const body = await request.json();
    const { uploadId } = body;

    if (!uploadId) {
      return NextResponse.json(
        { error: "uploadId is required" },
        { status: 400, headers: corsHeaders }
      );
    }

    // Decode fileKey (it comes URL-encoded)
    const decodedFileKey = decodeURIComponent(fileKey);

    // Update file status to COMPLETE
    const file = await prisma.uploadFile.updateMany({
      where: {
        uploadSessionId: uploadId,
        fileKey: decodedFileKey,
      },
      data: {
        status: "COMPLETE",
      },
    });

    if (file.count === 0) {
      return NextResponse.json(
        { error: "File not found" },
        { status: 404, headers: corsHeaders }
      );
    }

    return NextResponse.json(
      {
        success: true,
        fileKey: decodedFileKey,
        status: "COMPLETE",
      },
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
