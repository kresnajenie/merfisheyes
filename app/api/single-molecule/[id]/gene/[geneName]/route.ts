import { NextRequest, NextResponse } from "next/server";

import { generatePresignedDownloadUrl } from "@/lib/s3";

const corsHeaders = {
  "Access-Control-Allow-Origin": process.env.CORS_ORIGIN || "*",
  "Access-Control-Allow-Methods": "GET, OPTIONS",
  "Access-Control-Allow-Headers": "Content-Type",
};

export async function OPTIONS() {
  return NextResponse.json({}, { headers: corsHeaders });
}

/**
 * Generate presigned URL for a specific gene file
 * Returns URL to download the gene's .bin.gz file from S3
 */
export async function GET(
  request: NextRequest,
  { params }: { params: Promise<{ id: string; geneName: string }> },
) {
  try {
    const { id: datasetId, geneName } = await params;

    // Validate dataset ID format
    if (!datasetId || !datasetId.startsWith("sm_")) {
      return NextResponse.json(
        { error: "Invalid single molecule dataset ID format" },
        { status: 400, headers: corsHeaders },
      );
    }

    if (!geneName) {
      return NextResponse.json(
        { error: "Gene name is required" },
        { status: 400, headers: corsHeaders },
      );
    }

    // Sanitize gene name to match filename (same logic as processor)
    const sanitizedName = geneName
      .replace(/[^a-zA-Z0-9]/g, "_")
      .replace(/_+/g, "_")
      .replace(/^_|_$/g, "");

    // Generate S3 key
    const geneKey = `datasets/${datasetId}/genes/${sanitizedName}.bin.gz`;

    // Generate presigned URL (valid for 1 hour)
    const presignedUrl = await generatePresignedDownloadUrl(geneKey, 3600);

    return NextResponse.json(
      {
        geneName,
        sanitizedName,
        url: presignedUrl,
        expiresIn: 3600,
      },
      { headers: corsHeaders },
    );
  } catch (error: any) {
    console.error("Generate gene URL error:", error);

    return NextResponse.json(
      {
        error: "Internal server error",
        message: error.message,
      },
      { status: 500, headers: corsHeaders },
    );
  }
}
