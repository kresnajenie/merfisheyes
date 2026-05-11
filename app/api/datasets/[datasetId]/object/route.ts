// app/api/datasets/[datasetId]/object/route.ts
// Returns a presigned GET URL for a single object under datasets/{id}/.
// The zarr read path calls this on demand for each chunk it reads.
import { NextRequest, NextResponse } from "next/server";

import { generatePresignedDownloadUrl } from "@/lib/s3";

const corsHeaders = {
  "Access-Control-Allow-Origin": process.env.CORS_ORIGIN || "*",
  "Access-Control-Allow-Methods": "GET, OPTIONS",
  "Access-Control-Allow-Headers": "Content-Type",
};
const URL_TTL_SECONDS = 3600;

export async function OPTIONS() {
  return NextResponse.json({}, { headers: corsHeaders });
}

export async function GET(
  request: NextRequest,
  { params }: { params: Promise<{ datasetId: string }> },
) {
  try {
    const { datasetId } = await params;
    const relativeKey = request.nextUrl.searchParams.get("key");

    if (!datasetId.startsWith("ds_")) {
      return NextResponse.json(
        { error: "Invalid dataset ID format" },
        { status: 400, headers: corsHeaders },
      );
    }
    if (!relativeKey || relativeKey.includes("..")) {
      return NextResponse.json(
        { error: "Invalid or missing key parameter" },
        { status: 400, headers: corsHeaders },
      );
    }

    const fullKey = `datasets/${datasetId}/${relativeKey}`;
    const url = await generatePresignedDownloadUrl(fullKey, URL_TTL_SECONDS);

    return NextResponse.json(
      { url, expiresIn: URL_TTL_SECONDS },
      { headers: corsHeaders },
    );
  } catch (error: any) {
    console.error("/object error:", error);

    return NextResponse.json(
      { error: "Internal server error", message: error.message },
      { status: 500, headers: corsHeaders },
    );
  }
}
