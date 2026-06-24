// app/api/datasets/[datasetId]/list/route.ts
// Lists every S3 key under datasets/{datasetId}/. Used by the zarr read path
// to enumerate the store (zarrita's AsyncReadable has no LIST op).
import { ListObjectsV2Command } from "@aws-sdk/client-s3";
import { NextRequest, NextResponse } from "next/server";

import { S3_BUCKET, s3Client } from "@/lib/s3";

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
  { params }: { params: Promise<{ datasetId: string }> },
) {
  try {
    const { datasetId } = await params;

    if (!datasetId.startsWith("ds_")) {
      return NextResponse.json(
        { error: "Invalid dataset ID format" },
        { status: 400, headers: corsHeaders },
      );
    }
    if (!S3_BUCKET) {
      return NextResponse.json(
        { error: "AWS_S3_BUCKET not configured" },
        { status: 500, headers: corsHeaders },
      );
    }

    const prefix = `datasets/${datasetId}/`;
    const keys: string[] = [];
    let continuationToken: string | undefined = undefined;

    // Paginate — S3 returns up to 1000 keys per page; zarr stores can exceed.
    do {
      const cmd: ListObjectsV2Command = new ListObjectsV2Command({
        Bucket: S3_BUCKET,
        Prefix: prefix,
        ContinuationToken: continuationToken,
      });
      const res: any = await s3Client.send(cmd);

      for (const obj of res.Contents || []) {
        if (!obj.Key) continue;
        // Strip prefix; client works with paths relative to dataset root
        keys.push(obj.Key.slice(prefix.length));
      }
      continuationToken = res.IsTruncated ? res.NextContinuationToken : undefined;
    } while (continuationToken);

    return NextResponse.json({ keys }, { headers: corsHeaders });
  } catch (error: any) {
    console.error("/list error:", error);

    return NextResponse.json(
      { error: "Internal server error", message: error.message },
      { status: 500, headers: corsHeaders },
    );
  }
}
