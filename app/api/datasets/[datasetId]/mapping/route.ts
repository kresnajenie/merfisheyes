import { NextRequest, NextResponse } from "next/server";

import { getObjectBytes } from "@/lib/s3";

// GET — serve the dataset's mapping.json (public; the SC viewer fetches it to
// decide whether to overlay a linked single-molecule dataset). 404 when absent.
export async function GET(
  _req: NextRequest,
  { params }: { params: Promise<{ datasetId: string }> },
) {
  const { datasetId } = await params;
  const obj = await getObjectBytes(`datasets/${datasetId}/mapping.json`);

  if (!obj) {
    return NextResponse.json({ error: "Not found" }, { status: 404 });
  }

  return new Response(new Uint8Array(obj.body), {
    headers: {
      "Content-Type": "application/json",
      "Cache-Control": "public, max-age=60",
    },
  });
}
