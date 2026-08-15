import { NextRequest, NextResponse } from "next/server";

import { getObjectBytes } from "@/lib/s3";
import { prisma } from "@/lib/prisma";

// GET — serve the dataset's mapping.json (public; the SC viewer fetches it to
// decide whether to overlay a linked single-molecule dataset).
//
// Source of truth is the stored datasets/{id}/mapping.json (written by the
// ingest finalize step). If that file isn't there yet, we fall back to
// synthesizing the "__all__" config from the dataset's own
// processingParams.linkedSmDatasetId — set when the two were uploaded together.
// That makes the overlay work as soon as both datasets exist, without depending
// on the finalize file write having landed. 404 only when neither is present.
export async function GET(
  _req: NextRequest,
  { params }: { params: Promise<{ datasetId: string }> },
) {
  const { datasetId } = await params;
  const obj = await getObjectBytes(`datasets/${datasetId}/mapping.json`);

  if (obj) {
    return new Response(new Uint8Array(obj.body), {
      headers: {
        "Content-Type": "application/json",
        "Cache-Control": "public, max-age=60",
      },
    });
  }

  // Fallback: derive the mapping from the combined-upload link on the row.
  const dataset = await prisma.dataset.findUnique({
    where: { id: datasetId },
    select: { processingParams: true },
  });
  const linkedSmId = (
    dataset?.processingParams as { linkedSmDatasetId?: unknown } | null
  )?.linkedSmDatasetId;

  if (typeof linkedSmId === "string" && linkedSmId) {
    return NextResponse.json(
      {
        linkColumn: "__all__",
        links: { __all__: linkedSmId },
        source: "app",
      },
      { headers: { "Cache-Control": "public, max-age=60" } },
    );
  }

  return NextResponse.json({ error: "Not found" }, { status: 404 });
}
