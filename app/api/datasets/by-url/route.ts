import { NextRequest, NextResponse } from "next/server";

import { prisma } from "@/lib/prisma";
import { normalizeS3Url } from "@/lib/utils/viewer-config";

const corsHeaders = {
  "Access-Control-Allow-Origin": process.env.CORS_ORIGIN || "*",
  "Access-Control-Allow-Methods": "GET, OPTIONS",
  "Access-Control-Allow-Headers": "Content-Type",
};

export async function OPTIONS() {
  return NextResponse.json({}, { headers: corsHeaders });
}

/**
 * Look up a dataset registered from a raw S3 URL. Public — the viewer uses it
 * to decide whether to apply owner-saved defaults, whether to show the
 * "claim" banner, and which id to count views against. 404 = unregistered.
 */
export async function GET(request: NextRequest) {
  const url = request.nextUrl.searchParams.get("url");

  if (!url) {
    return NextResponse.json(
      { error: "Missing url" },
      { status: 400, headers: corsHeaders },
    );
  }

  const s3BaseUrl = normalizeS3Url(url);
  const dataset = await prisma.dataset.findUnique({
    where: { s3BaseUrl },
    select: {
      id: true,
      ownerId: true,
      adminOwned: true,
      viewerConfig: true,
      viewCount: true,
      title: true,
      // First project this dataset belongs to, so a viewer can show the rail
      // of its siblings without a second round trip.
      projects: {
        orderBy: { project: { createdAt: "asc" } },
        take: 1,
        select: { projectId: true },
      },
    },
  });

  if (!dataset) {
    return NextResponse.json(
      { registered: false },
      { status: 404, headers: corsHeaders },
    );
  }

  const { projects, ...rest } = dataset;

  return NextResponse.json(
    { registered: true, ...rest, projectId: projects[0]?.projectId ?? null },
    { headers: corsHeaders },
  );
}
