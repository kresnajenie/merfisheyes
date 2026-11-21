import { NextRequest, NextResponse } from "next/server";

import { prisma } from "@/lib/prisma";

// GET /api/presets/versions/[hash] — one version, published or not. Public,
// so a preview link (?preset=<hash>) works for whoever it is sent to.
export async function GET(
  _req: NextRequest,
  { params }: { params: Promise<{ hash: string }> },
) {
  const { hash } = await params;
  const version = await prisma.technologyPresetVersion.findUnique({
    where: { hash },
    select: { hash: true, technology: true, config: true },
  });

  if (!version) {
    return NextResponse.json({ error: "No such preset" }, { status: 404 });
  }

  return NextResponse.json(version);
}
