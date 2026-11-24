import { NextRequest, NextResponse } from "next/server";

import { RECOGNISED_TECHNOLOGIES } from "@/lib/config/technology-presets";
import { requirePresetEditor } from "@/lib/presets/permission";
import { prisma } from "@/lib/prisma";

// POST /api/presets/publish — put a saved version in force for its
// technology ({ technology, hash }), or turn the technology's preset off
// ({ technology, hash: null }). Rolling back is publishing an earlier version.
export async function POST(req: NextRequest) {
  const { error, session } = await requirePresetEditor();

  if (error) return error;

  const body = await req.json().catch(() => null);
  const technology = body?.technology;
  const hash = body?.hash ?? null;

  if (!RECOGNISED_TECHNOLOGIES.includes(technology)) {
    return NextResponse.json({ error: "Unknown technology" }, { status: 400 });
  }
  if (hash !== null) {
    const version = await prisma.technologyPresetVersion.findUnique({
      where: { hash: String(hash) },
      select: { technology: true },
    });

    if (version?.technology !== technology) {
      return NextResponse.json(
        { error: "No such version of this technology" },
        { status: 404 },
      );
    }
  }

  const publish = await prisma.technologyPresetPublish.create({
    data: { technology, versionHash: hash, publishedById: session.user.id },
  });

  return NextResponse.json({ publish });
}
