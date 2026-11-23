import { NextRequest, NextResponse } from "next/server";

import pkg from "@/package.json";
import { RECOGNISED_TECHNOLOGIES } from "@/lib/config/technology-presets";
import { requirePresetEditor } from "@/lib/presets/permission";
import { presetConfigError, presetHash } from "@/lib/presets/validate";
import { prisma } from "@/lib/prisma";

const who = { select: { name: true, email: true } };

// GET /api/presets/versions — every saved version and the publish log, newest
// first. Preset editors only.
export async function GET() {
  const { error } = await requirePresetEditor();

  if (error) return error;

  const [versions, publishes] = await Promise.all([
    prisma.technologyPresetVersion.findMany({
      orderBy: { createdAt: "desc" },
      include: { createdBy: who },
    }),
    prisma.technologyPresetPublish.findMany({
      orderBy: { createdAt: "desc" },
      include: { publishedBy: who },
    }),
  ]);

  return NextResponse.json({ versions, publishes });
}

// POST /api/presets/versions — save a config as a version (a draft until it
// is published). Saving a config that already exists returns that version.
export async function POST(req: NextRequest) {
  const { error, session } = await requirePresetEditor();

  if (error) return error;

  const body = await req.json().catch(() => null);
  const technology = body?.technology;

  if (!RECOGNISED_TECHNOLOGIES.includes(technology)) {
    return NextResponse.json(
      {
        error: `technology must be one of: ${RECOGNISED_TECHNOLOGIES.join(", ")}`,
      },
      { status: 400 },
    );
  }
  const problem = presetConfigError(body.config);

  if (problem) return NextResponse.json({ error: problem }, { status: 400 });

  const hash = presetHash(technology, body.config);
  const text = (v: unknown) =>
    typeof v === "string" && v.trim() ? v.trim().slice(0, 1000) : null;

  const version = await prisma.technologyPresetVersion.upsert({
    where: { hash },
    update: {},
    create: {
      hash,
      technology,
      config: body.config,
      note: text(body.note),
      testedOn: text(body.testedOn),
      appVersion: pkg.version,
      createdById: session.user.id,
    },
  });

  return NextResponse.json({ version });
}
