import { NextRequest, NextResponse } from "next/server";

import { prisma } from "@/lib/prisma";
import { requireAdmin } from "@/lib/admin-auth";

const includeEntries = { entries: { orderBy: { sortOrder: "asc" as const } } };

// POST /api/admin/community/[id] — approve or reject a community submission.
export async function POST(
  request: NextRequest,
  { params }: { params: Promise<{ id: string }> },
) {
  const { error } = await requireAdmin();

  if (error) return error;

  const { id } = await params;

  let body: Record<string, unknown>;

  try {
    body = await request.json();
  } catch {
    return NextResponse.json({ error: "Invalid JSON" }, { status: 400 });
  }

  const action = body.action;

  if (action !== "approve" && action !== "reject") {
    return NextResponse.json({ error: "Invalid action" }, { status: 400 });
  }

  const existing = await prisma.catalogDataset.findUnique({
    where: { id },
    select: { id: true, isCommunity: true },
  });

  if (!existing || !existing.isCommunity) {
    return NextResponse.json({ error: "Not found" }, { status: 404 });
  }

  const note = typeof body.note === "string" ? body.note.trim() || null : null;

  const item = await prisma.catalogDataset.update({
    where: { id },
    data:
      action === "approve"
        ? { reviewStatus: "APPROVED", isPublished: true, reviewNote: null }
        : { reviewStatus: "REJECTED", isPublished: false, reviewNote: note },
    include: includeEntries,
  });

  return NextResponse.json(item);
}
