import { NextRequest, NextResponse } from "next/server";
import { DatasetFaultCategory } from "@prisma/client";

import { requireAdmin } from "@/lib/admin-auth";
import { prisma } from "@/lib/prisma";

/**
 * PATCH /api/admin/processing/[id] — an admin overrides the fault tag on a
 * failed dataset (auto-classified on failure; humans get the final say).
 * Body: { faultCategory: "USER" | "PLATFORM" | "UNKNOWN" | null }.
 */
export async function PATCH(
  request: NextRequest,
  { params }: { params: Promise<{ id: string }> },
) {
  const { error } = await requireAdmin();

  if (error) return error;

  const { id } = await params;

  let body: { faultCategory?: unknown };

  try {
    body = await request.json();
  } catch {
    return NextResponse.json({ error: "Invalid JSON" }, { status: 400 });
  }

  const raw = body.faultCategory;
  const valid =
    raw === null ||
    (typeof raw === "string" && raw in DatasetFaultCategory);

  if (!valid) {
    return NextResponse.json(
      { error: "faultCategory must be USER, PLATFORM, UNKNOWN, or null" },
      { status: 400 },
    );
  }

  try {
    await prisma.dataset.update({
      where: { id },
      data: { faultCategory: (raw as DatasetFaultCategory) ?? null },
    });
  } catch {
    return NextResponse.json({ error: "Dataset not found" }, { status: 404 });
  }

  return NextResponse.json({ ok: true, faultCategory: raw });
}
