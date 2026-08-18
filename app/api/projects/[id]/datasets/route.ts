import { NextRequest, NextResponse } from "next/server";

import { requireUser } from "@/lib/admin-auth";
import { prisma } from "@/lib/prisma";

function isAdmin(session: { user: { role?: string | null } }): boolean {
  return session.user.role === "ADMIN" || session.user.role === "SUPER_ADMIN";
}

// POST /api/projects/[id]/datasets — add a dataset to the project.
export async function POST(
  request: NextRequest,
  { params }: { params: Promise<{ id: string }> },
) {
  const { error, session } = await requireUser();

  if (error) return error;

  const { id } = await params;
  const project = await prisma.project.findUnique({
    where: { id },
    select: { id: true, ownerId: true },
  });

  if (!project) {
    return NextResponse.json({ error: "Not found" }, { status: 404 });
  }

  if (project.ownerId !== session.user.id && !isAdmin(session)) {
    return NextResponse.json({ error: "Forbidden" }, { status: 403 });
  }

  let body: Record<string, unknown>;

  try {
    body = await request.json();
  } catch {
    return NextResponse.json({ error: "Invalid JSON" }, { status: 400 });
  }

  const datasetId = typeof body.datasetId === "string" ? body.datasetId : "";

  if (!datasetId) {
    return NextResponse.json(
      { error: "datasetId is required" },
      { status: 400 },
    );
  }

  const dataset = await prisma.dataset.findUnique({
    where: { id: datasetId },
    select: { id: true, ownerId: true, adminOwned: true },
  });

  if (!dataset) {
    return NextResponse.json({ error: "Dataset not found" }, { status: 400 });
  }

  // Only datasets belonging to the project owner (or admin-owned shared ones)
  // may be added — you can't file someone else's dataset into your project.
  const belongsToOwner =
    Boolean(dataset.ownerId) && dataset.ownerId === project.ownerId;

  if (!belongsToOwner && !dataset.adminOwned) {
    return NextResponse.json(
      { error: "Cannot add a dataset you don't own" },
      { status: 403 },
    );
  }

  // Idempotent: already a member → return ok.
  const existing = await prisma.projectDataset.findUnique({
    where: { projectId_datasetId: { projectId: id, datasetId } },
    select: { projectId: true },
  });

  if (existing) {
    return NextResponse.json({ ok: true });
  }

  const last = await prisma.projectDataset.findFirst({
    where: { projectId: id },
    orderBy: { sortOrder: "desc" },
    select: { sortOrder: true },
  });
  const sortOrder = (last?.sortOrder ?? -1) + 1;

  await prisma.projectDataset.create({
    data: { projectId: id, datasetId, sortOrder },
  });

  return NextResponse.json({ ok: true });
}
