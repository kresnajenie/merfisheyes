import { NextRequest, NextResponse } from "next/server";

import { requireUser } from "@/lib/admin-auth";
import { prisma } from "@/lib/prisma";
import {
  sanitizeMetadataPatch,
  sanitizeMetadataBag,
} from "@/lib/datasets/metadata";

function isAdmin(session: { user: { role?: string | null } }): boolean {
  return session.user.role === "ADMIN" || session.user.role === "SUPER_ADMIN";
}

const includeDatasets = {
  datasets: {
    orderBy: { sortOrder: "asc" as const },
    include: {
      dataset: {
        select: {
          id: true,
          title: true,
          datasetType: true,
          thumbnailUrl: true,
          status: true,
          numCells: true,
          numGenes: true,
        },
      },
    },
  },
};

// GET /api/projects/[id]
export async function GET(
  _req: NextRequest,
  { params }: { params: Promise<{ id: string }> },
) {
  const { error, session } = await requireUser();

  if (error) return error;

  const { id } = await params;
  const project = await prisma.project.findUnique({
    where: { id },
    include: includeDatasets,
  });

  if (!project) {
    return NextResponse.json({ error: "Not found" }, { status: 404 });
  }

  if (project.ownerId !== session.user.id && !isAdmin(session)) {
    return NextResponse.json({ error: "Forbidden" }, { status: 403 });
  }

  const submission = await prisma.catalogDataset.findFirst({
    where: { isCommunity: true, sourceProjectId: project.id },
    select: {
      id: true,
      reviewStatus: true,
      isPublished: true,
      reviewNote: true,
    },
  });

  return NextResponse.json({
    project: {
      ...project,
      datasetCount: project.datasets.length,
      datasets: project.datasets.map((pd) => ({
        ...pd.dataset,
        sortOrder: pd.sortOrder,
      })),
    },
    submission: submission
      ? {
          catalogId: submission.id,
          reviewStatus: submission.reviewStatus,
          isPublished: submission.isPublished,
          reviewNote: submission.reviewNote,
        }
      : null,
  });
}

// PATCH /api/projects/[id]
export async function PATCH(
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

  const data: Record<string, unknown> = { ...sanitizeMetadataPatch(body) };

  if (typeof body.title === "string") {
    const trimmed = body.title.trim();

    if (trimmed) data.title = trimmed.slice(0, 500);
  }

  if (typeof body.thumbnailUrl === "string" || body.thumbnailUrl === null) {
    data.thumbnailUrl =
      typeof body.thumbnailUrl === "string"
        ? body.thumbnailUrl.trim().slice(0, 1000) || null
        : null;
  }

  if ("metadata" in body) {
    data.metadata = sanitizeMetadataBag(body.metadata);
  }

  if (Object.keys(data).length === 0) {
    return NextResponse.json({ error: "Nothing to update" }, { status: 400 });
  }

  const updated = await prisma.project.update({ where: { id }, data });

  return NextResponse.json(updated);
}

// DELETE /api/projects/[id]
export async function DELETE(
  _req: NextRequest,
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

  await prisma.project.delete({ where: { id } });

  return NextResponse.json({ ok: true });
}
