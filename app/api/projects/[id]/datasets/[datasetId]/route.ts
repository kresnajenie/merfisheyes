import { NextRequest, NextResponse } from "next/server";

import { requireUser } from "@/lib/admin-auth";
import { prisma } from "@/lib/prisma";

function isAdmin(session: { user: { role?: string | null } }): boolean {
  return session.user.role === "ADMIN" || session.user.role === "SUPER_ADMIN";
}

// DELETE /api/projects/[id]/datasets/[datasetId] — remove a dataset membership.
export async function DELETE(
  _req: NextRequest,
  { params }: { params: Promise<{ id: string; datasetId: string }> },
) {
  const { error, session } = await requireUser();

  if (error) return error;

  const { id, datasetId } = await params;
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

  await prisma.projectDataset.deleteMany({
    where: { projectId: id, datasetId },
  });

  return NextResponse.json({ ok: true });
}
