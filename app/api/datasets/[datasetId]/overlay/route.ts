import { NextRequest, NextResponse } from "next/server";

import { auth } from "@/lib/auth";
import { prisma } from "@/lib/prisma";
import {
  getOverlaySmId,
  writeOverlay,
  removeOverlay,
  groupInProject,
} from "@/lib/ingest/overlay";

// Manage the single-molecule overlay on a single-cell dataset the user owns.
//   GET    → { smDatasetId | null }
//   POST   { smDatasetId, force? } → attach/replace (409 if one exists and !force)
//   DELETE → detach (leaves the project grouping intact)

async function ownedDataset(id: string, userId: string, isAdmin: boolean) {
  const d = await prisma.dataset.findUnique({
    where: { id },
    select: { id: true, ownerId: true, adminOwned: true, datasetType: true, title: true },
  });

  if (!d) return null;
  const owned = d.ownerId === userId || (isAdmin && d.adminOwned);

  return owned ? d : null;
}

export async function GET(
  _req: NextRequest,
  { params }: { params: Promise<{ datasetId: string }> },
) {
  const { datasetId } = await params;

  return NextResponse.json({ smDatasetId: await getOverlaySmId(datasetId) });
}

export async function POST(
  req: NextRequest,
  { params }: { params: Promise<{ datasetId: string }> },
) {
  const session = await auth();

  if (!session?.user?.id) {
    return NextResponse.json({ error: "Not signed in" }, { status: 401 });
  }
  const isAdmin =
    session.user.role === "ADMIN" || session.user.role === "SUPER_ADMIN";
  const { datasetId } = await params;
  const body = await req.json().catch(() => ({}));
  const smDatasetId = typeof body.smDatasetId === "string" ? body.smDatasetId : "";
  const force = body.force === true;

  if (!smDatasetId) {
    return NextResponse.json({ error: "smDatasetId is required" }, { status: 400 });
  }

  const sc = await ownedDataset(datasetId, session.user.id, isAdmin);

  if (!sc) {
    return NextResponse.json({ error: "Dataset not found" }, { status: 404 });
  }
  if (sc.datasetType !== "single_cell") {
    return NextResponse.json(
      { error: "Overlays attach to single-cell datasets" },
      { status: 400 },
    );
  }

  const sm = await ownedDataset(smDatasetId, session.user.id, isAdmin);

  if (!sm || sm.datasetType !== "single_molecule") {
    return NextResponse.json(
      { error: "Pick a single-molecule dataset you own" },
      { status: 400 },
    );
  }

  const existing = await getOverlaySmId(datasetId);

  if (existing && existing !== smDatasetId && !force) {
    return NextResponse.json(
      { error: "overlay_exists", smDatasetId: existing },
      { status: 409 },
    );
  }

  await writeOverlay(datasetId, smDatasetId);
  const projectId = await groupInProject(
    sc.ownerId ?? session.user.id,
    sc.title ?? "Overlay",
    datasetId,
    smDatasetId,
  );

  return NextResponse.json({ ok: true, smDatasetId, projectId });
}

export async function DELETE(
  _req: NextRequest,
  { params }: { params: Promise<{ datasetId: string }> },
) {
  const session = await auth();

  if (!session?.user?.id) {
    return NextResponse.json({ error: "Not signed in" }, { status: 401 });
  }
  const isAdmin =
    session.user.role === "ADMIN" || session.user.role === "SUPER_ADMIN";
  const { datasetId } = await params;
  const sc = await ownedDataset(datasetId, session.user.id, isAdmin);

  if (!sc) {
    return NextResponse.json({ error: "Dataset not found" }, { status: 404 });
  }
  await removeOverlay(datasetId);

  return NextResponse.json({ ok: true });
}
