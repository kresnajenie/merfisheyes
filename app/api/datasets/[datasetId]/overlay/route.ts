import { NextRequest, NextResponse } from "next/server";

import { auth } from "@/lib/auth";
import { prisma } from "@/lib/prisma";
import {
  getOverlaySmId,
  getOverlayImage,
  writeOverlay,
  writeOverlayImage,
  removeOverlay,
  groupInProject,
} from "@/lib/ingest/overlay";

// Manage the overlays on a single-cell dataset the user owns: a
// single-molecule dataset drawn over the cells, and an image drawn under them.
//   GET    → { smDatasetId | null, image | null }
//   POST   { smDatasetId, force? } → attach/replace (409 if one exists and !force)
//   PATCH  { image: url | null } → set or clear the image overlay
//   DELETE → detach the molecules (leaves the project grouping and image intact)

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

  return NextResponse.json({
    smDatasetId: await getOverlaySmId(datasetId),
    image: (await getOverlayImage(datasetId)) ?? null,
  });
}

export async function PATCH(
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
  let image: string | null = null;

  if (body.image !== null) {
    const raw = typeof body.image === "string" ? body.image.trim() : "";
    const url = URL.canParse(raw) ? new URL(raw) : null;

    if (!url || !/^https?:$/.test(url.protocol)) {
      return NextResponse.json(
        { error: "image must be an http(s) URL, or null to remove it" },
        { status: 400 },
      );
    }
    image = url.href.replace(/\/+$/, "");
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
  await writeOverlayImage(datasetId, image);

  return NextResponse.json({ ok: true, image });
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
