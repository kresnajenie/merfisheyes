import { NextRequest, NextResponse } from "next/server";

import { auth } from "@/lib/auth";
import { datasetKind } from "@/lib/datasets/kind";
import {
  getOverlaySmId,
  getOverlayUrl,
  getOverlayImage,
  hasMolecules,
  ownedDataset,
  writeOverlay,
  writeOverlayUrl,
  writeOverlayImage,
  imageOverlayFor,
  type ImageOverlay,
  removeOverlay,
  groupInProject,
} from "@/lib/ingest/overlay";
import { dataUrlOf, moleculeUrlFormat } from "@/lib/overlays/moleculeSource";

// Manage the overlays on a single-cell dataset the user owns: a
// single-molecule dataset drawn over the cells, and an image drawn under them.
//   GET    → { smDatasetId | null, smUrl | null, image | null }
//   POST   { smDatasetId, force? } | { url, force? } → attach/replace (409 if
//          one exists and !force); a url is a public store, a folder holding
//          one, or a molecule folder, overlaid by its link alone
//   PATCH  { imageDatasetId } | { image: url | null } → set or clear the image overlay
//   DELETE → detach the molecules (leaves the project grouping and image intact)

export async function GET(
  _req: NextRequest,
  { params }: { params: Promise<{ datasetId: string }> },
) {
  const { datasetId } = await params;

  const image = await getOverlayImage(datasetId);

  return NextResponse.json({
    smDatasetId: await getOverlaySmId(datasetId),
    smUrl: await getOverlayUrl(datasetId),
    image: image?.image ?? null,
    imageDatasetId: image?.imageDatasetId ?? null,
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
  let image: ImageOverlay | null = null;

  if (typeof body.imageDatasetId === "string" && body.imageDatasetId) {
    // One of the user's own image datasets
    const owned = await ownedDataset(
      body.imageDatasetId,
      session.user.id,
      isAdmin,
    );

    image = owned && (await imageOverlayFor(owned.id));
    if (!image) {
      return NextResponse.json(
        { error: "Pick an image you own" },
        { status: 400 },
      );
    }
  } else if (body.image !== null) {
    const raw = typeof body.image === "string" ? body.image.trim() : "";
    const url = URL.canParse(raw) ? new URL(raw) : null;

    if (!url || !/^https?:$/.test(url.protocol)) {
      return NextResponse.json(
        {
          error:
            "Give imageDatasetId, an http(s) image URL, or image: null to remove it",
        },
        { status: 400 },
      );
    }
    image = { image: url.href.replace(/\/+$/, "") };
  }

  const sc = await ownedDataset(datasetId, session.user.id, isAdmin);

  if (!sc) {
    return NextResponse.json({ error: "Dataset not found" }, { status: 404 });
  }
  if (datasetKind(sc.datasetType) !== "single_cell") {
    return NextResponse.json(
      { error: "Overlays attach to single-cell datasets" },
      { status: 400 },
    );
  }
  await writeOverlayImage(datasetId, image);

  return NextResponse.json({
    ok: true,
    image: image?.image ?? null,
    imageDatasetId: image?.imageDatasetId ?? null,
  });
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
  const smDatasetId =
    typeof body.smDatasetId === "string" ? body.smDatasetId : "";
  const pasted = typeof body.url === "string" ? dataUrlOf(body.url) : null;
  const force = body.force === true;

  if (!smDatasetId && !pasted) {
    return NextResponse.json(
      { error: "Give smDatasetId, or the url of public molecules" },
      { status: 400 },
    );
  }

  const sc = await ownedDataset(datasetId, session.user.id, isAdmin);

  if (!sc) {
    return NextResponse.json({ error: "Dataset not found" }, { status: 404 });
  }
  if (datasetKind(sc.datasetType) !== "single_cell") {
    return NextResponse.json(
      { error: "Overlays attach to single-cell datasets" },
      { status: 400 },
    );
  }

  // What is overlaid now, by dataset or by link: replaced only when asked
  const existing =
    (await getOverlaySmId(datasetId)) ?? (await getOverlayUrl(datasetId));

  // Public molecules by their link alone: checked where they live, the way
  // the viewer will open them, so a wrong or private link fails here
  if (pasted && !smDatasetId) {
    const linkFormat = await moleculeUrlFormat(pasted);

    if (!linkFormat) {
      return NextResponse.json(
        {
          error:
            "No molecules at that link. Give a public SpatialData store with a points element (or the folder holding it), or a merfisheyes molecule folder, that allows this site in CORS.",
        },
        { status: 400 },
      );
    }
    if (existing && existing !== pasted && !force) {
      return NextResponse.json(
        { error: "overlay_exists", smDatasetId: existing },
        { status: 409 },
      );
    }
    await writeOverlayUrl(datasetId, pasted, linkFormat);

    return NextResponse.json({ ok: true, url: pasted, linkFormat });
  }

  const sm = await ownedDataset(smDatasetId, session.user.id, isAdmin);

  if (!sm || smDatasetId === datasetId || !(await hasMolecules(sm))) {
    return NextResponse.json(
      { error: "Pick a dataset of yours with molecules" },
      { status: 400 },
    );
  }

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
