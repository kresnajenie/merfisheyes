import { NextRequest, NextResponse } from "next/server";

import { requireUser, ownerOrAdminError } from "@/lib/admin-auth";
import { prisma } from "@/lib/prisma";
import { upsertCommunitySubmission } from "@/lib/community/submit";

function isAdmin(session: { user: { role?: string | null } }): boolean {
  return session.user.role === "ADMIN" || session.user.role === "SUPER_ADMIN";
}

const datasetSelect = {
  id: true,
  title: true,
  datasetType: true,
  ingestSource: true,
  s3BaseUrl: true,
  thumbnailUrl: true,
  description: true,
  species: true,
  disease: true,
  tissue: true,
  platform: true,
  institute: true,
  tags: true,
  externalLink: true,
  publicationLink: true,
  metadata: true,
  numCells: true,
  numGenes: true,
  status: true,
  ownerId: true,
  adminOwned: true,
} as const;

// POST /api/community/submit — submit an owned dataset or project for review.
export async function POST(request: NextRequest) {
  const { error, session } = await requireUser();

  if (error) return error;

  let body: Record<string, unknown>;

  try {
    body = await request.json();
  } catch {
    return NextResponse.json({ error: "Invalid JSON" }, { status: 400 });
  }

  const type = body.type;
  const id = typeof body.id === "string" ? body.id : "";

  if (!id || (type !== "dataset" && type !== "project")) {
    return NextResponse.json({ error: "Invalid request" }, { status: 400 });
  }

  if (type === "dataset") {
    const dataset = await prisma.dataset.findUnique({
      where: { id },
      select: datasetSelect,
    });

    if (!dataset) {
      return NextResponse.json({ error: "Not found" }, { status: 404 });
    }

    const forbidden = ownerOrAdminError(dataset, session);

    if (forbidden) return forbidden;

    if (dataset.status !== "COMPLETE") {
      return NextResponse.json({ error: "Dataset not ready" }, { status: 400 });
    }

    const catalog = await upsertCommunitySubmission({
      userId: session.user.id,
      source: { kind: "dataset", dataset },
    });

    return NextResponse.json({
      ok: true,
      catalogId: catalog.id,
      reviewStatus: catalog.reviewStatus,
    });
  }

  // type === "project"
  const project = await prisma.project.findUnique({
    where: { id },
    include: {
      datasets: {
        orderBy: { sortOrder: "asc" },
        include: { dataset: { select: datasetSelect } },
      },
    },
  });

  if (!project) {
    return NextResponse.json({ error: "Not found" }, { status: 404 });
  }

  if (project.ownerId !== session.user.id && !isAdmin(session)) {
    return NextResponse.json({ error: "Forbidden" }, { status: 403 });
  }

  const datasets = project.datasets.map((pd) => pd.dataset);

  if (datasets.length === 0) {
    return NextResponse.json(
      { error: "Project has no datasets" },
      { status: 400 },
    );
  }

  const catalog = await upsertCommunitySubmission({
    userId: session.user.id,
    source: {
      kind: "project",
      project: {
        id: project.id,
        title: project.title,
        thumbnailUrl: project.thumbnailUrl,
        description: project.description,
        species: project.species,
        disease: project.disease,
        tissue: project.tissue,
        platform: project.platform,
        institute: project.institute,
        tags: project.tags,
        externalLink: project.externalLink,
        publicationLink: project.publicationLink,
        metadata: project.metadata,
        datasets,
      },
    },
  });

  return NextResponse.json({
    ok: true,
    catalogId: catalog.id,
    reviewStatus: catalog.reviewStatus,
  });
}

// DELETE /api/community/submit — withdraw a pending/live community submission.
export async function DELETE(request: NextRequest) {
  const { error, session } = await requireUser();

  if (error) return error;

  let type: unknown;
  let id = "";

  try {
    const body = await request.json();

    type = body.type;
    id = typeof body.id === "string" ? body.id : "";
  } catch {
    // Fall back to query params.
    type = request.nextUrl.searchParams.get("type");
    id = request.nextUrl.searchParams.get("id") ?? "";
  }

  if (!id || (type !== "dataset" && type !== "project")) {
    return NextResponse.json({ error: "Invalid request" }, { status: 400 });
  }

  const submission = await prisma.catalogDataset.findFirst({
    where: {
      isCommunity: true,
      ...(type === "dataset"
        ? { sourceDatasetId: id }
        : { sourceProjectId: id }),
    },
    select: { id: true, submittedById: true },
  });

  if (!submission) {
    return NextResponse.json({ error: "Not found" }, { status: 404 });
  }

  if (submission.submittedById !== session.user.id && !isAdmin(session)) {
    return NextResponse.json({ error: "Forbidden" }, { status: 403 });
  }

  await prisma.catalogDataset.delete({ where: { id: submission.id } });

  return NextResponse.json({ ok: true });
}
