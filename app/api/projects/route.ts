import { NextRequest, NextResponse } from "next/server";

import { requireUser } from "@/lib/admin-auth";
import { prisma } from "@/lib/prisma";
import { sanitizeMetadataPatch } from "@/lib/datasets/metadata";

// GET /api/projects — the caller's own projects, with member datasets and any
// community-submission status.
export async function GET() {
  const { error, session } = await requireUser();

  if (error) return error;

  const projects = await prisma.project.findMany({
    where: { ownerId: session.user.id },
    orderBy: { updatedAt: "desc" },
    include: {
      datasets: {
        orderBy: { sortOrder: "asc" },
        include: {
          dataset: {
            select: {
              id: true,
              title: true,
              datasetType: true,
              thumbnailUrl: true,
              status: true,
            },
          },
        },
      },
    },
  });

  const submissions = await prisma.catalogDataset.findMany({
    where: {
      isCommunity: true,
      sourceProjectId: { in: projects.map((p) => p.id) },
    },
    select: {
      id: true,
      sourceProjectId: true,
      reviewStatus: true,
      isPublished: true,
      reviewNote: true,
    },
  });
  const submissionByProject = new Map(
    submissions.map((s) => [s.sourceProjectId, s]),
  );

  return NextResponse.json({
    projects: projects.map((p) => {
      const sub = submissionByProject.get(p.id);

      return {
        id: p.id,
        title: p.title,
        description: p.description,
        species: p.species,
        disease: p.disease,
        tissue: p.tissue,
        platform: p.platform,
        institute: p.institute,
        tags: p.tags,
        thumbnailUrl: p.thumbnailUrl,
        externalLink: p.externalLink,
        publicationLink: p.publicationLink,
        createdAt: p.createdAt,
        updatedAt: p.updatedAt,
        datasetCount: p.datasets.length,
        datasets: p.datasets.map((pd) => pd.dataset),
        submission: sub
          ? {
              catalogId: sub.id,
              reviewStatus: sub.reviewStatus,
              isPublished: sub.isPublished,
              reviewNote: sub.reviewNote,
            }
          : null,
      };
    }),
  });
}

// POST /api/projects — create a project owned by the caller.
export async function POST(request: NextRequest) {
  const { error, session } = await requireUser();

  if (error) return error;

  let body: Record<string, unknown>;

  try {
    body = await request.json();
  } catch {
    return NextResponse.json({ error: "Invalid JSON" }, { status: 400 });
  }

  const title = typeof body.title === "string" ? body.title.trim() : "";

  if (!title) {
    return NextResponse.json({ error: "Title is required" }, { status: 400 });
  }

  const thumbnailUrl =
    typeof body.thumbnailUrl === "string"
      ? body.thumbnailUrl.trim().slice(0, 1000) || null
      : undefined;

  const project = await prisma.project.create({
    data: {
      title: title.slice(0, 500),
      ownerId: session.user.id,
      ...sanitizeMetadataPatch(body),
      ...(thumbnailUrl !== undefined ? { thumbnailUrl } : {}),
    },
  });

  return NextResponse.json(project, { status: 201 });
}
