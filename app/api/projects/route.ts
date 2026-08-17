import { NextRequest, NextResponse } from "next/server";

import { requireUser } from "@/lib/admin-auth";
import { prisma } from "@/lib/prisma";
import {
  sanitizeMetadataPatch,
  sanitizeMetadataBag,
} from "@/lib/datasets/metadata";

// GET /api/projects — the caller's own projects, with member datasets and any
// community-submission status. Supports optional search + pagination (opt-in:
// callers with no page/limit get the full list, unchanged).
export async function GET(request: NextRequest) {
  const { error, session } = await requireUser();

  if (error) return error;

  const params = request.nextUrl.searchParams;
  const search = params.get("search")?.trim() ?? "";
  const paginate = params.has("page") || params.has("limit");
  const page = Math.max(1, Number(params.get("page") ?? "1"));
  const limit = paginate
    ? Math.min(100, Math.max(1, Number(params.get("limit") ?? "20")))
    : undefined;
  const skip = limit ? (page - 1) * limit : undefined;

  const insensitive = "insensitive" as const;
  const where = {
    ownerId: session.user.id,
    ...(search
      ? {
          OR: [
            { title: { contains: search, mode: insensitive } },
            { description: { contains: search, mode: insensitive } },
            { tags: { has: search } },
          ],
        }
      : {}),
  };

  const [projects, total] = await Promise.all([
    prisma.project.findMany({
      where,
      orderBy: { updatedAt: "desc" },
      skip,
      take: limit,
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
    }),
    prisma.project.count({ where }),
  ]);

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
    total,
    page,
    limit: limit ?? total,
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
        metadata: p.metadata,
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
      ...("metadata" in body
        ? { metadata: sanitizeMetadataBag(body.metadata) }
        : {}),
    },
  });

  return NextResponse.json(project, { status: 201 });
}
