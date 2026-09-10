import { NextRequest, NextResponse } from "next/server";

import { prisma } from "@/lib/prisma";

/**
 * GET /api/projects/[id]/public — a project's datasets, without a session.
 *
 * The viewer's stage rail needs a project's members while signed out, but
 * `Project` has no public flag, so serving any project by id would expose a
 * user's private one to anyone who has the id.
 *
 * The rule instead: a project is publicly readable only when **every** dataset
 * in it is `adminOwned` — the flag that already means "shared, curated
 * content" rather than someone's personal upload. A project holding even one
 * personal dataset is not public, so this cannot leak private work.
 *
 * If projects later need to be published independently of their contents, that
 * wants an explicit `Project.isPublic` column rather than widening this.
 */
export async function GET(
  _req: NextRequest,
  { params }: { params: Promise<{ id: string }> },
) {
  const { id } = await params;

  const project = await prisma.project.findUnique({
    where: { id },
    select: {
      id: true,
      title: true,
      description: true,
      species: true,
      datasets: {
        orderBy: { sortOrder: "asc" },
        select: {
          sortOrder: true,
          dataset: {
            select: {
              id: true,
              title: true,
              datasetType: true,
              adminOwned: true,
              status: true,
              numCells: true,
              numGenes: true,
              s3BaseUrl: true,
              tags: true,
              metadata: true,
              viewerConfig: true,
            },
          },
        },
      },
    },
  });

  if (!project) {
    return NextResponse.json({ error: "Not found" }, { status: 404 });
  }

  const members = project.datasets.filter((d) => d.dataset.s3BaseUrl);

  if (members.length === 0 || !members.every((d) => d.dataset.adminOwned)) {
    // Deliberately 404, not 403: whether a private project exists is itself
    // not something an anonymous caller should learn.
    return NextResponse.json({ error: "Not found" }, { status: 404 });
  }

  return NextResponse.json({
    id: project.id,
    title: project.title,
    description: project.description,
    species: project.species,
    datasets: members.map(({ sortOrder, dataset }) => ({
      id: dataset.id,
      title: dataset.title,
      datasetType: dataset.datasetType,
      numCells: dataset.numCells,
      numGenes: dataset.numGenes,
      s3BaseUrl: dataset.s3BaseUrl,
      // The rail groups by stage; tags[0] is the fallback for rows registered
      // before the metadata bag carried it.
      stage:
        (dataset.metadata as { stage?: string } | null)?.stage ??
        dataset.tags[0] ??
        null,
      // So a preview can open on the orientation its owner chose.
      camera:
        (dataset.viewerConfig as { camera?: unknown } | null)?.camera ?? null,
      sortOrder,
    })),
  });
}
