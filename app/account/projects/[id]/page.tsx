import { redirect, notFound } from "next/navigation";

import { auth } from "@/lib/auth";
import { prisma } from "@/lib/prisma";
import { ProjectDetailClient } from "@/components/account/project-detail-client";

interface Props {
  params: Promise<{ id: string }>;
}

export const metadata = {
  title: "Project — MERFISHEYES",
};

export default async function ProjectPage({ params }: Props) {
  const { id } = await params;
  const session = await auth();

  if (!session?.user) {
    redirect(`/auth/signin?callbackUrl=/account/projects/${id}`);
  }

  const project = await prisma.project.findUnique({
    where: { id },
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
              numCells: true,
              numGenes: true,
              ingestSource: true,
              s3BaseUrl: true,
            },
          },
        },
      },
    },
  });

  if (!project) notFound();

  const isAdmin =
    session.user.role === "ADMIN" || session.user.role === "SUPER_ADMIN";

  if (project.ownerId !== session.user.id && !isAdmin) {
    notFound();
  }

  // The owner's completed datasets, for the "add datasets" picker.
  const ownDatasets = await prisma.dataset.findMany({
    where: {
      status: "COMPLETE",
      ...(isAdmin
        ? { OR: [{ ownerId: project.ownerId }, { adminOwned: true }] }
        : { ownerId: session.user.id }),
    },
    select: {
      id: true,
      title: true,
      datasetType: true,
      thumbnailUrl: true,
    },
    orderBy: { createdAt: "desc" },
    take: 500,
  });

  const submission = await prisma.catalogDataset.findFirst({
    where: { isCommunity: true, sourceProjectId: project.id },
    select: {
      id: true,
      reviewStatus: true,
      isPublished: true,
      reviewNote: true,
    },
  });

  return (
    <ProjectDetailClient
      initialProject={JSON.parse(JSON.stringify(project))}
      ownDatasets={JSON.parse(JSON.stringify(ownDatasets))}
      submission={submission ? JSON.parse(JSON.stringify(submission)) : null}
    />
  );
}
