import { redirect, notFound } from "next/navigation";

import { auth } from "@/lib/auth";
import { prisma } from "@/lib/prisma";
import { AccountDatasetDetailClient } from "@/components/account/dataset-detail-client";

interface Props {
  params: Promise<{ id: string }>;
}

export const metadata = {
  title: "Dataset — MERFISHEYES",
};

export default async function AccountDatasetPage({ params }: Props) {
  const { id } = await params;
  const session = await auth();

  if (!session?.user) {
    redirect(`/auth/signin?callbackUrl=/account/datasets/${id}`);
  }

  const dataset = await prisma.dataset.findUnique({
    where: { id },
    select: {
      id: true,
      title: true,
      ownerId: true,
      adminOwned: true,
      status: true,
      datasetType: true,
      s3BaseUrl: true,
      ingestSource: true,
      thumbnailUrl: true,
      numCells: true,
      numGenes: true,
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
    },
  });

  if (!dataset) notFound();

  const isAdmin =
    session.user.role === "ADMIN" || session.user.role === "SUPER_ADMIN";

  if (dataset.ownerId !== session.user.id && !(isAdmin && dataset.adminOwned)) {
    notFound();
  }

  // Any community submission for this dataset (for the Pending/Published badge).
  const submission = await prisma.catalogDataset.findFirst({
    where: { isCommunity: true, sourceDatasetId: dataset.id },
    select: {
      id: true,
      reviewStatus: true,
      isPublished: true,
      reviewNote: true,
    },
  });

  // Resolve the viewer URL the same way the list API does.
  const viewerPath =
    dataset.datasetType === "single_molecule" ? "sm-viewer" : "viewer";
  const viewerUrl =
    dataset.status !== "COMPLETE"
      ? null
      : dataset.ingestSource === "s3_registered" && dataset.s3BaseUrl
        ? `/${viewerPath}/from-s3?url=${encodeURIComponent(dataset.s3BaseUrl)}`
        : `/${viewerPath}/${dataset.id}`;

  return (
    <AccountDatasetDetailClient
      initialDataset={JSON.parse(JSON.stringify(dataset))}
      submission={submission ? JSON.parse(JSON.stringify(submission)) : null}
      viewerUrl={viewerUrl}
    />
  );
}
