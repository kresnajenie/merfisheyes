import { notFound } from "next/navigation";

import { prisma } from "@/lib/prisma";
import { DatasetDetailPage } from "@/components/explore/dataset-detail-page";
import { ExploreEditButton } from "@/components/explore/explore-edit-button";

interface Props {
  params: Promise<{ id: string }>;
}

export default async function DatasetPage({ params }: Props) {
  const { id } = await params;

  const dataset = await prisma.catalogDataset.findUnique({
    where: { id },
    include: { entries: { orderBy: { sortOrder: "asc" } } },
  });

  if (!dataset || !dataset.isPublished) {
    notFound();
  }

  return (
    <DatasetDetailPage
      dataset={JSON.parse(JSON.stringify(dataset))}
      headerActions={<ExploreEditButton catalogId={dataset.id} />}
    />
  );
}
