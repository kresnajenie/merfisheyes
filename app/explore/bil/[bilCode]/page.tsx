import { notFound } from "next/navigation";

import { prisma } from "@/lib/prisma";
import { DatasetDetailPage } from "@/components/explore/dataset-detail-page";

interface Props {
  params: Promise<{ bilCode: string }>;
}

export default async function BilDatasetPage({ params }: Props) {
  const { bilCode } = await params;

  const dataset = await prisma.catalogDataset.findUnique({
    where: { bilCode },
    include: { entries: { orderBy: { sortOrder: "asc" } } },
  });

  if (!dataset || !dataset.isPublished) {
    notFound();
  }

  return <DatasetDetailPage dataset={JSON.parse(JSON.stringify(dataset))} />;
}
