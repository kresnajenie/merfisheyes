import { notFound } from "next/navigation";

import { prisma } from "@/lib/prisma";
import { auth } from "@/lib/auth";
import { catalogRowWhere } from "@/lib/explore/visibility";
import { isAdminRole } from "@/lib/admin-auth";
import { DatasetDetailPage } from "@/components/explore/dataset-detail-page";
import { ExploreEditButton } from "@/components/explore/explore-edit-button";

interface Props {
  params: Promise<{ id: string }>;
}

export default async function DatasetPage({ params }: Props) {
  const { id } = await params;

  const session = await auth();
  const dataset = await prisma.catalogDataset.findFirst({
    where: { AND: [{ id }, catalogRowWhere(isAdminRole(session?.user?.role))] },
    include: { entries: { orderBy: { sortOrder: "asc" } } },
  });

  if (!dataset) {
    notFound();
  }

  return (
    <DatasetDetailPage
      dataset={JSON.parse(JSON.stringify(dataset))}
      headerActions={<ExploreEditButton catalogId={dataset.id} />}
    />
  );
}
