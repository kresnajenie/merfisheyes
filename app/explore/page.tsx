import { prisma } from "@/lib/prisma";
import { title } from "@/components/primitives";
import { ExplorePageClient } from "@/components/explore/explore-page-client";
import { ExploreBackground } from "@/components/explore/explore-background";

export const dynamic = "force-dynamic";

export default async function ExplorePage() {
  // SSR: fetch initial data directly from DB
  const [items, featured, speciesRaw, tissueRaw, platformRaw] = await Promise.all([
    prisma.catalogDataset.findMany({
      where: { isPublished: true },
      orderBy: [{ sortOrder: "asc" }, { createdAt: "desc" }],
      take: 20,
    }),
    prisma.catalogDataset.findMany({
      where: { isPublished: true, isFeatured: true },
      orderBy: { sortOrder: "asc" },
    }),
    prisma.catalogDataset.findMany({
      where: { isPublished: true, species: { not: null } },
      select: { species: true },
      distinct: ["species"],
      orderBy: { species: "asc" },
    }),
    prisma.catalogDataset.findMany({
      where: { isPublished: true, tissue: { not: null } },
      select: { tissue: true },
      distinct: ["tissue"],
      orderBy: { tissue: "asc" },
    }),
    prisma.catalogDataset.findMany({
      where: { isPublished: true, platform: { not: null } },
      select: { platform: true },
      distinct: ["platform"],
      orderBy: { platform: "asc" },
    }),
  ]);

  const total = await prisma.catalogDataset.count({ where: { isPublished: true } });

  const filters = {
    species: speciesRaw.map((r) => r.species).filter(Boolean) as string[],
    tissues: tissueRaw.map((r) => r.tissue).filter(Boolean) as string[],
    platforms: platformRaw.map((r) => r.platform).filter(Boolean) as string[],
  };

  // Serialize dates for client
  const serialize = (list: typeof items) =>
    JSON.parse(JSON.stringify(list));

  return (
    <>
      <ExploreBackground />
      <div className="w-full relative z-10">
        <div className="mb-8">
          <h1 className={title({ size: "lg" })}>
            <span className={title({ color: "violet", size: "lg" })}>
              Explore
            </span>{" "}
            Datasets
          </h1>
          <p className="text-default-500 mt-4">
            Browse and visualize our curated collection of spatial
            transcriptomics datasets
          </p>
        </div>

        <ExplorePageClient
          initialFeatured={serialize(featured)}
          initialFilters={filters}
          initialItems={serialize(items)}
          initialTotal={total}
        />
      </div>
    </>
  );
}
