import { prisma } from "@/lib/prisma";
import { auth } from "@/lib/auth";
import { title } from "@/components/primitives";
import { ExplorePageClient } from "@/components/explore/explore-page-client";
import { ExploreBackground } from "@/components/explore/explore-background";

export const dynamic = "force-dynamic";

const includeEntries = { entries: { orderBy: { sortOrder: "asc" as const } } };

export default async function ExplorePage() {
  const session = await auth();
  const isAdmin = session?.user?.role === "ADMIN" || session?.user?.role === "SUPER_ADMIN";

  // Base filter: published + not internal (for public queries)
  const publicBase = { isPublished: true, isInternal: false };

  // SSR: fetch initial data directly from DB
  const [items, featured, bil, speciesRaw, tissueRaw, platformRaw] = await Promise.all([
    prisma.catalogDataset.findMany({
      where: publicBase,
      include: includeEntries,
      orderBy: [{ sortOrder: "asc" }, { createdAt: "desc" }],
      take: 50,
    }),
    prisma.catalogDataset.findMany({
      where: { ...publicBase, isFeatured: true },
      include: includeEntries,
      orderBy: { sortOrder: "asc" },
    }),
    prisma.catalogDataset.findMany({
      where: { ...publicBase, isBil: true },
      include: includeEntries,
      orderBy: { sortOrder: "asc" },
    }),
    prisma.catalogDataset.findMany({
      where: { ...publicBase, species: { not: null } },
      select: { species: true },
      distinct: ["species"],
      orderBy: { species: "asc" },
    }),
    prisma.catalogDataset.findMany({
      where: { ...publicBase, tissue: { not: null } },
      select: { tissue: true },
      distinct: ["tissue"],
      orderBy: { tissue: "asc" },
    }),
    prisma.catalogDataset.findMany({
      where: { ...publicBase, platform: { not: null } },
      select: { platform: true },
      distinct: ["platform"],
      orderBy: { platform: "asc" },
    }),
  ]);

  const total = await prisma.catalogDataset.count({ where: publicBase });

  // Fetch internal datasets only for admins
  const internal = isAdmin
    ? await prisma.catalogDataset.findMany({
        where: { isPublished: true, isInternal: true },
        include: includeEntries,
        orderBy: [{ sortOrder: "asc" }, { createdAt: "desc" }],
      })
    : [];

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
          initialBil={serialize(bil)}
          initialFeatured={serialize(featured)}
          initialFilters={filters}
          initialInternal={serialize(internal)}
          initialItems={serialize(items)}
          initialTotal={total}
          isAdmin={isAdmin}
        />
      </div>
    </>
  );
}
