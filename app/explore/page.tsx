import { Suspense } from "react";

import { prisma } from "@/lib/prisma";
import { title } from "@/components/primitives";
import { ExplorePageClient } from "@/components/explore/explore-page-client";
import { ExploreBackground } from "@/components/explore/explore-background";

// Revalidate every 60 seconds — cached page served instantly, rebuilt in
// background. This only works because the page fetches public data only and
// never touches auth() (which would force per-request dynamic rendering).
// Admin-only internal datasets are loaded client-side via /api/auth/role +
// /api/explore?tab=internal, so they don't break the cache.
export const revalidate = 60;

const includeEntries = { entries: { orderBy: { sortOrder: "asc" as const } } };

// Base filter: published + not internal + not community + at least one
// viewable entry. Community submissions are curated separately and only appear
// on the Explore "Community" tab (fetched client-side), never the main grid.
const hasViewableEntry = { entries: { some: { s3BaseUrl: { not: null } } } };
const publicBase = {
  isPublished: true,
  isInternal: false,
  isCommunity: false,
  ...hasViewableEntry,
};

// Fetch all initial public data in a single DB round-trip.
async function loadPublicData() {
  const [items, total, featured, bil, speciesRaw, tissueRaw, platformRaw] =
    await Promise.all([
      prisma.catalogDataset.findMany({
        where: publicBase,
        include: includeEntries,
        orderBy: [{ sortOrder: "asc" }, { createdAt: "desc" }],
        take: 50,
      }),
      prisma.catalogDataset.count({ where: publicBase }),
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

  return {
    items,
    total,
    featured,
    bil,
    filters: {
      species: speciesRaw.map((r) => r.species).filter(Boolean) as string[],
      tissues: tissueRaw.map((r) => r.tissue).filter(Boolean) as string[],
      platforms: platformRaw.map((r) => r.platform).filter(Boolean) as string[],
    },
  };
}

export default async function ExplorePage() {
  // Degrade gracefully if the DB is briefly unreachable during the ISR
  // (build-time) render: serve an empty shell instead of failing the build.
  // The `ssrFailed` flag tells the client to hydrate from /api/explore, and
  // ISR regenerates a correct page within `revalidate` once the DB recovers.
  let data: Awaited<ReturnType<typeof loadPublicData>> | null = null;

  try {
    data = await loadPublicData();
  } catch (err) {
    // eslint-disable-next-line no-console
    console.error("[explore] initial data load failed; serving empty shell", err);
  }

  const ssrFailed = data === null;
  const items = data?.items ?? [];
  const total = data?.total ?? 0;
  const featured = data?.featured ?? [];
  const bil = data?.bil ?? [];
  const filters = data?.filters ?? { species: [], tissues: [], platforms: [] };

  // Serialize dates for client
  const serialize = (list: unknown) => JSON.parse(JSON.stringify(list));

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

        <Suspense>
          <ExplorePageClient
            initialBil={serialize(bil)}
            initialFeatured={serialize(featured)}
            initialFilters={filters}
            initialItems={serialize(items)}
            initialTotal={total}
            ssrFailed={ssrFailed}
          />
        </Suspense>
      </div>
    </>
  );
}
