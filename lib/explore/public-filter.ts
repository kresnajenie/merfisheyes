/**
 * What Explore shows publicly: curated content plus approved community
 * submissions. Community datasets are viewable by datasetId (no public
 * s3BaseUrl), so the viewable-entry check accepts either.
 */
const hasViewableEntryAny = {
  entries: {
    some: { OR: [{ s3BaseUrl: { not: null } }, { datasetId: { not: null } }] },
  },
};

export const publicCatalogWhere = {
  isPublished: true,
  isInternal: false,
  OR: [
    { isCommunity: false },
    { AND: [{ isCommunity: true }, { reviewStatus: "APPROVED" as const }] },
  ],
  ...hasViewableEntryAny,
};
