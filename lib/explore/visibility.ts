import type { Prisma } from "@prisma/client";

/**
 * Which catalog rows a caller may see — the one rule for Explore's list,
 * featured row, filters, detail route and pages, /api/catalog/facets and the
 * agent API (docs/specs/MCP-PROTOCOL.md §2).
 *
 * A row is visible when it is published, not internal (internal rows: admins
 * only), and curated or an approved community submission. A *list* also
 * needs one entry that opens: by s3BaseUrl, or — where community rows are in
 * scope — by datasetId (community datasets have no public URL).
 */

export type ExploreCategory =
  | "all"
  | "featured"
  | "bil"
  | "sddb"
  | "community"
  | "internal";

const CURATED_OR_APPROVED: Prisma.CatalogDatasetWhereInput = {
  OR: [
    { isCommunity: false },
    { AND: [{ isCommunity: true }, { reviewStatus: "APPROVED" }] },
  ],
};

/** One row, opened by id or BIL code. Admins also see internal rows. */
export function catalogRowWhere(
  isAdmin: boolean,
): Prisma.CatalogDatasetWhereInput {
  return {
    AND: [
      { isPublished: true },
      ...(isAdmin ? [] : [{ isInternal: false }]),
      CURATED_OR_APPROVED,
    ],
  };
}

/** A list in one Explore category. "internal" is the admins' tab; anyone
 * else asking for it gets "all". */
export function catalogListWhere(
  category: ExploreCategory,
  isAdmin: boolean,
): Prisma.CatalogDatasetWhereInput {
  const cat = category === "internal" && !isAdmin ? "all" : category;
  const includeCommunity = cat === "all" || cat === "community";
  const conditions: Prisma.CatalogDatasetWhereInput[] = [
    { isPublished: true },
    includeCommunity
      ? {
          entries: {
            some: {
              OR: [{ s3BaseUrl: { not: null } }, { datasetId: { not: null } }],
            },
          },
        }
      : { entries: { some: { s3BaseUrl: { not: null } } } },
  ];

  if (cat === "internal") {
    conditions.push({ isInternal: true }, { isCommunity: false });
  } else {
    conditions.push({ isInternal: false });
  }

  if (cat === "community") {
    conditions.push({ isCommunity: true }, { reviewStatus: "APPROVED" });
  } else if (cat === "all") {
    conditions.push(CURATED_OR_APPROVED);
  } else if (cat !== "internal") {
    // featured / bil / sddb: curated only
    conditions.push({ isCommunity: false });
  }

  if (cat === "featured") conditions.push({ isFeatured: true });
  if (cat === "bil") conditions.push({ isBil: true });
  // Rows synced from a SpatialData-DB catalog release (docs/specs/sddb-catalog.md).
  if (cat === "sddb") conditions.push({ source: "sddb" });

  return { AND: conditions };
}

/** The featured carousel: curated plus approved community rows flagged
 * featured (unlike the "featured" tab, which is curated only). */
export function catalogFeaturedWhere(
  isAdmin: boolean,
): Prisma.CatalogDatasetWhereInput {
  return {
    AND: [
      catalogRowWhere(isAdmin),
      {
        entries: {
          some: {
            OR: [{ s3BaseUrl: { not: null } }, { datasetId: { not: null } }],
          },
        },
      },
      { isFeatured: true },
    ],
  };
}

/** Rows whose values feed the filter dropdowns and facet lists: everything
 * "all" shows, internal rows included for admins. */
export function catalogFacetWhere(
  isAdmin: boolean,
): Prisma.CatalogDatasetWhereInput {
  return catalogRowWhere(isAdmin);
}

export function isAdminRole(role: string | null | undefined): boolean {
  return role === "ADMIN" || role === "SUPER_ADMIN";
}
