import { NextRequest, NextResponse } from "next/server";

import { auth } from "@/lib/auth";
import {
  catalogFacets,
  featuredCatalog,
  searchCatalog,
} from "@/lib/explore/search";
import { type ExploreCategory } from "@/lib/explore/visibility";
import { isAdminRole } from "@/lib/admin-auth";

// GET /api/explore — public endpoint for browsing published datasets.
//
// Response shape: { items, total, page, limit } — items are card-slim rows
// (no genes array; see lib/explore/query.ts). Pass `include=meta` to also get
// { featured, filters }; that payload changes rarely, so consumers fetch it
// once (first load / modal open), not on every keystroke.
//
// The search and the visibility rule live in lib/explore/ (search.ts,
// visibility.ts), shared with the agent API.
export async function GET(req: NextRequest) {
  const session = await auth();
  const isAdmin = isAdminRole(session?.user?.role);

  const url = new URL(req.url);
  const page = Math.max(1, Number(url.searchParams.get("page") ?? "1"));
  const limit = Math.min(
    100,
    Math.max(1, Number(url.searchParams.get("limit") ?? "20")),
  );

  // Category filter (was tabs): "all" | "featured" | "bil" | "sddb" |
  // "community" | "internal". "all" includes approved community submissions
  // alongside curated content.
  const { items, total } = await searchCatalog(
    {
      search: url.searchParams.get("search") ?? "",
      species: url.searchParams.get("species") ?? "",
      tissue: url.searchParams.get("tissue") ?? "",
      platform: url.searchParams.get("platform") ?? "",
      genes: url.searchParams.get("genes") ?? "",
      genesExact: url.searchParams.get("genesExact") ?? "",
      datasetType: url.searchParams.get("datasetType") ?? "",
      category: (url.searchParams.get("tab") || "all") as ExploreCategory,
      skip: (page - 1) * limit,
      take: limit,
    },
    isAdmin,
  );

  if (url.searchParams.get("include") !== "meta") {
    return NextResponse.json({ items, total, page, limit });
  }

  const [featured, filters] = await Promise.all([
    featuredCatalog(isAdmin),
    catalogFacets(isAdmin),
  ]);

  return NextResponse.json({ items, total, page, limit, featured, filters });
}
