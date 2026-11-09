import type { ExploreCategory } from "@/lib/explore/visibility";

import { NextResponse } from "next/server";

import { agentRoute } from "@/lib/agent/route";
import { page, paging, problem } from "@/lib/agent/respond";
import { cardOut } from "@/lib/agent/shape";
import { searchCatalog } from "@/lib/explore/search";

const CATEGORIES = ["all", "featured", "bil", "sddb", "community", "internal"];

// GET /api/agent/v1/catalog/search — meyes_search_datasets
// (docs/specs/MCP-PROTOCOL.md §3.1): Explore's own search.
export const GET = agentRoute(["catalog:read"], async ({ url, caller }) => {
  const q = url.searchParams;
  const category = (q.get("category") || "all") as ExploreCategory;

  if (!CATEGORIES.includes(category)) {
    return problem(
      400,
      "invalid_request",
      `category is one of: ${CATEGORIES.slice(0, 5).join(", ")}`,
    );
  }
  const p = paging(url);
  const { items, total } = await searchCatalog(
    {
      search: q.get("q") ?? "",
      species: q.get("species") ?? "",
      tissue: q.get("tissue") ?? "",
      platform: q.get("platform") ?? "",
      genes: q.get("genes") ?? "",
      genesExact: q.get("genes_exact") ?? "",
      datasetType: q.get("dataset_type") ?? "",
      category,
      skip: p.offset,
      take: p.limit,
    },
    caller.kind === "user" && caller.isAdmin,
  );

  return NextResponse.json(page(items.map(cardOut), total, p));
});
