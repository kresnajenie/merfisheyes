import { NextResponse } from "next/server";

import { agentRoute } from "@/lib/agent/route";
import { catalogFacets } from "@/lib/explore/search";

// GET /api/agent/v1/catalog/facets — meyes_list_facets: the spellings
// in use, so a filter matches.
export const GET = agentRoute(["catalog:read"], async ({ caller }) =>
  NextResponse.json(
    await catalogFacets(caller.kind === "user" && caller.isAdmin),
  ),
);
