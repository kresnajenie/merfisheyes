import { NextResponse } from "next/server";

import { agentRoute } from "@/lib/agent/route";
import { notFound } from "@/lib/agent/respond";
import { catalogDetailOut } from "@/lib/agent/shape";
import { findCatalogRow } from "@/lib/agent/visibility";

// GET /api/agent/v1/catalog/{id or BIL code} — meyes_get_dataset for a
// catalog row: metadata, entries with viewer URLs, and the store's sidecar
// files merged (docs/specs/MCP-PROTOCOL.md §4).
export const GET = agentRoute<{ id: string }>(
  ["catalog:read"],
  async ({ caller, params }) => {
    const row = await findCatalogRow(caller, params.id);

    if (!row) return notFound("catalog dataset");

    return NextResponse.json(await catalogDetailOut(row));
  },
);
