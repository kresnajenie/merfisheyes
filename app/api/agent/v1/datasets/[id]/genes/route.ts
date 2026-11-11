import { NextResponse } from "next/server";

import { agentRoute } from "@/lib/agent/route";
import { notFound, page, paging } from "@/lib/agent/respond";
import { genesPage } from "@/lib/agent/shape";
import { findDataset } from "@/lib/agent/visibility";

// GET /api/agent/v1/datasets/{id}/genes?q= — meyes_list_genes
export const GET = agentRoute<{ id: string }>(
  ["catalog:read"],
  async ({ url, caller, params }) => {
    const dataset = await findDataset(caller, params.id);

    if (!dataset) return notFound("dataset");
    const p = paging(url, 100, 1000);
    const { hits, total } = genesPage(
      dataset.genes,
      url.searchParams.get("q") ?? "",
      p,
    );

    return NextResponse.json(page(hits, total, p));
  },
);
