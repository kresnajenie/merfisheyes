import type { McpServer } from "@modelcontextprotocol/server";

import * as z from "zod/v4";

import { apiGet } from "../api.js";
import {
  type Page,
  join,
  many,
  n,
  one,
  paging,
  responseFormat,
} from "../format.js";

/**
 * The public catalog: Explore's search, its facets, one dataset in detail,
 * its genes (MCP-PROTOCOL §3.1). No sign-in needed; signed in, the person's
 * own datasets open too.
 */

interface Entry {
  label: string;
  kind: string;
  datasetType: string;
  datasetId: string | null;
  url: string | null;
  viewerUrl: string | null;
}

interface Card {
  id: string;
  title: string;
  species: string | null;
  tissue: string | null;
  platform: string | null;
  bilCode: string | null;
  cells: number | null;
  genes: number | null;
  source: string;
  license: string | null;
  matchedGenes?: string[];
  pageUrl: string;
  entries: Entry[];
}

interface Summary {
  cells: number | null;
  genes: number | null;
  molecules: number | null;
  sections: number | null;
  dims: number | null;
  images: string[];
  deColumns: string[];
  qc: {
    medianCountsPerCell: number | null;
    genesDetected: number | null;
    maskPercentile: number | null;
    pctArtifact: number | null;
  } | null;
  taxonomy: {
    method: string | null;
    levels: {
      level: string;
      groups: number | null;
      meanBootstrap: number | null;
    }[];
  } | null;
  investigator: string | null;
  license: string | null;
}

interface Detail {
  type: "catalog" | "dataset";
  id: string;
  title: string | null;
  species: string | null;
  tissue: string | null;
  platform: string | null;
  description: string | null;
  geneCount: number;
  cells?: number | null;
  genes?: number | null;
  status?: string;
  kind?: string;
  owned?: boolean;
  bilCode?: string | null;
  pageUrl?: string;
  viewerUrl?: string | null;
  entries?: Entry[];
  externalLink: string | null;
  publicationLink: string | null;
  sidecars: { folder: string; summary: Summary }[];
}

/** App datasets have prefixed ids; anything else is a catalog row or BIL code. */
export const isDatasetId = (id: string) => /^(ds|sm|s3)_/.test(id);

const cardMd = (c: Card) =>
  [
    `## ${c.title} (${c.id})`,
    `- ${join(c.species, c.tissue, c.platform, c.bilCode && `BIL ${c.bilCode}`) || "no species/tissue/platform given"}`,
    `- ${n(c.cells)} cells · ${n(c.genes)} genes${c.license ? ` · licence ${c.license}` : ""}`,
    ...(c.matchedGenes?.length
      ? [`- Matching genes: ${c.matchedGenes.join(", ")}`]
      : []),
    ...c.entries.map(
      (e) => `- ${e.label} (${e.kind}): ${e.viewerUrl ?? "no viewer link"}`,
    ),
  ].join("\n");

function detailMd(d: Detail): string {
  const lines = [`# ${d.title ?? d.id} (${d.id})`];

  lines.push(
    join(
      d.species,
      d.tissue,
      d.platform,
      d.bilCode && `BIL ${d.bilCode}`,
      d.status && `status ${d.status}`,
    ) || "",
  );
  if (d.type === "dataset")
    lines.push(
      `App dataset, ${d.kind}${d.owned ? ", yours" : ""}: ${d.viewerUrl ?? "not viewable yet"}`,
    );
  else lines.push(`Catalog page: ${d.pageUrl}`);
  for (const e of d.entries ?? [])
    lines.push(`- ${e.label} (${e.kind}): ${e.viewerUrl ?? "—"}`);
  lines.push(`Genes: ${n(d.geneCount)} (list them with meyes_list_genes)`);
  for (const s of d.sidecars) {
    const m = s.summary;

    lines.push(
      "",
      "## From the store's own files",
      `- ${n(m.cells)} cells · ${n(m.genes)} genes · ${n(m.molecules)} molecules · ${n(m.sections)} sections · ${m.dims ?? "—"}D`,
      `- DE stats for columns: ${m.deColumns.join(", ") || "none"} (not every column is a cell type)`,
    );
    if (m.taxonomy) {
      lines.push(
        `- Cell types (${m.taxonomy.method ?? "MapMyCells"}): ` +
          m.taxonomy.levels
            .map(
              (l) =>
                `${l.level} ${n(l.groups)}${l.meanBootstrap !== null ? ` (bootstrap ${l.meanBootstrap})` : ""}`,
            )
            .join(", "),
      );
    }
    if (m.qc) {
      lines.push(
        `- QC: median ${n(m.qc.medianCountsPerCell)} counts/cell · ${n(m.qc.genesDetected)} genes detected` +
          (m.qc.pctArtifact !== null
            ? ` · ${m.qc.pctArtifact}% masked as artifact (p${m.qc.maskPercentile})`
            : ""),
      );
    }
    if (m.images.length) lines.push(`- Images: ${m.images.join(", ")}`);
    if (m.investigator || m.license)
      lines.push(
        `- ${join(m.investigator && `Investigator ${m.investigator}`, m.license && `licence ${m.license}`)}`,
      );
  }
  if (d.description)
    lines.push(
      "",
      d.description.length > 1500
        ? `${d.description.slice(0, 1500)}…`
        : d.description,
    );
  if (d.publicationLink || d.externalLink)
    lines.push(
      "",
      join(
        d.publicationLink && `Publication: ${d.publicationLink}`,
        d.externalLink && `Source: ${d.externalLink}`,
      ),
    );

  return lines.filter((l) => l !== "").join("\n");
}

const READ = {
  readOnlyHint: true,
  destructiveHint: false,
  idempotentHint: true,
  openWorldHint: true,
};

export function registerCatalogTools(server: McpServer) {
  server.registerTool(
    "meyes_search_datasets",
    {
      title: "Search datasets",
      description: `Search the MERFISHEYES public catalog (Explore) of spatial transcriptomics datasets.
Filters combine: free text (title, description, BIL code, investigator, tags), species, tissue, platform, and genes the dataset measures. Returns one page of cards, each with its viewer links; no gene lists (use meyes_list_genes).
Spell species/tissue/platform as meyes_list_facets returns them.
Examples:
  - "mouse brain MERFISH datasets" -> species="Mouse", tissue="Brain", platform="MERFISH"
  - "which datasets measure Gad1 and Slc17a7" -> genes_exact=["Gad1","Slc17a7"]
  - "anything from the Allen Institute" -> query="Allen"
Use meyes_get_dataset for one dataset in detail.`,
      inputSchema: z.strictObject({
        query: z
          .string()
          .max(200)
          .optional()
          .describe("Free text; every word must match some field"),
        species: z
          .string()
          .max(200)
          .optional()
          .describe("Exact species, as meyes_list_facets spells it"),
        tissue: z
          .string()
          .max(200)
          .optional()
          .describe("Exact tissue, as meyes_list_facets spells it"),
        platform: z
          .string()
          .max(200)
          .optional()
          .describe("Exact platform, e.g. MERFISH, Xenium"),
        genes: z
          .array(z.string().min(1).max(50))
          .max(10)
          .optional()
          .describe("Gene name fragments; every one must match some gene"),
        genes_exact: z
          .array(z.string().min(1).max(50))
          .max(10)
          .optional()
          .describe("Gene symbols that must all be measured, case-insensitive"),
        category: z
          .enum(["all", "featured", "bil", "sddb", "community"])
          .default("all")
          .describe(
            "Which part of Explore: all, featured, BIL, SpatialData-DB, community submissions",
          ),
        ...paging,
        response_format: responseFormat,
      }),
      annotations: READ,
    },
    async (a, ctx) => {
      const page = await apiGet<Page<Card>>(
        "/catalog/search",
        {
          q: a.query,
          species: a.species,
          tissue: a.tissue,
          platform: a.platform,
          genes: a.genes?.join(","),
          genes_exact: a.genes_exact?.join(","),
          category: a.category,
          limit: a.limit,
          offset: a.offset,
        },
        ctx.mcpReq.signal,
      );

      return many(a.response_format, page, {
        heading: "MERFISHEYES datasets",
        item: cardMd,
        empty:
          "No dataset matches. Loosen a filter, check the spelling with meyes_list_facets, or search genes with `genes` (fragments) instead of `genes_exact`.",
      });
    },
  );

  server.registerTool(
    "meyes_list_facets",
    {
      title: "List species, tissues and platforms",
      description:
        'The species, tissue and platform values used in the MERFISHEYES catalog, spelled exactly as meyes_search_datasets filters on them. Call before filtering when unsure of a spelling ("Mouse" or "mouse"? "Brain" or "brain"?).',
      inputSchema: z.strictObject({ response_format: responseFormat }),
      annotations: READ,
    },
    async ({ response_format }, ctx) => {
      const f = await apiGet<{
        species: string[];
        tissues: string[];
        platforms: string[];
      }>("/catalog/facets", {}, ctx.mcpReq.signal);

      return one(response_format, f, () =>
        [
          `# Catalog values`,
          `**Species:** ${f.species.join(", ")}`,
          `**Tissues:** ${f.tissues.join(", ")}`,
          `**Platforms:** ${f.platforms.join(", ")}`,
        ].join("\n\n"),
      );
    },
  );

  server.registerTool(
    "meyes_get_dataset",
    {
      title: "Get a dataset",
      description: `One dataset in detail: what it is, its viewer links, and — when its store carries them — the numbers its own files record (cells, genes, molecules, sections, which columns have DE stats, cell-type levels, QC).
Takes a catalog id or BIL code (from meyes_search_datasets) or an app dataset id (ds_…, sm_…, s3_…, from meyes_list_my_datasets or a project).
Examples: id="ace-dip-vat", id="cm1x…", id="ds_a1B2c3D4e5".`,
      inputSchema: z.strictObject({
        id: z
          .string()
          .min(1)
          .max(200)
          .describe("Catalog id, BIL code, or app dataset id"),
        response_format: responseFormat,
      }),
      annotations: READ,
    },
    async ({ id, response_format }, ctx) => {
      const path = isDatasetId(id)
        ? `/datasets/${encodeURIComponent(id)}`
        : `/catalog/${encodeURIComponent(id)}`;
      const d = await apiGet<Detail>(path, {}, ctx.mcpReq.signal);

      return one(response_format, d, () => detailMd(d));
    },
  );

  server.registerTool(
    "meyes_list_genes",
    {
      title: "List a dataset's genes",
      description: `The genes a dataset measures, filtered by a name fragment and paged (100 a page, up to 1,000). Use it to check a gene is measured before asking about it.
Examples: id="ace-dip-vat", query="Gad" -> Gad1, Gad2.`,
      inputSchema: z.strictObject({
        id: z
          .string()
          .min(1)
          .max(200)
          .describe("Catalog id, BIL code, or app dataset id"),
        query: z
          .string()
          .max(50)
          .optional()
          .describe("Case-insensitive fragment of a gene name"),
        limit: z
          .number()
          .int()
          .min(1)
          .max(1000)
          .default(100)
          .describe("How many to return (1–1000)"),
        offset: z
          .number()
          .int()
          .min(0)
          .default(0)
          .describe("How many to skip; use next_offset"),
        response_format: responseFormat,
      }),
      annotations: READ,
    },
    async ({ id, query, limit, offset, response_format }, ctx) => {
      const path = isDatasetId(id)
        ? `/datasets/${encodeURIComponent(id)}/genes`
        : `/catalog/${encodeURIComponent(id)}/genes`;
      const page = await apiGet<Page<string>>(
        path,
        { q: query, limit, offset },
        ctx.mcpReq.signal,
      );

      if (response_format === "json")
        return many("json", page, { heading: "", item: String, empty: "" });
      if (page.total === 0)
        return {
          content: [
            {
              type: "text" as const,
              text: query
                ? `No gene of ${id} contains "${query}".`
                : `${id} lists no genes.`,
            },
          ],
        };
      const next = page.offset + page.items.length;

      return {
        content: [
          {
            type: "text" as const,
            text: `${page.total} gene(s)${query ? ` matching "${query}"` : ""}, showing ${page.items.length}:\n${page.items.join(", ")}${page.has_more ? `\n\nMore: call again with offset=${next}.` : ""}`,
          },
        ],
      };
    },
  );
}
