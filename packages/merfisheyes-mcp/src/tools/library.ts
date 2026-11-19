import type { McpServer } from "@modelcontextprotocol/server";

import * as z from "zod/v4";

import { CHARACTER_LIMIT } from "../config.ts";
import { join, n, one, responseFormat, text } from "../format.ts";
import { libraryIndex } from "../library/data.ts";
import {
  type Filters,
  type Hit,
  type SearchResult,
  VALUE_COLUMNS,
  hit,
  search,
  values,
} from "../library/search.ts";

/**
 * Research the library (MCP.md pillar B): every live MERFISHEYES dataset —
 * ours and SpatialData-DB's — searched by its metadata and genes, from the
 * published index. One call answers a count, a breakdown and the first
 * rows; no sign-in, and no read of any dataset's own files.
 */

const READ = {
  readOnlyHint: true,
  destructiveHint: false,
  idempotentHint: true,
  openWorldHint: true,
};

const list = z
  .union([z.string(), z.array(z.string())])
  .transform((v) => (Array.isArray(v) ? v : [v]))
  .optional();

const filterShape = {
  organism: list.describe(
    'Any of these organisms, e.g. "mouse", "human", "macaque", or an NCBITaxon id',
  ),
  tissue: list.describe(
    'Any of these tissues, by name or UBERON id, each with everything below it: "brain" finds the hippocampus. Organ and system groupings work too ("nervous system")',
  ),
  disease: list.describe(
    'Any of these diseases, by name or MONDO id, with everything below it; also a class or group ("cancer", "neurodegenerative disease") or "normal"',
  ),
  technology: list.describe(
    'Any of these technologies, by id or name: "10x_xenium", "MERSCOPE", "Visium HD", "CosMx"… ("merfish" is MERFISH outside Vizgen\'s MERSCOPE; ask for both when the person means either)',
  ),
  source: list.describe(
    'Where the dataset comes from: "sddb" (SpatialData-DB), "bil" (Brain Image Library), "curated", "community"',
  ),
  kind: list.describe('"single_cell" (cells) or "single_molecule" (molecules)'),
  metadata_from: list.describe(
    'How far to trust the metadata: "store" (in the dataset itself, checked against the ontologies), "sddb" (SpatialData-DB\'s curation), "explore" (the site\'s free text, not checked)',
  ),
  modality: z
    .string()
    .optional()
    .describe('"RNA", "protein", "Metabolomics", "Epigenomics"'),
  sex: z.enum(["male", "female"]).optional(),
  age_days_min: z
    .number()
    .optional()
    .describe(
      "Youngest age in days (8 weeks = 56); rows with no age are left out",
    ),
  age_days_max: z.number().optional().describe("Oldest age in days"),
  genotype: z
    .string()
    .optional()
    .describe('Words the genotype contains, e.g. "APP", "C57BL/6"'),
  condition: z.string().optional().describe("Words the condition contains"),
  treatment: z.string().optional().describe("Words the treatment contains"),
  text: z
    .string()
    .optional()
    .describe(
      "Words that must all appear in the title, description, study, DOI, BIL code, investigator or sample",
    ),
  min_cells: z
    .number()
    .int()
    .optional()
    .describe("At least this many cells (or spots, bins)"),
  genes_measured: list.describe(
    'Genes the dataset measures (all of them), by symbol or Ensembl id, case-insensitive: ["Gad1","Slc17a7"]. Covers every dataset',
  ),
  genes_expressed: z
    .array(
      z.strictObject({
        gene: z.string(),
        min_fraction: z
          .number()
          .min(0)
          .max(1)
          .default(0.05)
          .describe("Fraction of cells detecting it, 0–1 (default 0.05)"),
      }),
    )
    .optional()
    .describe(
      "Genes detected in at least a fraction of the dataset's cells. Counts are recorded for MERFISHEYES's own datasets only for now; the answer says how many it could not check",
    ),
  collapse_rereleases: z
    .boolean()
    .default(true)
    .describe(
      "Fold re-releases of the same data (a newer Space Ranger run, say) into the newest",
    ),
};

function filtersOf(a: Record<string, unknown>): Filters {
  const f = { ...a } as Filters & {
    genes_measured?: string[];
    genes_expressed?: { gene: string; min_fraction: number }[];
  };
  const genes = [
    ...(f.genes_measured ?? []).map((gene) => ({ gene, min_fraction: 0 })),
    ...(f.genes_expressed ?? []),
  ];

  delete f.genes_measured;
  delete f.genes_expressed;

  return genes.length ? { ...f, genes } : f;
}

const TRUST: Record<string, string> = {
  store: "metadata checked, in the store",
  sddb: "SpatialData-DB's curation",
  explore: "Explore's text, not checked",
};

const hitMd = (h: Hit) =>
  [
    `### ${h.title ?? h.id} (${h.id})`,
    `- ${join(h.organism, h.tissue, h.disease && h.disease !== "unknown" && h.disease, h.technology) || "—"}`,
    `- ${n(h.cells)} ${h.obs_unit ?? "cells"} · ${n(h.genes_measured)} features` +
      (h.timepoint ? ` · ${h.timepoint}` : "") +
      (h.sex ? ` · ${h.sex}` : ""),
    ...(h.genotype || h.condition || h.treatment
      ? [
          `- ${join(h.genotype && `genotype ${h.genotype}`, h.condition && `condition ${h.condition}`, h.treatment && `treatment ${h.treatment}`)}`,
        ]
      : []),
    ...(h.gene_stats?.length
      ? [
          `- ${h.gene_stats
            .map(
              (g) =>
                `${g.gene}: ${g.fraction_detected === null ? "measured" : `${(g.fraction_detected * 100).toFixed(1)} % of cells`}${g.total_counts !== null ? `, ${n(Math.round(g.total_counts))} counts` : ""}`,
            )
            .join(" · ")}`,
        ]
      : []),
    `- ${join(h.source, TRUST[h.metadata_from ?? ""] ?? h.metadata_from, h.study_id, h.doi && `doi ${h.doi}`, h.bil_code && `BIL ${h.bil_code}`)}`,
    `- View: ${h.viewer_url ?? "—"} · Page: ${h.page_url ?? "—"}`,
    ...(h.also_published_as?.length
      ? [
          `- Also published as: ${h.also_published_as.map((o) => `${o.id}${o.release ? ` (${o.release})` : ""}`).join(", ")}`,
        ]
      : []),
  ].join("\n");

function render(
  format: "markdown" | "json",
  res: SearchResult,
  offset: number,
  limit: number,
) {
  const build = (k: number) => {
    const items = res.rows.slice(offset, offset + k).map((r) => hit(r, res));
    const next = offset + items.length;
    const page = {
      release: res.release,
      total: res.total,
      notes: res.notes,
      breakdown: res.breakdown,
      count: items.length,
      offset,
      has_more: res.total > next,
      ...(res.total > next ? { next_offset: next } : {}),
      items,
    };

    if (format === "json") return { page, body: JSON.stringify(page, null, 2) };
    const lines = [
      `# ${n(res.total)} dataset${res.total === 1 ? "" : "s"}`,
      `Index ${res.release}.`,
      ...res.notes.map((x) => `- ${x}`),
    ];

    if (res.total) {
      lines.push("", "## Breakdown");
      for (const [facet, vals] of Object.entries(res.breakdown))
        lines.push(
          `- ${facet.replace(/^meyes_/, "")}: ${vals.map((v) => `${v.value} ${n(v.count)}`).join(" · ")}`,
        );
      lines.push(
        "",
        `## Datasets ${offset + 1}–${next} of ${n(res.total)}`,
        "",
        items.map(hitMd).join("\n\n"),
      );
      if (page.has_more)
        lines.push("", `More: call again with offset=${next}.`);
    }

    return { page, body: lines.join("\n") };
  };
  let k = limit;
  let r = build(k);

  while (k > 1 && r.body.length > CHARACTER_LIMIT) {
    k = Math.ceil(k / 2);
    r = build(k);
  }

  return format === "json"
    ? {
        content: [{ type: "text" as const, text: r.body }],
        structuredContent: r.page as unknown as Record<string, unknown>,
      }
    : text(r.body);
}

function failure(e: unknown) {
  return {
    ...text(e instanceof Error ? e.message : String(e)),
    isError: true,
  };
}

export function registerLibraryTools(server: McpServer) {
  server.registerTool(
    "meyes_find_datasets",
    {
      title: "Find datasets in the library",
      description: `Search every public dataset on MERFISHEYES (its own and SpatialData-DB's, ~12,700) by metadata and genes, and answer how many there are, how they break down, and the first ones — in one call. Use it for "how many…", "which…", "find datasets with…" questions; meyes_get_dataset then gives one in detail.
Filters combine with AND; a list inside one filter is OR. Tissue and disease include everything below them in the ontology. Filters on rarely recorded fields (age, sex, genotype, condition, treatment) leave out rows without them and the notes say how many.
Each row says where its metadata comes from (checked store metadata, SpatialData-DB's curation, or unchecked Explore text) and has a viewer link to give the person.
Examples:
  - "mouse brain datasets measuring Gad1 and Slc17a7" -> organism="mouse", tissue="brain", genes_measured=["Gad1","Slc17a7"]
  - "human Alzheimer's datasets on Xenium" -> organism="human", disease="Alzheimer disease", technology="10x_xenium"
  - "adult mouse (P56 or older) cortex expressing Rab3b in >5% of cells" -> organism="mouse", tissue="cerebral cortex", age_days_min=56, genes_expressed=[{gene:"Rab3b", min_fraction:0.05}]
Spell values as meyes_index_values lists them when unsure; an unknown value returns the closest ones.`,
      inputSchema: z.strictObject({
        ...filterShape,
        sort: z
          .enum(["newest", "cells", "title"])
          .default("newest")
          .describe("Order of the rows"),
        limit: z.number().int().min(1).max(50).default(10),
        offset: z.number().int().min(0).default(0),
        response_format: responseFormat,
      }),
      annotations: READ,
    },
    async (a) => {
      try {
        const { sort, limit, offset, response_format, ...rest } = a;
        const index = await libraryIndex();
        const res = await search(index, filtersOf(rest), sort);

        return render(response_format, res, offset, limit);
      } catch (e) {
        return failure(e);
      }
    },
  );

  server.registerTool(
    "meyes_index_values",
    {
      title: "List the values of a field",
      description: `The values one field takes across the library, most common first, with counts — optionally among the datasets that match filters (the same filters as meyes_find_datasets). Use it to spell a filter value, or to answer "which tissues/technologies/diseases are there for X".`,
      inputSchema: z.strictObject({
        column: z.enum(VALUE_COLUMNS).describe("The field to list"),
        containing: z
          .string()
          .optional()
          .describe("Only values containing this text, case-insensitive"),
        ...filterShape,
        limit: z.number().int().min(1).max(200).default(50),
        response_format: responseFormat,
      }),
      annotations: READ,
    },
    async (a) => {
      try {
        const { column, containing, limit, response_format, ...rest } = a;
        const index = await libraryIndex();
        const f = filtersOf(rest);
        const filtered = Object.values(f).some(
          (v) => v !== undefined && v !== true,
        );
        const res = filtered
          ? await search(index, f)
          : { rows: index.rows, total: index.rows.length, notes: [] };
        const all = values(res.rows, column, containing);
        const items = all.slice(0, limit);
        const data = {
          release: index.release,
          column,
          datasets: res.total,
          distinct: all.length,
          notes: res.notes,
          items,
        };

        return one(response_format, data, () =>
          [
            `# ${column}: ${n(all.length)} value${all.length === 1 ? "" : "s"} over ${n(res.total)} datasets`,
            ...res.notes.map((x) => `- ${x}`),
            "",
            ...items.map((v) => `- ${v.value} (${n(v.count)})`),
            ...(all.length > items.length
              ? [
                  `- … ${n(all.length - items.length)} more; raise limit or use containing`,
                ]
              : []),
          ].join("\n"),
        );
      } catch (e) {
        return failure(e);
      }
    },
  );
}
