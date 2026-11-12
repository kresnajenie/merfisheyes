import type { McpServer } from "@modelcontextprotocol/server";

import * as z from "zod/v4";

import { apiGet } from "../api.js";
import { openBrowser } from "../auth.js";
import { SITE } from "../config.js";
import { type Page, text } from "../format.js";
import { isDatasetId } from "./catalog.js";

/**
 * meyes_open_view: open a dataset in the person's browser already set up —
 * a gene coloured, a cell-type column chosen, 2D or 3D. No bridge to a tab:
 * the viewers read their settings from the link (`v=`, base64url JSON,
 * lib/utils/url-visualization-state.ts), so the tool builds that link, after
 * correcting gene names to the dataset's own spelling (the viewer only applies
 * an exact match), and opens it.
 */

interface Entry {
  label: string;
  kind: string;
  viewerUrl: string | null;
}

interface Detail {
  id: string;
  title: string | null;
  kind?: string;
  viewerUrl?: string | null;
  entries?: Entry[];
}

/** Which settings a viewer reads from `v=`. */
type ViewerKind = "cells" | "molecules" | "none";

export function viewerKind(url: URL): ViewerKind {
  const p = url.pathname;

  if (p.startsWith("/sm-viewer")) return "molecules";
  if (p.startsWith("/viewer") || p.startsWith("/spatialdata-viewer"))
    return "cells";

  return "none";
}

/** Distinct and bright: the viewer's background is black. */
const GENE_COLORS = [
  "#00d8ff",
  "#ff3da5",
  "#ffd400",
  "#7cff4f",
  "#ff8a00",
  "#b18cff",
  "#ff5c5c",
  "#40ffc0",
  "#ffffff",
  "#4e9bff",
];

/** base64url of the JSON, as toBase64Url in lib/utils/url-visualization-state.ts. */
export function toBase64Url(obj: object): string {
  return Buffer.from(JSON.stringify(obj), "utf8").toString("base64url");
}

export interface ViewSettings {
  genes: string[];
  column?: string;
  celltypes?: string[];
  mode3d?: boolean;
}

/** The `v=` state a viewer of this kind reads (CellVizUrlState /
 * SMVizUrlState in lib/utils/url-visualization-state.ts). */
export function viewState(
  kind: ViewerKind,
  s: ViewSettings,
): Record<string, unknown> | null {
  const state: Record<string, unknown> = {};

  if (kind === "cells") {
    const [g, g2] = s.genes;

    if (g) {
      state.g = g;
      state.m = ["gene"];
    }
    if (g2) {
      state.g2 = g2;
      state.ce = true;
    }
    if (s.column) state.c = s.column;
    if (s.celltypes?.length) {
      state.ct = s.celltypes;
      if (!g) state.m = ["celltype"];
    }
  } else if (kind === "molecules") {
    if (s.genes.length)
      state.genes = s.genes.map((g, i) => [
        g,
        GENE_COLORS[i % GENE_COLORS.length],
        1,
        true,
      ]);
  } else return null;
  if (s.mode3d) state.vm = "3D";

  return Object.keys(state).length ? state : null;
}

/** The dataset's own spelling of a gene, or the near matches when there is
 * none. */
async function resolveGene(
  id: string,
  gene: string,
  signal: AbortSignal,
): Promise<{ name: string } | { near: string[] }> {
  const path = isDatasetId(id)
    ? `/datasets/${encodeURIComponent(id)}/genes`
    : `/catalog/${encodeURIComponent(id)}/genes`;
  const page = await apiGet<Page<string>>(path, { q: gene, limit: 20 }, signal);
  const exact = page.items.find((g) => g.toLowerCase() === gene.toLowerCase());

  return exact ? { name: exact } : { near: page.items.slice(0, 8) };
}

/** Pick the viewer link to open: the one whose label matches `part`, else
 * the first that can show what was asked (cells first). */
function pickEntry(
  d: Detail,
  part: string | undefined,
): { entry: Entry | null; choices: Entry[] } {
  const all: Entry[] = d.entries?.length
    ? d.entries.filter((e) => e.viewerUrl)
    : d.viewerUrl
      ? [{ label: d.title ?? d.id, kind: d.kind ?? "", viewerUrl: d.viewerUrl }]
      : [];

  if (part) {
    const p = part.toLowerCase();
    const hit =
      all.find((e) => e.label.toLowerCase() === p) ??
      all.find((e) => e.label.toLowerCase().includes(p)) ??
      all.find((e) =>
        p.startsWith("mol")
          ? e.kind === "single_molecule"
          : p.startsWith("cell")
            ? e.kind === "single_cell"
            : false,
      );

    return { entry: hit ?? null, choices: all };
  }

  return {
    entry: all.find((e) => e.kind === "single_cell") ?? all[0] ?? null,
    choices: all,
  };
}

export function registerViewTools(server: McpServer) {
  server.registerTool(
    "meyes_open_view",
    {
      title: "Open a dataset in the viewer",
      description: `Open a dataset in the person's browser, already showing what was asked: genes coloured, a cell-type column chosen and some of its types picked out, 2D or 3D. Returns the link too (it can be shared; it reopens the same view).
Gene names are corrected to the dataset's own spelling; a gene it does not measure is an error that lists the near matches. Column and cell-type names are passed to the viewer as given (it falls back to its default column when the name is not one of the dataset's) — use the names meyes_get_dataset reports.
Cells viewer: one gene colours the cells; two genes show co-expression. Molecule viewer: each gene is a colour (up to 10).
Examples:
  - "show me Slc17a7 in ace-den-fog" -> id="ace-den-fog", genes=["Slc17a7"]
  - "molecules of Gad1 and Slc17a7 in section 609889HK01" -> id="ace-den-fog", part="609889HK01", genes=["Gad1","Slc17a7"]
  - "my dataset ds_… in 3D coloured by subclass" -> id="ds_…", column="subclass", view="3D"`,
      inputSchema: z.strictObject({
        id: z
          .string()
          .min(1)
          .max(200)
          .describe("Catalog id, BIL code, or app dataset id"),
        part: z
          .string()
          .max(200)
          .optional()
          .describe(
            'Which of the dataset\'s viewer links, by label (e.g. a section name) or "cells" / "molecules"; default: its cells',
          ),
        genes: z
          .array(z.string().min(1).max(50))
          .max(10)
          .default([])
          .describe("Genes to show"),
        column: z
          .string()
          .max(200)
          .optional()
          .describe("Cell-type (obs) column to colour by, cells viewer only"),
        celltypes: z
          .array(z.string().min(1).max(200))
          .max(50)
          .optional()
          .describe("Values of that column to pick out, cells viewer only"),
        view: z.enum(["2D", "3D"]).default("2D").describe("2D or 3D"),
        open: z
          .boolean()
          .default(true)
          .describe("Open it in the browser (false: only return the link)"),
      }),
      annotations: {
        readOnlyHint: true,
        destructiveHint: false,
        idempotentHint: true,
        openWorldHint: true,
      },
    },
    async (a, ctx) => {
      const signal = ctx.mcpReq.signal;
      const d = await apiGet<Detail>(
        isDatasetId(a.id)
          ? `/datasets/${encodeURIComponent(a.id)}`
          : `/catalog/${encodeURIComponent(a.id)}`,
        {},
        signal,
      );
      const { entry, choices } = pickEntry(d, a.part);
      const fail = (msg: string) => ({
        isError: true,
        content: [{ type: "text" as const, text: msg }],
      });

      if (!entry?.viewerUrl) {
        return fail(
          choices.length
            ? `No viewer link of ${a.id} matches "${a.part}". Its links: ${choices
                .slice(0, 15)
                .map((e) => `"${e.label}"`)
                .join(
                  ", ",
                )}${choices.length > 15 ? `, … (${choices.length} in all)` : ""}.`
            : `${a.id} has no viewer link yet (still processing, or nothing to show).`,
        );
      }
      const url = new URL(entry.viewerUrl);

      // Only ever open the site this server talks to.
      if (url.origin !== new URL(SITE).origin)
        return fail(`The viewer link is not on ${SITE}; not opened.`);
      const kind = viewerKind(url);

      if (kind === "none" && (a.genes.length || a.column || a.celltypes))
        return fail(
          `"${entry.label}" opens in a viewer that takes no gene or column settings (an image). Call again without them, or pick another part.`,
        );
      if (kind === "molecules" && (a.column || a.celltypes?.length))
        return fail(
          `"${entry.label}" is a molecule viewer: it has no cell-type column. Pass only genes, or pick the cells part.`,
        );
      if (kind === "cells" && a.genes.length > 2)
        return fail(
          "The cells viewer shows one gene, or two as co-expression. Pass at most two, or use a molecule part for more.",
        );

      const genes: string[] = [];

      for (const g of a.genes) {
        const r = await resolveGene(d.id ?? a.id, g, signal);

        if ("near" in r)
          return fail(
            `${a.id} does not measure "${g}".${r.near.length ? ` Near matches: ${r.near.join(", ")}.` : ""} Check with meyes_list_genes.`,
          );
        genes.push(r.name);
      }

      const state = viewState(kind, {
        genes,
        column: a.column,
        celltypes: a.celltypes,
        mode3d: a.view === "3D",
      });

      if (state) url.searchParams.set("v", toBase64Url(state));
      const link = url.toString();

      if (a.open) openBrowser(link);

      const shown = [
        genes.length && `genes ${genes.join(", ")}`,
        genes.length === 2 && kind === "cells" && "as co-expression",
        a.column && `column ${a.column}`,
        a.celltypes?.length && `cell types ${a.celltypes.join(", ")}`,
        a.view === "3D" && "3D",
      ].filter(Boolean);

      return text(
        [
          `${a.open ? "Opened" : "Link to"} ${d.title ?? a.id} — ${entry.label}${shown.length ? `, with ${shown.join(", ")}` : ""}:`,
          "",
          link,
          ...(a.column || a.celltypes?.length
            ? [
                "",
                "Column and cell-type names were passed as given; if the viewer shows another column, the name is not one of this dataset's.",
              ]
            : []),
          ...(a.open
            ? [
                "",
                "If no tab opened (no browser on this computer), give the person the link.",
              ]
            : []),
        ].join("\n"),
      );
    },
  );
}
