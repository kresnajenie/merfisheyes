import type { McpServer } from "@modelcontextprotocol/server";

import * as z from "zod/v4";

import { apiGet } from "../api.js";
import { openBrowser } from "../auth.js";
import { SITE } from "../config.js";
import { type Page, text } from "../format.js";
import { isDatasetId } from "./catalog.js";

/**
 * meyes_open_view: a link to a dataset in the viewer, already set up — genes
 * coloured, a cell-type column chosen, 2D or 3D, dot size and opacity, and
 * optionally a second dataset beside it in split screen. No bridge to a tab:
 * the viewers read their settings from the link (`v=` / `rv=`, base64url
 * JSON, lib/utils/url-visualization-state.ts; split screen from `split` /
 * `splitS3Url`, `splitType`, `sync`, lib/hooks/useSplitUrlParams.ts), so the
 * tool builds that link, after correcting gene names to the dataset's own
 * spelling (the viewer applies only an exact match). It returns the link;
 * it opens a browser tab only when asked, one per call.
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

/** The viewers whose pages read the split-screen parameters. */
const SPLIT_PAGES = ["/viewer", "/sm-viewer"];
const splitCapable = (url: URL) =>
  SPLIT_PAGES.some(
    (p) => url.pathname === p || url.pathname.startsWith(`${p}/`),
  );

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
  /** Dot size multiplier, 1 = the viewer's default. */
  dotSize?: number;
  /** Opacity of coloured cells, 0–1 (cells viewer). */
  opacity?: number;
  /** Opacity of cells not among the picked-out types, 0–1 (cells viewer). */
  unselectedOpacity?: number;
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
    if (s.dotSize !== undefined) state.sz = s.dotSize;
    const av: Record<string, number> = {};

    if (s.opacity !== undefined) {
      av.expressionAlphaMin = s.opacity;
      av.expressionAlphaMax = s.opacity;
    }
    if (s.unselectedOpacity !== undefined)
      av.greyedOutAlpha = s.unselectedOpacity;
    if (Object.keys(av).length) state.av = av;
  } else if (kind === "molecules") {
    if (s.genes.length)
      state.genes = s.genes.map((g, i) => [
        g,
        GENE_COLORS[i % GENE_COLORS.length],
        1,
        true,
      ]);
    if (s.dotSize !== undefined) state.gs = s.dotSize;
  } else return null;
  if (s.mode3d) state.vm = "3D";

  return Object.keys(state).length ? state : null;
}

/** Where a viewer link points the split screen's right panel: an app dataset
 * by id, or a store by its URL (the `url` parameter of a from-s3 link). */
export function splitParams(right: URL): Record<string, string> | null {
  const type = viewerKind(right) === "molecules" ? "sm" : "cell";

  if (!splitCapable(right)) return null;
  const s3 = right.searchParams.get("url");

  // useSplitUrlParams reads splitS3Url with decodeURIComponent after
  // URLSearchParams has decoded it once: it is written encoded twice.
  if (right.pathname.endsWith("/from-s3") && s3)
    return { splitS3Url: encodeURIComponent(s3), splitType: type };
  const id = right.pathname.split("/").filter(Boolean)[1];

  if (id && id !== "from-local" && id !== "from-s3")
    return { split: id, splitType: type };

  return null;
}

class ViewError extends Error {}

/** The dataset's own spelling of a gene, or a ViewError listing near
 * matches. */
async function resolveGene(
  id: string,
  gene: string,
  signal: AbortSignal,
): Promise<string> {
  const path = isDatasetId(id)
    ? `/datasets/${encodeURIComponent(id)}/genes`
    : `/catalog/${encodeURIComponent(id)}/genes`;
  const page = await apiGet<Page<string>>(path, { q: gene, limit: 20 }, signal);
  const exact = page.items.find((g) => g.toLowerCase() === gene.toLowerCase());

  if (exact) return exact;
  const near = page.items.slice(0, 8);

  throw new ViewError(
    `${id} does not measure "${gene}".${near.length ? ` Near matches: ${near.join(", ")}.` : ""} Check with meyes_list_genes.`,
  );
}

/** Pick the viewer link to open: the one whose label matches `part`, else
 * the first that can show what was asked (cells first). */
function pickEntry(
  d: Detail,
  id: string,
  part: string | undefined,
): Entry & { viewerUrl: string } {
  const all = (
    d.entries?.length
      ? d.entries
      : [
          {
            label: d.title ?? d.id,
            kind: d.kind ?? "",
            viewerUrl: d.viewerUrl ?? null,
          },
        ]
  ).filter((e): e is Entry & { viewerUrl: string } => !!e.viewerUrl);
  let hit: (typeof all)[number] | undefined;

  if (part) {
    const p = part.toLowerCase();

    hit =
      all.find((e) => e.label.toLowerCase() === p) ??
      all.find((e) => e.label.toLowerCase().includes(p)) ??
      all.find((e) =>
        p.startsWith("mol")
          ? e.kind === "single_molecule"
          : p.startsWith("cell")
            ? e.kind === "single_cell"
            : false,
      );
  } else hit = all.find((e) => e.kind === "single_cell") ?? all[0];
  if (hit) return hit;
  throw new ViewError(
    !all.length
      ? `${id} has no viewer link yet (still processing, or nothing to show).`
      : `No viewer link of ${id} matches "${part}". Its links: ${all
          .slice(0, 15)
          .map((e) => `"${e.label}"`)
          .join(", ")}${all.length > 15 ? `, … (${all.length} in all)` : ""}.`,
  );
}

const sideSchema = {
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
    .describe(
      "Genes to show: 1, or 2 as co-expression, in a cells viewer; up to 10 coloured in a molecule viewer",
    ),
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
};

type Side = {
  id: string;
  part?: string;
  genes: string[];
  column?: string;
  celltypes?: string[];
};

/** One side of the view: its dataset, viewer link and state. */
async function buildSide(
  s: Side,
  shared: Omit<ViewSettings, "genes" | "column" | "celltypes">,
  signal: AbortSignal,
) {
  const d = await apiGet<Detail>(
    isDatasetId(s.id)
      ? `/datasets/${encodeURIComponent(s.id)}`
      : `/catalog/${encodeURIComponent(s.id)}`,
    {},
    signal,
  );
  const entry = pickEntry(d, s.id, s.part);
  const url = new URL(entry.viewerUrl);

  // Only ever link to the site this server talks to.
  if (url.origin !== new URL(SITE).origin)
    throw new ViewError(`The viewer link is not on ${SITE}.`);
  const kind = viewerKind(url);

  if (kind === "none" && (s.genes.length || s.column || s.celltypes))
    throw new ViewError(
      `"${entry.label}" of ${s.id} opens in a viewer that takes no gene or column settings (an image). Leave them out, or pick another part.`,
    );
  if (kind === "molecules" && (s.column || s.celltypes?.length))
    throw new ViewError(
      `"${entry.label}" of ${s.id} is a molecule viewer: it has no cell-type column. Pass only genes, or pick the cells part.`,
    );
  if (kind === "cells" && s.genes.length > 2)
    throw new ViewError(
      `The cells viewer shows one gene, or two as co-expression. Pass at most two for ${s.id}, or use a molecule part (part="molecules") for more.`,
    );
  const genes: string[] = [];

  for (const g of s.genes) genes.push(await resolveGene(d.id, g, signal));

  return {
    title: d.title ?? s.id,
    entry,
    url,
    kind,
    genes,
    state: viewState(kind, {
      ...shared,
      genes,
      column: s.column,
      celltypes: s.celltypes,
    }),
  };
}

const describe = (
  side: Awaited<ReturnType<typeof buildSide>>,
  s: Side,
): string =>
  [
    `${side.title} (${s.id}) — ${side.entry.label}`,
    side.genes.length &&
      `genes ${side.genes.join(", ")}${side.genes.length === 2 && side.kind === "cells" ? " (co-expression)" : ""}`,
    s.column && `column ${s.column}`,
    s.celltypes?.length && `cell types ${s.celltypes.join(", ")}`,
  ]
    .filter(Boolean)
    .join(", ");

export function registerViewTools(server: McpServer) {
  server.registerTool(
    "meyes_open_view",
    {
      title: "Link to a dataset in the viewer",
      description: `A link to a dataset in the MERFISHEYES viewer, already showing what was asked: genes coloured, a cell-type column with some of its types picked out, 2D or 3D, dot size, opacity — and, with \`compare\`, a second dataset beside it in split screen (two slices, two ages, two conditions; camera synced with sync=true).
Returns the link; the link reopens the same view and can be shared. Set open=true only when the person asked to open it now — it opens one browser tab. Never call this several times to open several tabs: give links, or use split screen for two.
Gene names are corrected to the dataset's own spelling; a gene it does not measure is an error listing near matches. Column and cell-type names are passed as given (the viewer falls back to its default column when the name is not one of the dataset's) — use names meyes_get_dataset reports.
Split screen works for cell and molecule viewers (not SpatialData stores or images). Settings that are not per side (view, dot_size, opacity) apply to both.
Examples:
  - "show me Slc17a7 in ace-den-fog" -> id="ace-den-fog", genes=["Slc17a7"]
  - "molecules of Gad1 and Slc17a7 in section 609889HK01, bigger dots" -> id="ace-den-fog", part="609889HK01", genes=["Gad1","Slc17a7"], dot_size=2
  - "compare the P0 and P4 brains side by side, coloured by class" -> id="ace-let-law", column="class_name", compare={id:"ace-let-lax", column:"class_name"}, sync=true`,
      inputSchema: z.strictObject({
        ...sideSchema,
        compare: z
          .strictObject(sideSchema)
          .optional()
          .describe(
            "A second dataset (or another part of the same one) for the right panel of a split screen",
          ),
        sync: z
          .boolean()
          .default(true)
          .describe("Split screen: move both panels' cameras together"),
        view: z.enum(["2D", "3D"]).default("2D").describe("2D or 3D"),
        dot_size: z
          .number()
          .min(0.1)
          .max(10)
          .optional()
          .describe("Dot size multiplier: 1 is the default, 2 twice as big"),
        opacity: z
          .number()
          .min(0)
          .max(1)
          .optional()
          .describe("Opacity of coloured cells, 0–1 (cells viewer)"),
        unselected_opacity: z
          .number()
          .min(0)
          .max(1)
          .optional()
          .describe(
            "Opacity of cells outside the picked-out cell types, 0–1 (cells viewer)",
          ),
        open: z
          .boolean()
          .default(false)
          .describe(
            "Open it in a browser tab now; only when the person asked to open it",
          ),
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
      const shared = {
        mode3d: a.view === "3D",
        dotSize: a.dot_size,
        opacity: a.opacity,
        unselectedOpacity: a.unselected_opacity,
      };

      try {
        const left = await buildSide(a, shared, signal);
        const right = a.compare
          ? await buildSide(a.compare, shared, signal)
          : null;
        const url = left.url;

        if (left.state) url.searchParams.set("v", toBase64Url(left.state));
        if (right) {
          if (!splitCapable(url))
            throw new ViewError(
              `"${left.entry.label}" opens in a viewer without split screen (a SpatialData store or an image). Give links to the two views instead.`,
            );
          const split = splitParams(right.url);

          if (!split)
            throw new ViewError(
              `"${right.entry.label}" of ${a.compare!.id} cannot be the right panel of a split screen (a SpatialData store or an image). Give links to the two views instead.`,
            );
          for (const [k, v] of Object.entries(split))
            url.searchParams.set(k, v);
          if (a.sync) url.searchParams.set("sync", "1");
          if (right.state) url.searchParams.set("rv", toBase64Url(right.state));
        }
        const link = url.toString();

        if (a.open) openBrowser(link);

        const extras = [
          a.view === "3D" && "3D",
          a.dot_size !== undefined && `dot size ×${a.dot_size}`,
          a.opacity !== undefined && `opacity ${a.opacity}`,
          a.unselected_opacity !== undefined &&
            `other cells at opacity ${a.unselected_opacity}`,
          right && (a.sync ? "cameras synced" : "cameras independent"),
        ].filter(Boolean);

        return text(
          [
            right
              ? `Split screen — left: ${describe(left, a)}; right: ${describe(right, a.compare!)}${extras.length ? `; ${extras.join(", ")}` : ""}.`
              : `${describe(left, a)}${extras.length ? `; ${extras.join(", ")}` : ""}.`,
            "",
            link,
            "",
            a.open
              ? "Opened in a browser tab. If none opened (no browser on this computer), give the person the link."
              : "Give the person this link; it was not opened.",
            ...(a.column || a.celltypes?.length || a.compare?.column
              ? [
                  "Column and cell-type names were passed as given; if the viewer shows another column, the name is not one of this dataset's.",
                ]
              : []),
          ].join("\n"),
        );
      } catch (e) {
        if (e instanceof ViewError)
          return {
            isError: true,
            content: [{ type: "text" as const, text: e.message }],
          };
        throw e;
      }
    },
  );
}
