import type { McpServer } from "@modelcontextprotocol/server";

import { mkdirSync, readFileSync, writeFileSync } from "node:fs";
import { dirname, join, resolve } from "node:path";

import * as z from "zod/v4";

import { text } from "../format.ts";
import { cacheRoot } from "../local/convert.ts";
import { local } from "../local/serve.ts";
import { idColumn, parseTable } from "../local/table.ts";

/**
 * The live session (MCP.md pillar D3): the agent reads, changes and
 * pictures the view open in the person's browser — the tab
 * meyes_open_local's link opened, served by this package. The page runs the
 * commands (lib/agent/view-commands.ts in the viewer).
 */

const fail = (e: unknown) => ({
  isError: true,
  content: [{ type: "text" as const, text: (e as Error).message }],
});

const ACT = {
  readOnlyHint: false,
  destructiveHint: false,
  idempotentHint: true,
  openWorldHint: false,
};

export function registerSessionTools(server: McpServer) {
  server.registerTool(
    "meyes_view_state",
    {
      title: "What the open view shows",
      description: `What the MERFISHEYES view open in the person's browser shows right now: the dataset, cells and genes, what it is coloured by, what is drawn (counted from the cells: visible_cells), dot size, 2D/3D, the embeddings it has (UMAPs …) and which is shown, every column with its type, and the colour column's values with cells per value.
With column, that column's values and their cell counts (it is loaded if need be) — filter with containing ("Epen" → every ependymal type), page with limit/offset. Use it to learn exact names before meyes_view_set, and to count before concluding (a cell type with 2 cells is not a niche).
Needs a view opened with meyes_open_local and still open in a tab.`,
      inputSchema: z.strictObject({
        column: z
          .string()
          .min(1)
          .max(200)
          .optional()
          .describe("List this column's values with their cell counts"),
        containing: z
          .string()
          .min(1)
          .max(100)
          .optional()
          .describe("Only values containing this text, case-insensitive"),
        limit: z.number().int().min(1).max(500).default(100),
        offset: z.number().int().min(0).default(0),
      }),
      annotations: { ...ACT, readOnlyHint: true },
    },
    async (a) => {
      try {
        return text(
          JSON.stringify(await local.command("state", a, 120_000), null, 2),
        );
      } catch (e) {
        return fail(e);
      }
    },
  );

  server.registerTool(
    "meyes_view_set",
    {
      title: "Change the open view",
      description: `Change the MERFISHEYES view open in the person's browser — cells: colour by a gene or by a column (cell type, cluster, pseudotime, fov, …), show only some values of the current categorical column or hide some (e.g. a bad FOV), show all again, dot size, 2D or 3D, reset the camera, zoom to what is kept (fit); show an embedding (a UMAP) in the embedding panel. Answers briefly with what is coloured and drawn afterwards (meyes_view_state for the rest).
A gene while values are picked out is drawn on those cells only; show_all colours it everywhere.
Molecules (the molecule viewer, or molecules drawn over cells): molecules={show:["Gad1","Slc17a7"], colors:{...}, opacity, size}.
Names must be the dataset's own (meyes_view_state lists them, with column=… and containing=… for long lists); gene names are matched without regard to case.
Examples: "colour by Slc17a7" -> gene="Slc17a7"; "colour by cell type and hide fov 5" -> column="fov", hide=["5"] then column="celltype"; "only radial glia and immature neurons" -> column="celltype", show_only=["radial glia","immature"].`,
      inputSchema: z.strictObject({
        gene: z
          .string()
          .min(1)
          .max(100)
          .optional()
          .describe("Colour by this gene"),
        column: z
          .string()
          .min(1)
          .max(200)
          .optional()
          .describe(
            "Colour by this column (and the column show_only/hide refer to)",
          ),
        show_only: z
          .array(z.string().max(200))
          .max(200)
          .optional()
          .describe("Show only these values of the column; hide the rest"),
        hide: z
          .array(z.string().max(200))
          .max(200)
          .optional()
          .describe("Hide these values of the column"),
        show_all: z.boolean().optional().describe("Show every value again"),
        treat_as: z
          .enum(["categorical", "numerical"])
          .optional()
          .describe(
            "Treat the column as categories or as numbers (the viewer's switch); a numeric column is treated as categories automatically when you show or hide its values",
          ),
        dot_size: z
          .number()
          .min(0.1)
          .max(10)
          .optional()
          .describe("Dot size multiplier: 1 is the default"),
        view: z.enum(["2D", "3D"]).optional(),
        reset_camera: z
          .boolean()
          .optional()
          .describe("Fit the data in view again"),
        molecules: z
          .strictObject({
            show: z
              .array(z.string().min(1).max(100))
              .max(30)
              .optional()
              .describe(
                "Exactly these genes' molecules drawn (others removed)",
              ),
            add: z.array(z.string().min(1).max(100)).max(30).optional(),
            remove: z.array(z.string().min(1).max(100)).max(30).optional(),
            colors: z
              .record(z.string(), z.string().regex(/^#[0-9a-fA-F]{6}$/))
              .optional()
              .describe('Colours per gene, e.g. {"Gad1": "#ff3da5"}'),
            opacity: z.number().min(0.01).max(1).optional(),
            size: z
              .number()
              .min(0.1)
              .max(10)
              .optional()
              .describe("Molecule size multiplier"),
          })
          .optional()
          .describe(
            "The molecules on show — in the molecule viewer, or molecules drawn over the cells",
          ),
        fit: z
          .boolean()
          .optional()
          .describe(
            "Zoom the view to the cells drawn in colour (after show_only / hide, the ones kept)",
          ),
        embedding: z
          .string()
          .max(200)
          .optional()
          .describe(
            'Show this embedding in the embedding panel, e.g. "umap" or "umap_DG" (meyes_view_state lists them); "" closes the panel',
          ),
      }),
      annotations: ACT,
    },
    async (a) => {
      try {
        return text(
          JSON.stringify(
            await local.command(
              "set",
              a,
              120_000, // the tab may still be loading its data
            ),
            null,
            2,
          ),
        );
      } catch (e) {
        return fail(e);
      }
    },
  );

  server.registerTool(
    "meyes_view_snapshot",
    {
      title: "Picture of the open view",
      description: `A PNG of the MERFISHEYES view open in the person's browser, framed on the cells drawn (or on focus values of the colour column — one niche or cluster —, or exactly as the person sees it, frame=current), with a legend beside it (the visible categories and their colours, or the gene and its colour scale). target="embedding" pictures the embedding panel (show one first with meyes_view_set embedding=…). Returned as an image you can look at, and saved to a file — save_to (a .png path, e.g. in the folder of a slide deck being made) or MERFISHEYES's snapshots folder. Look at it before describing the view; use it for figures.`,
      inputSchema: z.strictObject({
        save_to: z
          .string()
          .max(4096)
          .optional()
          .describe(
            "Where to save the PNG; default: MERFISHEYES's snapshots folder",
          ),
        width: z
          .number()
          .int()
          .min(400)
          .max(3000)
          .optional()
          .describe(
            "Width of the scene in pixels (the legend is added beside it)",
          ),
        frame: z
          .enum(["visible", "current"])
          .default("visible")
          .describe(
            "visible: framed on the cells drawn in colour (the person's camera is put back after); current: exactly the person's view",
          ),
        focus: z
          .array(z.string().max(200))
          .max(200)
          .optional()
          .describe(
            "Frame on these values of the colour column (a niche, a cluster), ignoring stray cells",
          ),
        target: z
          .enum(["scene", "embedding"])
          .default("scene")
          .describe(
            "scene: the cells in space; embedding: the embedding panel",
          ),
      }),
      annotations: ACT,
    },
    async ({ save_to, width, frame, focus, target }) => {
      try {
        const shot = await local.command<{
          png_base64: string;
          width: number;
          height: number;
          legend?: { value: string; color: string }[];
          gene?: string;
          column?: string;
          molecules?: { gene: string; color: string }[];
          embedding?: string;
        }>("snapshot", { width, frame, focus, target }, 120_000);
        const file = save_to
          ? resolve(save_to)
          : join(
              dirname(cacheRoot()),
              "snapshots",
              `snapshot-${new Date().toISOString().replace(/[:.]/g, "-")}.png`,
            );

        mkdirSync(dirname(file), { recursive: true });
        writeFileSync(file, Buffer.from(shot.png_base64, "base64"));

        return {
          content: [
            {
              type: "image" as const,
              data: shot.png_base64,
              mimeType: "image/png",
            },
            {
              type: "text" as const,
              text: [
                `Saved ${shot.width}×${shot.height} PNG to ${file}.`,
                shot.embedding ? `Embedding ${shot.embedding}.` : "",
                shot.gene
                  ? `Coloured by gene ${shot.gene}.`
                  : shot.column
                    ? `Coloured by ${shot.column}: ${(shot.legend ?? []).map((l) => l.value).join(", ")}.`
                    : "",
                shot.molecules?.length
                  ? `Molecules: ${shot.molecules.map((m) => m.gene).join(", ")}.`
                  : "",
              ].join("\n"),
            },
          ],
        };
      } catch (e) {
        return fail(e);
      }
    },
  );
}

export function registerAddColumnTool(server: McpServer) {
  server.registerTool(
    "meyes_view_add_column",
    {
      title: "Add your results to the open view",
      description: `Add columns you computed (pseudotime, new clusters, labels such as "radial glia" / "immature", scores) to the MERFISHEYES view open in the person's browser, without changing their data: they then behave like the dataset's own columns — colour by them, show only or hide values, snapshot.
Write the results to a CSV or TSV first: one row per cell, a header, an id column that matches the dataset's cell ids (pandas' df.to_csv() writes the index as the first column, which is used) and one column per result. Without an id column, give exactly one row per cell in the dataset's cell order. Cells missing from the file get an empty value.
Example: df[["pseudotime","fate"]].to_csv("/tmp/results.csv") then file="/tmp/results.csv", color_by=true.
An embedding you computed (a UMAP of one niche): two numeric columns, given as embedding={name:"umap_DG", x:"umap1", y:"umap2"}; cells without a row (outside the niche) are left out of it. It then shows in the embedding panel (meyes_view_set embedding=…).`,
      inputSchema: z.strictObject({
        file: z.string().min(1).max(4096).describe("The CSV or TSV"),
        id_column: z
          .string()
          .max(200)
          .optional()
          .describe(
            "The column with cell ids (default: an unnamed first column, or cell_id / obs_names / barcode / index)",
          ),
        columns: z
          .array(z.string().max(200))
          .max(50)
          .optional()
          .describe("Which columns to add (default: all but the id column)"),
        color_by: z
          .boolean()
          .default(false)
          .describe("Colour by the first added column"),
        replace: z
          .boolean()
          .default(false)
          .describe("Replace a column of the same name already in the view"),
        embedding: z
          .strictObject({
            name: z.string().min(1).max(100),
            x: z.string().min(1).max(200),
            y: z.string().min(1).max(200),
          })
          .optional()
          .describe(
            "Add the x and y columns as an embedding (shown in the embedding panel) instead of as columns",
          ),
      }),
      annotations: ACT,
    },
    async ({ file, id_column, columns, color_by, replace, embedding }) => {
      try {
        const path = resolve(file);
        const { header, rows } = parseTable(readFileSync(path, "utf8"));
        const id = idColumn(header, id_column);

        if (embedding) {
          const at = (name: string) => {
            const i = header.indexOf(name);

            if (i < 0)
              throw new Error(
                `No column "${name}" in ${path}; its columns: ${header.join(", ")}.`,
              );

            return i;
          };
          const xi = at(embedding.x);
          const yi = at(embedding.y);
          const coords = rows.map((r) => {
            const x = Number(r[xi]);
            const y = Number(r[yi]);

            return r[xi] === "" ||
              r[yi] === "" ||
              !Number.isFinite(x) ||
              !Number.isFinite(y)
              ? null
              : ([x, y] as [number, number]);
          });
          const res = await local.command<{
            embedding: string;
            cells: number;
            covered: number;
          }>(
            "add_embedding",
            {
              name: embedding.name,
              ids: id === null ? undefined : rows.map((r) => r[id] ?? ""),
              coords,
              replace,
              show: true,
            },
            120_000,
          );

          return text(
            `Added embedding ${res.embedding} to the open view: ${res.covered} of ${res.cells} cells have coordinates; it is shown in the embedding panel (meyes_view_snapshot target="embedding" pictures it).`,
          );
        }
        const wanted = (columns ?? header.filter((_, i) => i !== id)).map(
          (name) => {
            const i = header.indexOf(name);

            if (i < 0)
              throw new Error(
                `No column "${name}" in ${path}; its columns: ${header.join(", ")}.`,
              );

            return { name, i };
          },
        );

        if (!wanted.length)
          throw new Error(`${path} has no column to add besides the ids.`);
        const result = await local.command<{
          added: { name: string; type: string; n_values: number }[];
          cells: number;
          matched: number;
          unmatched: number;
          cells_without_a_value: number;
        }>(
          "add_column",
          {
            ids: id === null ? undefined : rows.map((r) => r[id] ?? ""),
            columns: wanted.map(({ name, i }) => ({
              name,
              values: rows.map((r) => r[i] ?? ""),
            })),
            color_by,
            replace,
          },
          120_000,
        );

        return text(
          [
            `Added to the open view: ${result.added.map((c) => `${c.name} (${c.type}, ${c.n_values} distinct values)`).join("; ")}.`,
            id === null
              ? `Values taken in cell order (${result.cells} cells).`
              : `Matched ${result.matched} of ${rows.length} rows to the dataset's ${result.cells} cells by "${header[id] || "the first column"}"${result.unmatched ? `; ${result.unmatched} ids matched no cell` : ""}${result.cells_without_a_value ? `; ${result.cells_without_a_value} cells have no value` : ""}.`,
            color_by
              ? `Coloured by ${result.added[0].name}.`
              : "Colour by them with meyes_view_set.",
          ].join("\n"),
        );
      } catch (e) {
        return fail(e);
      }
    },
  );
}
