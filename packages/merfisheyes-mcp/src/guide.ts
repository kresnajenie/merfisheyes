import type { McpServer } from "@modelcontextprotocol/server";

import * as z from "zod/v4";

/**
 * The primer an agent reads before working with MERFISHEYES (meyes_guide):
 * what the tools are for, in what order, what to ask before acting, and the
 * traps that cost earlier agents their time (MCP.md, the 2026-11-24 dry runs
 * on a P28 mouse brain). Server instructions stay short (clients cut them);
 * this is the long form.
 */
export const GUIDE = `# Working with MERFISHEYES

MERFISHEYES is a viewer for spatial transcriptomics (cells, molecules, outlines, images). With this package the person's own data opens in the viewer on their computer — served from their machine, never uploaded — and you can read, change and picture the view, and run analysis steps on the data.

## Ask before you assume
The person's words rarely fix everything an analysis needs. Before running anything that takes more than a moment, ask (a few questions at once, each with your recommendation):
- Which dataset, and which sections/samples (multi-section files often lay sections side by side).
- Which labels define a region or population. Cell-type labels are not anatomy: "DG neurons" marks the dentate gyrus, ependymal/astrocyte labels line every ventricle, olfactory-bulb neuroblasts are not the SVZ.
- Distances (a niche radius in µm), and whether to measure them within each section.
- Whether sections/batches need correcting (not built in: say clusters made of one section may be batch effects).
- Where results go: a copy of their file. Ask which format — .zarr (a SpatialData store: the format MERFISHEYES reads best, opens at once, labels carry DE stats) or .h5ad (for scanpy). The original is never changed.
Say what you chose when you choose something yourself (a root cell, a cluster name).

## Order of work on local data
1. meyes_open_local(path) → give the person the link (it opens the viewer on their computer; an .h5ad or raw instrument output is converted once — the first time needs a one-time Python setup: ask, then meyes_setup_converter).
2. meyes_view_state → columns, the colour column's values with cell counts. For long lists: column=…, containing=… ("Epen", "IMN", "DG").
3. meyes_profile → the panel and the data: genes=[markers you expect] says which the panel measures (panels are often 300–500 genes: name the missing markers), breakdown={column, by: section column} shows where cells are. Coordinates' spacing tells you whether they are microns.
4. Count before you conclude: a type with 2 cells in a section is not a niche there; pick another section or dataset and say why.
5. Analysis: meyes_niche (cells within radius of anchor labels, per section) → meyes_recluster (that subset: UMAP + Leiden at a few resolutions; each cluster's top genes and main existing labels) → meyes_markers (per cluster; genes=[literature markers], heatmap=true) → name clusters from evidence, "?" when unsure → meyes_save_analysis (the format they chose).
6. Show and check: meyes_view_set (colour by the new labels, hide "" to show the niche alone, embedding="umap_<name>") and meyes_view_snapshot (focus=[values] frames one niche; target="embedding" pictures the UMAP). Look at every snapshot before describing it.

## A raw MERFISH run (the lab pipeline's .npz output: final_spots/ + final_mats/)
1. Check it before anything else: meyes_raw_check(path) → show the person the plots and read out the flags with what they mean (a field with a high blank fraction decoded poorly; a gene no commoner than the blanks is not measured; a set far from the others' depth may be another library). genes=[…] draws every molecule of each in one set (set=…).
2. Ask before converting: the whole run or a few sets first (sets=[…], minutes); an .h5ad of the same run whose labels to add (labels_from; look for one in the run folder); where the output goes (output) — it takes about 13 GB per billion molecules and 1–3 minutes per 100 million (meyes_raw_check gives conversion_estimate).
3. meyes_open_local(path, sets, labels_from, output) converts it as a job: tell the person, then meyes_job(id) until it gives the link. Sets are sections side by side (colour by "section"); "field" names each cell's field (hide a bad one); n_blank and gfp_count are per cell.
4. Then as any dataset: meyes_view_state, meyes_view_snapshot, the analysis tools.

Long steps run as jobs: a tool answers within about a minute; if it hands back a job id, tell the person it is running and call meyes_job(id).

## Knowledge you need but do not have
Use web search and the literature (cite them) for marker genes of a lineage, the names in a taxonomy, what a region is made of. Then check them against the data (meyes_profile genes=…, meyes_markers genes=…). Allen/ABC atlas labels (MapMyCells) look like "045 OB-STR-CTX Inh IMN": number, region, type; "IMN" = immature neurons; levels class > subclass > supertype > cluster. Markers are not drivers: say "marked by", not "driven by", unless a source shows it.

## Your own computations (where you can run Python)
Prefer the analysis tools: they run in MERFISHEYES's own Python, need no environment of yours, and survive long runs. If you compute something yourself, put it in the view with meyes_view_add_column (a CSV keyed by cell id; an embedding as embedding={name, x, y}).

## Questions about the library
"Which/how many datasets …" → meyes_find_datasets (organism, tissue with everything below it, disease, technology, age, genes measured/expressed) — one call answers count, breakdown and rows; meyes_index_values spells values. Report its notes (rows left out, genes it could not check).

## Traps
- Several sections in one file sit side by side in one coordinate frame: distances across them mean nothing — use within=<section column>.
- A view command before the person opens the link waits for the tab; if it never opens, ask them to open it.
- Coloured by a gene, the scene still draws only the picked-out values; show_all colours it everywhere.
- A store with molecules may draw some genes' molecules over the cells by default, hiding them in pictures: molecules={show:[]} clears them.
- meyes_view_set answers briefly; meyes_view_state has the detail.
- Results live in the workspace (and the open view) until meyes_save_analysis writes the copy.
`;

/** meyes_guide (any client) and the same text as a resource (clients that
 * read resources themselves, e.g. Claude Code). */
export function registerGuide(server: McpServer) {
  server.registerTool(
    "meyes_guide",
    {
      title: "How to work with MERFISHEYES",
      description:
        "Read this once before working with the person's data or the viewer: the order of the tools, what to ask the person before acting, where to look things up, and the traps earlier agents fell into. Short; saves many calls.",
      inputSchema: z.strictObject({}),
      annotations: {
        readOnlyHint: true,
        destructiveHint: false,
        idempotentHint: true,
        openWorldHint: false,
      },
    },
    async () => ({ content: [{ type: "text" as const, text: GUIDE }] }),
  );
  server.registerResource(
    "guide",
    "merfisheyes://guide",
    {
      title: "How to work with MERFISHEYES",
      description: "The primer meyes_guide returns.",
      mimeType: "text/markdown",
    },
    async (uri) => ({
      contents: [{ uri: uri.href, mimeType: "text/markdown", text: GUIDE }],
    }),
  );
}
