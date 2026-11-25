/**
 * The server instructions — the manual an agent reads when it connects
 * (MCP.md pillar E3). One text for both servers: `local` adds what only the
 * package on the person's computer can do (sign in there, open a browser
 * tab, serve their own folders).
 */
export function instructions(site: string, local: boolean): string {
  return [
    `MERFISHEYES (${site}) is a web viewer for spatial transcriptomics: cells, molecules, cell outlines and microscopy images, at up to billions of molecules. These tools find datasets and read the person's own projects; the site does the heavy work.`,
    "",
    local
      ? "- Public catalog tools (meyes_find_datasets, meyes_index_values, meyes_search_datasets, meyes_list_facets, meyes_get_dataset, meyes_list_genes) need no sign-in. The person's own data needs meyes_sign_in once; ask only for the scopes the task needs."
      : "- The person signed in to MERFISHEYES when they added this connector; meyes_whoami says as whom and which projects this agent may reach.",
    "- Questions about the library — how many datasets, which ones, of what organism, tissue, disease, technology, age, genotype, measuring or expressing which genes — go to meyes_find_datasets: one call answers the count, a breakdown and the first rows, over every public dataset (MERFISHEYES's own and SpatialData-DB's). Report its notes (rows left out for lacking a field, genes it could not check) with the answer. Use meyes_index_values to spell a value; use meyes_search_datasets only for Explore's cards.",
    "- Prefer what MERFISHEYES already records over computing it yourself: meyes_get_dataset returns the counts, cell-type levels, DE columns and QC each store's own files hold. Never fetch raw per-cell data to answer something a tool answers.",
    local
      ? "- To show the person something, call meyes_open_view: it returns a link to the viewer already set up (genes, column, cell types, 2D/3D, dot size, opacity). Give them the link. Pass open=true only when they ask to open it now, and never open more than one tab per request: for several views give links, and for two views to compare use one split screen (compare=…). Do not describe a scene you cannot see."
      : "- To show the person something, call meyes_open_view: it returns a link to the viewer already set up (genes, column, cell types, 2D/3D, dot size, opacity). Give them the link; for two views to compare, use one split screen (compare=…). Do not describe a scene you cannot see.",
    ...(local
      ? [
          "- Data on the person's own computer (a .zarr store, a MERFISHEYES folder, an OME-Zarr image) opens with meyes_open_local: it is served from their machine to the viewer, never uploaded. Use it, not an upload, when they ask to see their local data.",
          "- Once the person has the view from meyes_open_local open, drive it: meyes_view_state (what it shows, column and value names), meyes_view_set (colour by a gene or column, show only / hide values such as a bad FOV, dot size, 2D/3D), meyes_view_snapshot (a PNG with legend, saved to a file and shown to you — look at it before describing the view, and use it for figures). Results you compute (pseudotime, labels, scores) go into the view with meyes_view_add_column from a CSV keyed by cell id; then colour by them, hide values and snapshot as with the dataset's own columns.",
          "- When the person asks to sign in, or asks for their own projects or datasets while signed out, call meyes_sign_in: it opens the sign-in page and waits for them.",
        ]
      : [
          "- Data on the person's own computer cannot be reached from here: the @merfisheyes/mcp package (npx -y @merfisheyes/mcp@latest view <folder>) opens it in the viewer without uploading it.",
        ]),
    "- Lists are paged: use limit/offset and next_offset rather than repeating a call. Calls are rate-limited; a 429 says how long to wait.",
    "- Not every column with DE stats is a cell type (some stores carry QC columns too); group by annotation columns.",
  ].join("\n");
}
