# Plan: an agent that drives MERFISHEYES (MCP)

Let a person ask an agent — Claude Code first, then any MCP client — to run
an analysis on a dataset, and have the agent do it **through a MERFISHEYES
browser tab**: open datasets, colour, select, split panes, and show the
results — plots and statistics — in a panel beside the scene. Public datasets
and the person's own projects are both in reach, and the agent may do
anything the person could do on the site, within computational reason.

This file is the brief for every session working on it: the facts, the
decisions taken, the steps with their gates. **Status (2026-11-14): steps
1–4 done and released; 0.1.0 on npm; D1a built, not committed. After the
meeting with Tim Treis: metadata has its own plan
([METADATA.md](METADATA.md)), and the site itself becomes the MCP server
(pillar E).** Keep it current: a step that changes a fact or a decision edits this file in the
same commit. Branch: `feature/mcp` (from `bil/pipeline-dry-run`). The
contract the steps build is [docs/specs/MCP-PROTOCOL.md](../specs/MCP-PROTOCOL.md).

## 0. The pillars

**Where 0.1 stands.** `@merfisheyes/mcp` 0.1.0 is on npm (2026-11-13):
an agent can search the catalog, read a person's projects, and hand them
viewer links set up as asked (genes, columns, 2D/3D, split screen). That
is the site used the way a person uses it, faster. It is not yet what a
researcher would install it for.

**The two goals (2026-11-13, after a demo with a collaborator).**

1. **Research the library.** One question across everything MERFISHEYES
   holds — *"Alzheimer's is associated with gene X: compare every mouse
   brain sample aged 1 month to 1 year that has X"* — answered by a query
   over standardised metadata and per-dataset statistics on the server,
   returning a collection the person can compare, view and keep. Never by
   an agent walking S3 paths and computing on files.
2. **MERFISHEYES on your own data.** A researcher analysing data on their
   own machine — in Python, R, a notebook, with an agent — opens it in the
   MERFISHEYES viewer on `localhost`: nothing uploaded, the same viewer,
   driven by the agent as part of the analysis (show the spatial
   projection of cells and molecules without the bad FOV, colour by the
   pseudotime it just computed, capture the figures for the deck). Not a
   Python package: an npm package and a browser.

**What the ecosystem already has.** SpatialData-DB (Tim Treis, Theis lab;
`theislab/spatialdata-db-client`) publishes, per release: a catalog in a
versioned schema (`_catalog_schema.py`, schema 3) with ontology ids and
hierarchies (`organism_id`, `tissue_id` + `tissue_system`/`tissue_organ`,
`disease_id` + `disease_class`/`disease_group`), `development_stage`,
`timepoint_days`, `genotype`, `treatment`, `condition`, a technology
vocabulary, `tier`, `review_status`; a **gene index** (`gene_index.parquet`:
one row per dataset × gene, `total_counts`, `fraction_obs_detected`); a
feature index; cohort "primers" (shared embeddings, UMAPs, QC per cohort,
e.g. `mouse-brain-imaging`); and **its own MCP server** (Python: query by
facets and ranges, datasets expressing genes, describe, cite, download a
cohort). Goal 1 is that, for MERFISHEYES's library, plus what only
MERFISHEYES does: show it. So the metadata model is theirs unless the
conversation with Tim says otherwise (D20), and the two MCP servers are
meant to sit side by side in one agent: theirs to fetch and analyse a
cohort in Python, ours to find across both libraries and to look.

**After the meeting with Tim Treis (2026-11-14).** SpatialData-DB is the
gold standard for what the site and the MCP should be; what MERFISHEYES
adds is the viewer. Three things follow:

- **Metadata** — SDDB's schema and ontologies, validated with bionty,
  written into every store, recorded in Lamin with verification, and
  published as a query index of all public datasets that filters and
  agents start from. That is a plan of its own:
  [METADATA.md](METADATA.md). Pillar B queries what it builds.
- **The MCP lives on the site.** Tim's model: no install — an agent adds
  the site's `/mcp` address as a connector, the server's own
  instructions teach it the site, and long work runs asynchronously with
  the agent checking back. (Their published client still ships a stdio
  server; the hosted form is what Tim described.) So the site serves the
  same tools at `www.merfisheyes.com/mcp` (pillar E), and the npm package
  becomes the **local companion**: what needs the person's own machine —
  files on disk, the local viewer (pillar D).
- **Nothing changes for local use** except that a store carries its own
  metadata, so a store opened from disk shows it.

### Pillar A — Publish 0.1 — done

`@merfisheyes/mcp@0.1.0` on npm 2026-11-13; install is one line
(`claude mcp add --scope user merfisheyes -- npx -y
@merfisheyes/mcp@latest`). Left: PR #210 (README); 0.1.1 says "needs Node
20.3" instead of failing silently on older Node; the release
`develop` → `main` (molecule opacity in links, account page wording);
publishing from GitHub Actions later.

### Pillar B — Research the library (goal 1)

Everything here runs on the server over the **query index** that
[METADATA.md](METADATA.md) builds (steps 6–7: SDDB's catalog format for
every public dataset, plus a gene index); a tool call is one query, its
answer a table, a count and a breakdown — no file reads. Items 1–2 below
are METADATA.md's to build; 3–6 are this plan's.

1. **One metadata model for every dataset** — SDDB's catalog schema as
   columns (not a JSON blob) on the catalog: BIL, community and personal
   datasets carry the same fields as SDDB rows. The SDDB sync imports all
   of them (today it maps a subset and keeps the rest in `metadata`).
   Region queries follow the UBERON hierarchy ("brain" finds hippocampus);
   age queries use `timepoint_days`.
2. **One gene index** — dataset × gene, total counts and fraction of cells
   detected, from SDDB's `gene_index.parquet` for theirs and from our
   pipeline's statistics for ours (pillar C). "Has gene X" means measured;
   "expresses X" means a fraction above a threshold, and the tool says
   which it used.
3. **`meyes_find`** — structured filters (organism, region with
   descendants, age range in days, sex, disease and its class, genotype,
   treatment, technology, genes measured / expressed with a minimum
   fraction, tier) plus text; answers the count, a breakdown by the
   facets that matter, and the first rows. Replaces paging through search
   to count.
4. **Collections** — a found set saved with its query (so it can be
   re-run as the library grows), opened on the site, shared, and
   compared. Stored as a project for now; a collection page later.
5. **Compare a collection** — one table: for each dataset, gene X's
   fraction and counts, cells, age, genotype, region; per cell type when
   the dataset records it. Views: a split screen for two, a grid of links
   for more.
6. **Per-dataset statistics served** (the analysis half of step 6, now
   inside this pillar): structure (one slice or combined, slices, obs
   columns and values), gene stats, composition, gene by cell type,
   markers — read from the files the pipeline writes; "not recorded for
   this dataset" when they are missing, never a guess.

Gate: the Alzheimer's-style question above answered in one `meyes_find`
plus one compare call on the live site, with the numbers checked against
the stores; no S3 read by the agent.

### Pillar C — Metadata and statistics: [METADATA.md](METADATA.md)

Its own plan since 2026-11-14: SDDB's schema (M1), metadata in the store,
in Lamin and in the index (M2), bionty validation and verification (M3),
the query index (M4), the gene index (M5), the backfill of what exists,
and a site built to filter. This plan keeps one part: **the sweep** — a
scripted run of the read tools over the whole catalog that fills §10 and
checks every row against the schema, instead of waiting for an agent to
trip over each error.

### Pillar D — MERFISHEYES on your own data (goal 2)

The package grows a local mode. Rule 1 holds: **one viewer** — the local
viewer is the site's viewer code, never a fork.

1. **Serve and open (prototype, days).** `merfisheyes-mcp view <path>`
   and a tool `meyes_local_open(path)`: the package serves the folder on
   `127.0.0.1` with CORS and byte ranges (as `e2e/helpers/fixture-server.ts`
   does for the tests) and opens the site's viewer on that URL. The
   browser reads the files from the person's machine; nothing is
   uploaded. Opens what the viewer opens by URL: SpatialData Zarr v3
   stores, MERFISHEYES folders, OME-Zarr images. Known limits, to measure:
   needs the internet for the page; Chrome asks once to let the site
   reach `localhost` (local network access); Safari to be tried; an
   `.h5ad` needs a conversion first (find the shortest one from an
   analysis: `adata.write_zarr`, or `prepare.py`).

   **D1a done (2026-11-13; not yet released).** Site: `/open?url=` picks
   the viewer for any URL with the logic the "Load from URL" dialog had
   (`lib/viewer/viewer-for-url.ts`, now its one implementation); a URL on
   this computer (`lib/utils/local-url.ts`) is never sent to the site —
   the five `by-url` lookups go through `lib/datasets/lookup-by-url.ts`,
   which skips it, the claim banner hides, and analytics drops it
   (`components/site-analytics.tsx`). Package: `src/local/serve.ts`
   (127.0.0.1 only, a random token per folder, inside the folder only —
   symlinks included —, GET/HEAD with byte ranges, CORS for `SITE` only,
   the Host header checked, Chrome's private-network preflight answered),
   `merfisheyes-mcp view <folder>`, the tool `meyes_open_local`. Gate:
   package unit tests (serving and every refusal); `e2e/tests/viewer/open-url.spec.ts`
   (`/open` routes a store, a cells folder and a molecules folder served
   on 127.0.0.1, passes `v=` on, makes no `by-url` request, shows no claim
   banner); the package's server and `/open` together on a local site, all
   three in the right viewer. Measured on the live site before `/open`
   existed: the viewer reads a store from 127.0.0.1 once Chrome's
   local-network-access permission is granted; without it Chrome blocks
   the request ("Permission was denied … loopback address space").
   Left for D1b: `.h5ad` and raw instrument output; Safari and Firefox
   tried on a Mac.
2. **A standalone viewer.** A second build target of the same viewer
   components — static files, no server, no account — shipped as its own
   npm package and served by the local mode next to the data: offline,
   no permission prompt, same origin.
3. **A live session.** Because the package serves the page, page and
   agent share an origin: a WebSocket without pairing codes. Tools:
   `meyes_view_set` (colour by gene or column, hide categories — the bad
   FOV —, filter, camera, split), `meyes_view_state` (what is shown and
   selected), `meyes_snapshot` (a PNG at a chosen size with its legend —
   for figures and decks). This is the bridge of step 5, built where it
   is easiest first.
4. **The agent's results in the view.** Columns the agent computed —
   pseudotime, clusters, labels — pushed by cell id (a CSV or parquet
   beside the data) and shown without rewriting the store; cell sets
   highlighted or hidden.
5. **Fit an analysis.** The reference workflow (a collaborator's, 2026-11-13):
   define radial glia vs immature neurons from pseudotime, UMAP and a
   molecule gradient; output UMAP, spatial projection of cells and
   molecules, dot plots and pseudotime, the bad FOV removed, gradient
   genes plotted, all in a `.pptx`. The agent does the statistics in
   Python; MERFISHEYES draws the spatial views and hands back snapshots
   for the deck.

Gate: that workflow end to end on a real dataset, the spatial figures in
the deck coming from `meyes_snapshot`, nothing leaving the machine.

### Pillar E — The site is the MCP server (a connector, no install)

1. **`www.merfisheyes.com/mcp`** — a remote MCP server (Streamable HTTP,
   stateless, so it runs on Vercel's functions) serving the same tools as
   the package. Signing in is the OAuth server we already run (steps
   2–3): the connector discovers it from the site, registers through our
   client metadata document or registration, and the person approves on
   the consent page they know; the token's resource becomes the `/mcp`
   address. Works in claude.ai, Claude Desktop and Claude Code connectors
   and any client that takes a URL — nothing installed.
2. **One implementation of every tool.** The package and `/mcp` must not
   grow two copies of the tool descriptions and formatting: the tool
   definitions move to one module that both import, each with its own
   way of reaching the data (the package over the agent API; `/mcp` in
   process).
3. **Taught by the site.** The server instructions are the manual: when
   to use which tool, what the index holds, what the viewer can show.
   Versioned with the site, so every connector learns a change at once.
4. **Asynchronous work.** A sweep, a collection compared across hundreds
   of datasets, a processing job: started by one call that returns at
   once, with the agent checking back (MCP's tasks, or a job id and a
   status tool), and progress the person can see on the site.
5. **A tab on the site, driven from the connector.** The live session of
   pillar D3, relayed through the site to a tab the person has open (the
   old step 5 bridge, without localhost or pairing codes), and the Agent
   panel showing B's numbers in MERFISHEYES.

6. **An agent on the site.** What Tim described (confirmed by the user,
   2026-11-14): the person's agent, through the connector, talks to an
   agent that lives on the site and knows it — asks it a research
   question, and the site's agent works through the index, the stores and
   the viewer and answers, asynchronously. The person installs nothing.
   Its shape (one MCP tool that hands a request to the site's agent, or
   agent-to-agent messages), the model it runs, and who pays for it are
   open; E1–E4 are its plumbing either way.

Order inside the pillar: 1–3 first (they carry the 0.1 tools as they
are), 4 with B's long queries, 5 after D3, 6 once 1–4 work.

### Pillar F — Compute and writes

What no file records — DE between arbitrary selections, neighbourhoods,
spatial statistics — on the server (Batch, with a confirm and a cost) or
in the local session; write tools for projects, collections and dataset
metadata. Step 7, then steps 8–9 as written.

### Order

A is done; D1a is built. Next, in parallel:
- **METADATA.md steps 1–4** (schema on paper, bionty validation, metadata
  in the store, Lamin verification) — the foundation for B.
- **E1–E3** — the site as the MCP server, carrying the 0.1 tools, so a
  person adds a URL instead of installing anything.
- **D1b** — `.h5ad` and raw instrument output for the local viewer.

Then **B** on METADATA.md's index (steps 6–7), **D2–D5**, **E4–E5**, **F**.

**Rule for every pillar:** a new tool is added beside the 0.1 tools,
never by changing what they return; each ships with its tests and a live
check, and the package's minor version goes up with it.

## 1. The principle: the browser does the work, the agent directs it

The agent sends **intents** ("colour by `subclass`", "markers of cluster 12",
"violin of *Gad1* by `class`") and gets back **summaries** (tens of rows, a
link, a status). The data stays where it already is — in the tab, read
through the same readers the viewer uses. The agent never holds a cell ×
gene matrix in its context.

Consequences:

- Every capability is a MERFISHEYES feature first and a tool second. A tool
  that needs an analysis the site cannot do is a missing site feature, not a
  reason to ship data to the model.
- Tool results are bounded: a hard row cap, paging, and a `truncated` flag.
- When the agent must compute something itself (§4, tier 3), the data goes to
  a **file on the person's disk**, not into the conversation, and the result
  comes back to the panel as a plot spec or a small table.
- Nothing in the core is specific to one agent client. Claude-specific
  packaging (a skill, a plugin) wraps the server; it never replaces it.

## 2. Shape

```
agent (any MCP client; Claude Code first)
  │  MCP, stdio
merfisheyes-mcp  (local process, `npx`)          ── HTTPS + OAuth token ──▶ merfisheyes.com API
  │  bridge: WebSocket on 127.0.0.1, paired          (catalog, projects, writes, sidecar JSON)
browser tab on merfisheyes.com
  └─ executes commands against its stores; renders an Agent panel (split pane)
```

- **Local MCP + localhost bridge first** (D1). The tab connects to the
  local process after the person clicks *Connect agent* and confirms a pairing
  code. No new server infrastructure; exported data lands on local disk; the
  Vercel limits of §6 do not apply.
- **The command protocol is transport-independent.** A remote MCP at
  `merfisheyes.com/api/mcp` with a server relay (for claude.ai and other
  hosted clients) is the same commands over a different pipe — step 9, not a
  redesign.
- The tab is signed in already, so it opens the person's private datasets
  with its own session. The OAuth token (§5) is for what the MCP process
  asks the server directly: search, projects, writes, status.

## 3. Data sources the tools read

| What | Source | Notes |
|---|---|---|
| Dataset list, filters | Explore DB (`CatalogDataset`, `Dataset`, `Project`) via `lib/explore/search.ts` and `lib/explore/visibility.ts` | One source of truth for visibility rules (`isPublished`, `reviewStatus`, owner). |
| One BIL dataset's detail | `catalog.json` + `stats.json` + `qc.json` beside the store, over public HTTPS | Three small GETs, no Zarr opened. Written by every BIL pipeline run since 2026-11-05; see [LAMIN.md §2](LAMIN.md); schemas in [specs/schemas/](../specs/schemas/). **Missing files are normal** (legacy `bil-psc-data2/`, non-BIL). |
| Genes | `Dataset.genes` / `CatalogDataset.genes` | Already filled server-side (`lib/datasets/read-genes.ts`). |
| DE stats (markers, gene × cluster) | Legacy: `de/{col}.bin.gz`. SpatialData: `varm/de_mean__<col>`, `varm/de_pct__<col>`, listed in `uns/merfisheyes_de` (= `stats.json` `de_stats`) | Both already read by the viewer (`ChunkedDataAdapter`, `SpatialDataTableAdapter`, `lib/utils/de-stats.ts`). Tools go through these, never a second parser. |
| Everything else | The tab | Expression, coordinates, selections, plots. |

**Lamin is out of scope.** `merfisheyes/eyes-catalog` is a dev-side copy of
the sidecar JSON; the MCP reads the site's sources only.

Caveats from the sidecar files: `qc.json` keys `n_<level>` /
`mean_boot_<level>` vary with `mmc_levels` — iterate, never hard-code;
`qc.json` `n_genes` counts blanks and controls, `stats.json` `genes` does not;
the marmoset Xenium rows say `platform: "other"`.

Measured on the 14 BIL folders (2026-11-06, `scripts/validate_sidecars.py`):
all four files present and valid everywhere. `n_cluster` appears in 11 of
14 `qc.json`, `n_class`/`n_subclass`/`n_supertype` in 6, `mean_boot_*`
sparsely, `mean_corr` in none. **`de_stats` does not follow the taxonomy**:
ace-low-\* ran MapMyCells to `cluster` but carry no cluster DE; and the
three ace-owl-\* (Xenium) carry DE stats for QC columns
(`segmentation_method`, `nucleus_count`, `control_*_counts`, …) — which
`prepare.py` treats as categorical. Tools must not offer those as groupings
(MCP-PROTOCOL §3.4); the stores stay as they are (D13).

### What the site does today (audit, 2026-11-06)

Every route under `app/api` was read (58 files, 91 handlers); the
classification of each site action as a tool or excluded is
MCP-PROTOCOL §3.2. Facts that shaped the design:

- **A dataset has no visibility flag**: any `Dataset` opens by id
  (`/api/datasets/[id]`, `/object`, `/list`, `/single-molecule/[id]`) —
  the id is the only protection. The agent API does not inherit that
  (MCP-PROTOCOL §2, `lib/agent/visibility.ts`).
- **The catalog rule was implemented five times and disagreed**: the
  facets route and the Explore detail pages (`/explore/[id]`,
  `/explore/bil/[code]`) showed internal rows to anyone. Now one rule,
  `lib/explore/visibility.ts`, used by all of them and the agent API.
- **Admin power differs**: on datasets it covers admin-owned ones only;
  on projects any admin opens any project. An agent gets neither — it sees
  the person's own projects.
- **No route took a bearer token**; every helper calls `auth()` with no
  request. The agent API resolves its caller in one place
  (`lib/agent/caller.ts`). Write tools (step 7) must move each route's
  logic into `lib/` so the route and the tool share it.

## 4. When the site works and when the agent computes

This is the instruction the MCP server ships (server `instructions`, tool
descriptions, and a Claude Code skill), so an agent does not burn tokens
re-deriving what MERFISHEYES already has.

| Tier | Use when | Examples | Returns |
|---|---|---|---|
| 1 · view | Always, for anything the person should see | open dataset, colour by gene/column, filter cell types, split pane, camera, 2D/3D | ok + viewer URL |
| 2 · built-in analysis | A precomputed or in-tab statistic answers it | markers per cluster (DE), gene × cluster mean/pct, cluster composition, histogram, co-expression, cells in a selection | ≤ N rows, + a panel plot |
| 3 · agent compute | No tier-2 tool answers it | a custom test, a correlation the site does not offer, combining datasets | `export_table` → file path on disk; agent runs code; `show_plot` / `show_table` back to the panel |

Rules the instructions state: prefer the lowest tier; never ask for raw
per-cell data when a tier-2 tool exists; exports are capped and go to disk;
plots go to the panel, not the chat. Every tier-3 request that recurs is a
candidate tier-2 feature.

## 5. Authentication: OAuth 2.1

MERFISHEYES becomes an OAuth authorization server. NextAuth stays the login
("who are you", Google); a new layer issues **grants** ("what may this
agent do"): a consent page, scoped short-lived access tokens, refresh
tokens, a *Connected agents* list in `/account` with revoke.

- **Grant reach is the person's choice** (D3): at consent, *all my projects*
  or *these projects* (a picker). Both are built; the grant stores which.
  Public datasets are always readable.
- Local MCP: the loopback flow (as `gh auth login`): the process opens the
  browser, the person signs in and consents, the code returns to
  `127.0.0.1` with PKCE, the token is stored locally (OS keychain, else a
  0600 file).
- Remote MCP (step 9): the same server plus the MCP-spec discovery
  documents (`/.well-known/oauth-protected-resource`, authorization-server
  metadata) and client registration — what claude.ai and other hosted
  clients need.
- Scopes: `catalog:read` (public, no token needed), `projects:read`,
  `projects:write`, `datasets:write`, `viewer:control`, `export`;
  later `upload` (P1).
- **Every write is the site's own route logic** called with the grant's
  identity — ownership checks (`lib/datasets/owner.ts`) and project
  restriction applied in one place. The MCP never gets a privileged path.

## 6. Constraints

**Hosting (Vercel).** A request is a short-lived function: a duration cap, no
memory shared between requests or instances, no long-lived WebSocket server,
a ~4.5 MB response body. So: the MCP-facing API is stateless and returns
summaries and links; long work runs in the tab or as a Batch job, never in a
function; a server relay to a tab (remote MCP) needs an external pub/sub
(Redis, or a realtime service) — the reason the local bridge comes first.

**Contracts that must stay stable.** Three, each versioned:

1. *Tool names, inputs, outputs.* Models read them and people bake them into
   skills and prompts. Additive changes within a major version; deprecate
   before removing.
2. *The bridge protocol* between the `npx` package and the deployed site.
   They update independently (an old package against today's site), so the
   tab and the process exchange a protocol version and a capability list on
   connect, and refuse cleanly on mismatch.
3. *The data the tools read* — DE layouts, the sidecar JSON. These become
   an API the moment a tool reads them: JSON Schemas for `catalog.json`,
   `stats.json`, `qc.json` in [specs/schemas/](../specs/schemas/), and tools
   read through the viewer's readers so a format change is made once.

## 7. Decisions

| # | Decision | Status |
|---|---|---|
| D1 | Local MCP + localhost bridge first; remote MCP + relay later (step 9) | **taken** |
| D2 | Full OAuth 2.1 (loopback flow now, MCP-spec discovery for remote later); no personal tokens | **taken** |
| D3 | The person chooses at consent: all projects or selected ones; both built | **taken** |
| D4 | Panel plots: built-in plot components by name, Vega-Lite for custom ones, data by reference to the tab where possible | **taken** |
| D5 | The person uploads; the agent works on uploaded data. Agent upload is pinned (P1) | **taken** |
| D6 | Any MCP client eventually; Claude Code first for testing; nothing client-specific in the core | **taken** |
| D7 | Write tools in v1: anything the person can do on the site, within computational reason; destructive actions confirmed | **taken** |
| D8 | Catalog from the Explore DB + sidecar JSON; no Lamin | **taken** |
| D9 | OAuth server hand-built on Prisma (`lib/oauth/`), not a library: the MCP profile is narrow (public clients, code + PKCE S256, refresh rotation, CIMD + DCR), no maintained Node library fits Next route handlers without an Express shim, and every rule is pinned by `tests/unit/agent-oauth-*.test.ts` | **taken** (step 3) |
| D10 | The package lives in this repo, `packages/merfisheyes-mcp`, as its own folder with its own `package.json` and CI job (not an npm workspace: the site's install and Vercel build stay as they are). The package, the site's bridge executor and the protocol change together, in one PR; the package is published to npm from that folder under its own version | **taken** 2026-11-09 |
| D11 | Tool names carry the `meyes_` prefix (`meyes_search_datasets`, …): short, already the repo's shorthand, still distinct beside other servers | **taken** 2026-11-08 |
| D12 | The tool list of MCP-PROTOCOL §3 is approved for v1; more tier-2 analysis and tier-3 compute tools come later (pinned P2) | **taken** 2026-11-08 |
| D13 | Gaps in the BIL stores' DE stats (QC columns on ace-owl-\*, no cluster DE on ace-low-\*) are accepted as they are; the tools account for them, the pipeline is not changed for this | **taken** 2026-11-08 |
| D14 | Rate limits, in shared storage (serverless instances share nothing), checked in `lib/agent/route.ts` and the OAuth routes, answering 429 with the wait: agent API 60/min per IP without a token, 300/min per grant and 10,000/day per user with one; S3- or Batch-touching calls (detail, status) 30/min per caller; `/api/oauth/register` 10/hour per IP; `/api/oauth/token` 30/min per client and IP; consent decisions 20/hour per user. Built at the start of step 4, before the pipeline branch merges | **taken** 2026-11-08 |
| D15 | Personal datasets stay "anyone with the link" (any id opens on the site); the agent API keeps its stricter rule (MCP-PROTOCOL §2) | **taken** 2026-11-08 |
| D16 | The package is `@merfisheyes/mcp` on npm, under the `merfisheyes` npm organisation (created 2026-11-10, two-factor on) so it outlives any one maintainer; run with `npx -y @merfisheyes/mcp` (the command it installs is `merfisheyes-mcp`). First published by hand once the live site has the agent API (after `develop` → `main` and `OAUTH_ISSUER` in Production); later releases from GitHub Actions (npm trusted publishing). Until then it runs from the repository (`node packages/merfisheyes-mcp/dist/index.js`) | **taken** 2026-11-10 (publishing open) |
| D17 | Tokens are kept in a file only the person can read (`~/.config/merfisheyes-mcp/credentials.json`, 0600, written atomically), keyed by issuer — not the OS keychain, which needs native modules `npx` cannot rely on | **taken** 2026-11-09 |
| D18 | The work after step 4 is organised by five pillars (§0): A publish 0.1, B analyses from what a store records, C data, D bridge and Agent panel, E compute and writes. B — the analysis half of step 6, served by the agent API from files the pipeline writes — comes **before** the bridge (step 5): it works in every client with no tab open, and the third manual test showed it is what agents lacked. New tools are added beside the 0.1 tools, never by changing what those return | **taken** 2026-11-12 |
| D19 | The MCP serves two goals: research across the whole library (queries over standardised metadata and per-dataset statistics, on the server, returning collections) and MERFISHEYES on a researcher's own data (the viewer on `localhost`, driven by their agent inside an analysis). The pillars of §0 are rewritten around them; D18's order is replaced by §0's | **taken** 2026-11-13 |
| D20 | MERFISHEYES's metadata model is SpatialData-DB's catalog schema (ontology ids and hierarchies, `timepoint_days`, genotype, treatment, condition, technology vocabulary) and its gene index, for every dataset, not a schema of our own | **taken** 2026-11-14, after the meeting with Tim Treis; carried out by [METADATA.md](METADATA.md) (M1–M7) |
| D21 | Local use is an npm package and a browser, not a Python package; it serves the site's viewer code (one viewer, never a fork): first the hosted viewer reading `127.0.0.1`, then a standalone build of the same components shipped with the package | **taken** 2026-11-13 (the standalone build's shape open) |
| D22 | The site serves the MCP at `www.merfisheyes.com/mcp` (remote, Streamable HTTP, our OAuth server) as the main way agents reach MERFISHEYES; the npm package is the local companion (files on disk, the local viewer) and a fallback for clients without remote MCP | **taken** 2026-11-14 |
| D23 | One implementation of each tool's definition (name, description, schema, formatting) shared by the package and `/mcp` | **taken** 2026-11-14 (the module's shape open, pillar E2) |

## 8. Steps

Each step ends at a gate; the next does not start until the gate passes.

### Step 1 — Contracts on paper

Write `docs/specs/MCP-PROTOCOL.md`: every tool (name, tier, scope, input,
output, row cap), the bridge messages (handshake, command, result, event,
error), versioning rules. JSON Schemas in `docs/specs/schemas/` for the
bridge messages and for `catalog.json`, `stats.json`, `qc.json`. Audit
`app/api/**` and list each site action as a tool or as excluded, with the
reason (D7).

**Gate:** the schemas validate the sidecar files of all 14 BIL datasets in
`bil-psc-zarr/`; the tool list is reviewed by the user.

**Done (2026-11-06).** [MCP-PROTOCOL.md](../specs/MCP-PROTOCOL.md);
`docs/specs/schemas/{stats,qc,catalog}.schema.json`;
`scripts/validate_sidecars.py` (reuses `validate_store.schema_errors`);
`tests/python/test_sidecars.py` with two committed fixtures
(`tests/fixtures/sidecars/ace-dip-vat`, `ace-owl-cup`). Gate: all 14
folders pass over HTTPS, 70 checks, no FAIL or WARN; 15 pytest cases.
The tool list was reviewed and approved on 2026-11-08 (D12).

### Step 2 — One visibility query, one agent API

Extract the "datasets this caller may see" rule from Explore into one
function used by Explore and every agent route. Add the read routes the
MCP calls: search, facets, dataset detail (DB row + sidecar JSON merged,
missing files tolerated), genes, my projects. Stateless, bounded, cached.

**Gate:** unit tests prove unpublished, unapproved and other people's
private rows never come back; Explore's results are unchanged.

**Done (2026-11-06).** `lib/explore/visibility.ts` (the catalog rule) and
`lib/explore/search.ts` (Explore's search, moved out of the route) — now
used by `/api/explore`, `/api/explore/[id]`, `/api/catalog/facets`, the
SSR Explore page and both detail pages. `lib/agent/` (caller, visibility,
sidecars, shapes, responses) and eleven routes under `/api/agent/v1`
(MCP-PROTOCOL §4). `ownedDatasetsWhere` moved to `lib/datasets/owner.ts`
and shared with `/api/ingest/mine`. Sidecar files are fetched through
`lib/net/limited-fetch.ts` (HTTPS, S3 hosts only, 256 KB, 5 s, no private
addresses). Gate: `tests/unit/agent-oauth-db.test.ts` on a real
PostgreSQL — the catalog, dataset and project visibility cases, Explore's
list, detail and facets against the same rule; a deliberately broken
reach rule and a deliberately broken internal-row rule each fail it. It
runs in CI: the unit-tests job starts a Postgres service, migrates it and
sets `TEST_DATABASE_URL` (`.github/workflows/test.yml`).
Behaviour that changed on purpose: internal rows no longer appear in
`/api/catalog/facets` or on `/explore/[id]`, `/explore/bil/[code]` for
non-admins.

### Step 3 — OAuth authorization server

Prisma: `OAuthClient`, `OAuthGrant` (user, client, scopes, all-projects or
a project list), hashed access and refresh tokens. Routes: authorize
(NextAuth sign-in → consent page with the project choice), token (PKCE),
revoke, metadata. `/account` → *Connected agents* with revoke. One bearer
helper that every agent route uses. Settle D9 at the start.

**Gate:** tests for PKCE, expiry, refresh, revoke, scope denial, project
restriction; an e2e run of the consent flow.

**Done (2026-11-06), following the MCP authorization draft read that
day** (CIMD preferred, DCR kept, PKCE S256 only, RFC 8707 resource
binding, RFC 9207 `iss`, rotating refresh tokens). Migration
`20261106211012_agent_oauth` (four tables, additive). `lib/oauth/`
(config, crypto, clients, flow, http); `/oauth/authorize` (sign-in, then
consent: client name and host, redirect host, a warning for loopback-only
clients, scopes in words, all projects or chosen ones);
`/api/oauth/{authorize,token,register,revoke}`; discovery at
`/.well-known/oauth-authorization-server` and
`/.well-known/oauth-protected-resource[/api/agent/v1]` (rewrites in
`next.config.js`); *Connected agents* on `/account`
(`/api/account/agents`). Each deployment sets `OAUTH_ISSUER` (else
`NEXT_PUBLIC_APP_URL`) to its own origin. Gate: 22 database tests
(PKCE, single-use codes with replay revocation, resource binding,
rotation with reuse detection, narrowing-only refresh, expiry, revoke,
disconnect, scope challenges 401/403, the consent POST's origin,
signature, expiry and person checks, discovery fields) and 13 pure-rule
tests; type-check, lint and `next build` clean; a smoke run of the built
app over HTTP (discovery through the rewrites, the 401 challenge,
registration, the signed-out redirect to sign-in, an untrusted redirect
refused). And in a browser, `e2e/tests/upload/agent-consent.spec.ts` (2026-11-08):
sign in, register an agent, the consent page names it, its return host and
the loopback warning; allow one project; swap the code; the token reaches
that project and not another; Disconnect on `/account` refuses it at once.
It runs in CI's e2e-upload job (tag `@upload`) and locally on any migrated
database with `E2E_AGENT=1` — no Docker needed. The e2e server now sets
`OAUTH_ISSUER` to its own origin (`playwright.config.ts`).

### Step 4 — The MCP package (read tools)

Stdio server on the MCP TypeScript SDK v2 (`@modelcontextprotocol/server`);
loopback OAuth; token storage; the catalog and project tools of step 2; the
§4 instructions. In `packages/merfisheyes-mcp` (D10).

**Done first, before the pipeline branch merged (2026-11-09):** the rate
limits of D14 (`lib/rate-limit.ts`, migration
`20261109081522_agent_rate_limits`, MCP-PROTOCOL §6.5), tested in
`tests/unit/agent-oauth-db.test.ts`.

**Gate:** MCP Inspector lists the tools; in Claude Code (`claude mcp add`)
an agent finds a public dataset and one of the user's own.

**Done (2026-11-09).** `packages/merfisheyes-mcp` (npm name
`@merfisheyes/mcp`, not yet published — D16): stdio server on
`@modelcontextprotocol/server` 2.4; eleven tools (the nine of
MCP-PROTOCOL §3.1 plus `meyes_sign_in`, `meyes_sign_out`); server
instructions carrying §4's rules for what exists so far; sign-in as an
OAuth public client per the MCP spec (discovery with issuer checks, the
site's Client ID Metadata Document `app/oauth/merfisheyes-mcp.json` over
HTTPS, registration over http, PKCE, `resource`, `iss` checked, refresh
rotation) answered on a loopback listener — the tool returns the link at
once and the listener finishes in the background; tokens in a 0600 file
(D17); `merfisheyes-mcp login | logout | whoami` for a terminal. Root
`tsconfig` excludes `packages/`; CI: an `MCP package` job (build + tests),
and the upload E2E job builds it. Gate:
- MCP Inspector (CLI) lists the eleven tools with their annotations and
  returns an out-of-range argument as a readable tool error.
- `e2e/tests/upload/agent-mcp.spec.ts`: the built package over stdio
  against the site — public search with no sign-in; told to sign in;
  signed in through a browser allowing one project; reaches that project
  and its dataset and no other; signs out.
- In Claude Code (`claude -p --mcp-config … --strict-mcp-config`, a local
  site): asked to find a public mouse hippocampus dataset and the person's
  own datasets, it called `meyes_search_datasets`, `meyes_whoami` and
  `meyes_list_my_datasets` and answered both from the site.
- Package unit tests (sign-in against a stand-in server: forged state and
  foreign `iss` refused, PKCE, rotation, revoke, file mode; paging under
  the character limit); a database test that the authorization server
  accepts the site's own client document.

Released to the live site on 2026-11-10 (PR #202): the package works
against www.merfisheyes.com from a checkout; npm publishing (D16) is
still to do.

**First manual test, on staging (2026-11-10 → 11).** A person on a Mac,
the package from branch `feature/mcp` (commit `93dcf87` added
`MERFISHEYES_VERCEL_BYPASS`, since staging sits behind Vercel's Deployment
Protection), driven from Claude Code in VS Code.

Worked:
- The package built on the Mac (Node ≥ 20.3) and Claude Code listed the
  `meyes_*` tools.
- Catalog search through Claude worked on staging's real catalog
  ("mouse brain": 1,332 with `species`/`tissue` filters, 1,349 free text),
  with paging.

Failed — fix before the next test:
- **F1 · Terminal sign-in stopped at the consent page** with "the
  redirect_uri is not one this client registered". The code on both sides
  agrees (registration lists `http://127.0.0.1/callback`, the site accepts
  any loopback port, the E2E passes), so the cause is not found yet. The
  leading suspect is the link losing `redirect_uri` on its way through the
  Vercel and Google logins. `validateAuthorize` gives the same message for
  a missing `redirect_uri` as for a wrong one, so the page should tell
  them apart.
- **F2 · Claude Code's package was not signed in, and registered a new
  client.** The package registers only when its credentials file holds
  nothing for the issuer, so that process found no saved sign-in. This is
  consistent with F1, but the terminal and Claude Code reading different
  files (or different issuers) is not ruled out. `whoami` should print the
  credentials path and the issuer it looked under.
- **F3 · The scopes surprise.** `meyes_sign_in` with
  `scopes=["projects:read"]` asked for `catalog:read projects:read`. Adding
  `catalog:read` is deliberate (DEFAULT_SCOPES), but the tool's answer
  should say so.
- **F4 · Setup friction.** Commands with a `# comment` broke in zsh
  (interactive comments are off; the apostrophe opened a quote). The
  README's commands must have no inline comments, and should give
  `node -v` (≥ 20.3) as the first step.

Problems in the site's data that the agent ran into are logged in §10.

Decided after it: manual testing moves off staging, whose Vercel login,
bypass secret and separate database add failure points that are not the
package's. Since `develop` was released on 2026-11-10 (PR #202, which
carries #205), the live site has the agent API, and testing moved there.

**Second manual test, on the live site (2026-11-11).** Same Mac, package
from `feature/mcp`, no settings (the default, www.merfisheyes.com; the
site's Client ID Metadata Document, no registration). Passed: sign-in,
then in Claude Code the person's two projects (and their datasets),
catalog search ("mouse brain": 1,335) and one dataset's detail — F1 and
F2 did not recur, so both were staging's.

Found:
- **F5 · A search page of 5, not 20.** `many()` halves a page until it is
  under `CHARACTER_LIMIT`; an item that lists every molecule section
  (58 for ace-den-fog) makes the page too long. The search item should
  say how many sections and leave their links to `meyes_get_dataset`.
- Data problems: logged in §10.
- The account page's agent list was reworded to be short: name,
  permissions, which projects, last used.

**Polish for the demo (2026-11-11).** F3–F5 fixed in the package: search
cards keep one viewer link per kind with the counts (`compactCard`), so a
page holds 20; `meyes_sign_in` says what it asks for and that
`catalog:read projects:read` are always included; `whoami` names the
credentials file it read; the README has no inline comments and starts
with `node -v`; on Windows the browser opens through `rundll32` (`cmd
/c start` cut links at the first `&`). New tool **`meyes_open_view`**
(MCP-PROTOCOL §3.1): opens a dataset in the person's browser with genes,
column, cell types and 2D/3D set through the viewer's `v=` link state —
the part of step 5 that needs no bridge. Gate: unit tests
(`test/view.test.ts`); against the live site, ace-den-fog's cells link
came back from the viewer with `Slc17a7` applied and its expression range
filled in, and section 609889HK01 drew Gad1 and Slc17a7 in 3D.

**Third manual test, the live site (2026-11-12), and what followed.** The
tools answered search, details, "show me gene X" and "are you signed in"
in seconds. Three of five questions had no tool, and Claude Code left the
tools for the shell: it read molecule manifests and DE files from S3 and
decoded them by hand ("most abundant genes", twice) and read slice lists
("the combined atlas") — 1–1½ minutes each against seconds, against §1's
principle, and impossible in a client without a shell. What it lacked —
gene abundance, dataset structure (one slice or combined, slice count,
obs columns and their values), counts per facet, and the plot panel's
numbers (cell-type counts, a gene's expression by cell type) — is step
6's analysis tools, served by the agent API from the files the site
already writes. Data problems it met are in §10 (E10–E14).

Changed after it (package only): `meyes_open_view` returns a link and
opens a tab only with `open=true`, which the instructions reserve for an
explicit ask, one tab per request (three genes had opened three tabs);
it takes `compare` for split screen (`split`/`splitS3Url`, `splitType`,
`sync`, `rv=`), `dot_size` and `opacity` / `unselected_opacity`;
`meyes_sign_in` waits up to two minutes for the person and answers who
is signed in. Opening the plot panel from a link was dropped: the
person wants the agent to use what the panel computes, which is step 6.
Gate: unit tests; on the live site a split screen of two ace-den-fog
sections came up with two panels, both genes on each, cameras synced.

### Step 5 — The bridge and tier-1 view control

Tab side: *Connect agent*, pairing code, WebSocket to `127.0.0.1`, the
version + capability handshake, an "agent connected" indicator and
disconnect. A command executor mapped onto store actions through the panel
context (not `pathname`). Tier-1 commands: open dataset, colour by
gene/column, filter, split pane, camera, 2D/3D.

**Gate:** e2e on the tiny fixtures with a scripted bridge; works in
Chrome through its local-network-access prompt; Firefox and Safari tried
and the result recorded here.

### Step 6 — The Agent panel and tier-2 analysis

A panel in the split pane rendering built-in plots by name and Vega-Lite
specs. Tier-2 commands for both formats through the existing readers:
markers, gene × cluster, composition, histogram, co-expression, selection
summary. Row caps everywhere.

**Gate:** each command returns the same numbers the viewer shows, on a
legacy and a SpatialData fixture; no result exceeds its cap.

### Step 7 — Tier 3 and write tools

`export_table` (capped, to disk), `show_plot`, `show_table`. The write tools
from step 1's audit, each calling the site's route logic; destructive
actions (delete, unpublish) require an explicit confirm argument.

**Gate:** every write tool honours scope and project restriction in tests;
a destructive call without confirm is refused.

### Step 8 — Instructions, packaging, a token budget

Final tool descriptions and server instructions; a Claude Code skill or
plugin wrapping the server; the user guide page `/docs/agents`. A
reference analysis run end to end, its token count recorded in
[BUDGETS.md](../specs/BUDGETS.md).

**Gate:** the reference analysis completes within the recorded budget,
with the result shown in the panel.

### Step 9 — Other clients and remote MCP

*Now pillar E (2026-11-14): the remote server comes first, not last.*

Claude Desktop and another MCP client against the stdio server. Then
`/api/mcp` (stateless Streamable HTTP), the MCP-spec OAuth discovery and
client registration, and a relay to the tab over an external pub/sub —
the same bridge protocol.

**Gate:** the reference analysis of step 8 runs from a second client and
from claude.ai.

### Pinned

- **P1 · Agent upload.** The local process can read the person's files, so
  an `upload` tool can pass only a path and use the existing presigned and
  multipart routes (`lib/upload/`). Needs: resume, progress events, the
  `upload` scope, and a confirm before a Batch job starts.
- **P2 · More analysis and compute tools.** The tier-2 and tier-3 lists of
  MCP-PROTOCOL §3.4–3.5 are the v1 set (D12); grow them from what agents
  actually ask for (every recurring tier-3 request is a tier-2 candidate,
  §4).

### After step 4

The work from here is organised by the pillars at the top of this file
(§0): steps 5–9 below are their detail.

## 9. Risks

- **Visibility leak**: one shared query (step 2), tested.
- **Token drain**: row caps, paging, §4 tiers in the instructions, a
  measured budget (step 8).
- **Localhost from a public page**: browsers gate a public site reaching
  `127.0.0.1` (Chrome's local-network-access prompt; Safari and Firefox to be
  tested in step 5). The pairing step must survive that prompt.
- **A paired tab is a remote control**: pairing code, origin check,
  `viewer:control` scope, a visible indicator and a disconnect button.
- **Writes on the person's behalf**: same route logic as the site, scoped
  grants, confirm on destructive actions, revoke in `/account`.
- **Cost**: server reads of S3 are cached and rate-limited; Batch jobs only
  start with a confirm. (Cached: yes, per instance. Rate-limited: not yet —
  the only limiter on the site is the sign-in one; add one before launch.)
- **Holes the audit found, older than this work, to close before agents
  are announced** (not fixed on this branch): `POST /api/send-email` and
  `/api/send-email-single-molecule` send mail to any address with no
  sign-in; `POST /api/datasets/[id]/complete`, `files/complete`,
  `files/[key]/complete`, `upload-progress`,
  `/api/single-molecule/[id]/complete` and its `files/[key]/complete` need
  no sign-in or ownership; `/api/datasets/check-duplicate/[fingerprint]`
  reveals any dataset's id and title.
- **CIMD fetch and DNS rebinding**: `lib/net/limited-fetch.ts` checks the
  resolved addresses, then `fetch` resolves again — a narrowed, not closed,
  window. Acceptable for public metadata; revisit if anything secret is
  ever fetched this way.

## 10. Site errors the agent finds (running log)

Until the agent API, these errors could only be found by reading the
database or S3 directly. An agent sees the site as a visitor does, one
dataset at a time, so it notices typos, wrong names, duplicates and
inconsistencies in passing. **Add each one here when it is met**: date,
where, what, and what fixing it touches. Mark it fixed with the commit
instead of deleting it. A deliberate sweep (an agent walking the whole
catalog) is a later step; this log is what such a sweep would add to.

| # | Found | Where | What | Fix touches | Status |
|---|---|---|---|---|---|
| E1 | 2026-11-10 | Catalog cards | No cell and/or gene count on some cards: two SpatialData-DB stores have neither; the BIL MERSCOPE atlas entries (ace-low-bar, ace-low-zip, ace-den-fog) have no cell count, though ace-low-bar has a cells viewer | the stores' `stats.json`, or how the card reads it | open |
| E2 | 2026-11-10 | Search ranking | For "mouse brain" the BIL MERFISH atlases come 4th and 5th, behind CUT&Tag and ATAC entries with no counts | `lib/explore/search.ts` ordering | open |
| E3 | 2026-11-10 | Platform facet | `xenium` and `Xenium` both listed, splitting the filter | catalog rows; canonical spellings in `docs/specs/conventions.md` | open |
| E4 | 2026-11-10 | Species facet | `Drosophila` / `Drosophila melanogaster`; `Macaque` / `Rhesus Macaque` | as E3 | open |
| E5 | 2026-11-10 | Tissue facet | `otak` (Malay/Indonesian for "brain"); `Whole embryo, heart` (two values in one) | as E3 | open |
| E6 | 2026-11-11 | Species facet | `Mouse` and `Mice` both listed; `Mice` + `Brain` finds nothing | as E3 | open |
| E7 | 2026-11-11 | Catalog entry "Re-activation of neurogenic niches" | A "Xenium Human Breast Cancer (molecules)" entry is attached that does not belong to it | that catalog row's entries | open |
| E8 | 2026-11-11 | ace-den-fog | 702 genes in the store; the description says 500 were profiled | check which is right; the description or the store | to check |
| E9 | 2026-11-11 | Catalog entry "Mouse Brain Saggital 10 Slices" | Typo in the title: "Saggital" for "Sagittal" | that catalog row's title | open |
| E10 | 2026-11-12 | Developing mouse atlas (P28) | No combined dataset: its sections are ~57 single-slice entries (ace-low-\*, ace-lot-\*), e.g. ace-low-bag holds one slice, 140,542 cells; the combined developing datasets that exist (ace-let-law/lax/lay, ace-met-war) are P0/P4 | combine the P28 slices (`combine_slices_v3.py`) and register it | open |
| E11 | 2026-11-12 | ace-let-law, ace-let-lax, ace-let-lay, ace-met-war | Combined datasets with no cell count recorded and no molecule views | their `stats.json`; molecules per slice | open |
| E12 | 2026-11-12 | ace-met-war | Titled "741564", unlike its siblings' naming | its catalog title | open |
| E13 | 2026-11-12 | BIL legacy datasets | No sidecar files (`stats.json`, `catalog.json`, `qc.json`), so no counts, slice numbers or cell-type levels reach the agent or the cards | reprocessing with the BIL pipeline writes them | open |
| E14 | 2026-11-12 | Tissue facet | Brain regions are separate tissues (Hippocampal formation, Cerebral cortex, Striatum, …), so tissue="Brain" misses them | a tissue hierarchy, or "Brain" on region datasets too | open |

