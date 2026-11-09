# Plan: an agent that drives MERFISHEYES (MCP)

Let a person ask an agent — Claude Code first, then any MCP client — to run
an analysis on a dataset, and have the agent do it **through a MERFISHEYES
browser tab**: open datasets, colour, select, split panes, and show the
results — plots and statistics — in a panel beside the scene. Public datasets
and the person's own projects are both in reach, and the agent may do
anything the person could do on the site, within computational reason.

This file is the brief for every session working on it: the facts, the
decisions taken, the steps with their gates. **Status (2026-11-08): steps
1–3 done.** Keep it
current: a step that changes a fact or a decision edits this file in the
same commit. Branch: `feature/mcp` (from `bil/pipeline-dry-run`). The
contract the steps build is [docs/specs/MCP-PROTOCOL.md](../specs/MCP-PROTOCOL.md).

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
