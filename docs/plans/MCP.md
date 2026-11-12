# Plan: an agent that drives MERFISHEYES (MCP)

Let a person ask an agent — Claude Code first, then any MCP client — to run
an analysis on a dataset, and have the agent do it **through a MERFISHEYES
browser tab**: open datasets, colour, select, split panes, and show the
results — plots and statistics — in a panel beside the scene. Public datasets
and the person's own projects are both in reach, and the agent may do
anything the person could do on the site, within computational reason.

This file is the brief for every session working on it: the facts, the
decisions taken, the steps with their gates. **Status (2026-11-09): steps
1–4 done.** Keep it
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
| D16 | The package is `@merfisheyes/mcp` on npm, under the `merfisheyes` npm organisation (created 2026-11-10, two-factor on) so it outlives any one maintainer; run with `npx -y @merfisheyes/mcp` (the command it installs is `merfisheyes-mcp`). First published by hand once the live site has the agent API (after `develop` → `main` and `OAUTH_ISSUER` in Production); later releases from GitHub Actions (npm trusted publishing). Until then it runs from the repository (`node packages/merfisheyes-mcp/dist/index.js`) | **taken** 2026-11-10 (publishing open) |
| D17 | Tokens are kept in a file only the person can read (`~/.config/merfisheyes-mcp/credentials.json`, 0600, written atomically), keyed by issuer — not the OS keychain, which needs native modules `npx` cannot rely on | **taken** 2026-11-09 |

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

### Next, from 2026-11-12: publish, then real analyses

Steps 1–4 and the link-built views let an agent explore the site the way
a person can. The goal is analyses a person could not do by clicking.
The third manual test showed what stands between the two: the numbers an
analysis needs exist in files the pipeline writes, but no tool serves
them, so the agent fetched and decoded them itself. The order below
serves them first, through the agent API, before any bridge — it works
in every client, with no tab open, in seconds.

**A · Publish 0.1 (by the 2026-11-14 demo).**
1. PR `feature/mcp` → `develop`; CI green; merge.
2. Release `develop` → `main`: brings the molecule-opacity link field
   (until then the live molecule viewer ignores `op`) and the account
   page's agent list.
3. `npm publish` from `packages/merfisheyes-mcp` by hand (`npm login`,
   two-factor); then on a clean machine `claude mcp add --scope user
   merfisheyes -- npx -y @merfisheyes/mcp` lists the 12 tools. The
   README's "until published" section goes.
4. A demo script: prompts the tools answer directly (search, detail,
   views, split screen, sign-in, projects). Not yet: "most abundant",
   "which dataset combines the slices" — those are B.
5. Later: publishing from GitHub Actions (npm trusted publishing) on an
   `mcp-v*` tag; semver, a minor version per batch of new tools.

**B · Analyses from what a store records (step 6, first half; no bridge).**
New agent-API routes (`/api/agent/v1/datasets/{id}/…`, `expensive` rate
tier, cached per store) reading the files the pipeline already writes,
and a tool for each. Every answer says where its numbers come from, and
"not recorded for this dataset" when the file is missing — never a guess.

| Tool | Answers | Read from |
|---|---|---|
| `meyes_dataset_structure` | one slice or combined; slices and their ids; obs columns with type and number of values; embeddings; images; molecule sections | `obs/metadata.json`, `_sample_id`, manifest, `stats.json`, store metadata |
| `meyes_gene_stats` | per gene: total counts, cells expressing, molecules (in and out of cells); "most abundant", "is X expressed" | DE stats (counts × means), molecule `manifest.json` `molecule_counts` |
| `meyes_composition` | cells per category of a column, overall or per slice — the plot panel's bar plot | DE stats category sizes, obs columns |
| `meyes_gene_by_group` | a gene's mean and % expressing per category — the plot panel's gene × cell type | DE stats |
| `meyes_markers` | top genes per category against the rest | DE stats (fold change and % computed on the server) |
| `meyes_compare` | two datasets, slices or conditions: composition and gene differences | the above, per side |

Charts: each analysis tool can return a Vega-Lite spec with its rows
(D4), which Claude Desktop and claude.ai draw. Gate: each number equals
the viewer's (plot panel, legend) on a legacy and a SpatialData fixture;
an answer stays under the character limit; a store missing a file gets
"not recorded", tested.

**C · Data the analyses can stand on (in parallel; BIL-CLUSTER.md).**
Most of B's "not recorded" answers are E1–E14: legacy BIL datasets
without sidecars, missing counts, no combined P28 atlas, facet
spellings. Reprocessing with `prepare.py` / the BIL pipeline writes
`stats.json`, `catalog.json`, `qc.json` and DE for every store; B is
built against stores that have them and falls back on the others. An
agent sweep — a scripted run of the read tools over the whole catalog —
fills §10 instead of waiting for an agent to stumble on each error.

**D · The bridge (step 5) and the Agent panel (step 6, second half).**
What links cannot do: change a tab already open, read what the person is
looking at ("what's selected?", "this region"), analyses of an arbitrary
selection the browser holds, and plots drawn inside MERFISHEYES. After B,
so the panel shows B's numbers rather than inventing its own.

**E · Compute and writes (step 7).** What no file records: DE between
arbitrary selections, neighbourhoods, spatial statistics — on the server
(Batch, with a confirm and a cost) or in the tab. Write tools for
projects and dataset metadata. Then steps 8–9 as written.

Rule for all of it: a new tool is added beside the 0.1 tools, never by
changing what they return; each ships with its tests and a live check,
and the package's minor version goes up with it.

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

