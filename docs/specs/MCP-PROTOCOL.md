# MCP protocol — what an agent may see and do, and how it asks

**Status: draft 0.1 (2026-11-06).** The contract between three parties:
the MCP server an agent runs (`merfisheyes-mcp`, local, stdio), the
MERFISHEYES API it calls (`/api/agent/v1`, OAuth), and the browser tab it
drives (the bridge). The plan that builds it is
[docs/plans/MCP.md](../plans/MCP.md); the principle behind every choice
here is its §1 — *the browser does the work, the agent directs it*.

Three things here are public contracts the moment a client depends on
them, and change only under §8: the **tool catalogue** (§3), the **agent
API** (§4) and the **bridge messages** (§7).

## 1. Who the caller is

Every agent request resolves to one **caller** (`lib/agent/caller.ts`):

| Field | Meaning |
|---|---|
| `userId`, `role` | The person, as the site knows them. `null` for an anonymous caller. |
| `via` | `oauth` (a bearer token from §6) or `session` (the site's own cookie — the browser tab). |
| `scopes` | What the grant allows (§6.3). A session caller holds every scope. |
| `projects` | `"all"`, or the project ids the person chose at consent. A session caller holds `"all"`. |

An anonymous caller holds `catalog:read` only.

## 2. What a caller may see

The site has two visibility rules, and the agent uses each without a copy:

1. **Catalog rows** (Explore): `lib/explore/visibility.ts`. A row is
   visible when it is published, not internal (internal: admins only, and
   only when they ask for it), and curated or an *approved* community
   submission; a list additionally needs one entry that can be opened.
   Explore's list, featured row, filters, detail route, detail pages and
   `/api/catalog/facets` all use it.
2. **Datasets** (`Dataset` rows): the site has **no visibility flag** —
   any dataset is readable by id, and the id is the only protection. The
   agent does not inherit that. A dataset is visible to a caller when:
   - an entry of a catalog row visible under rule 1 names it (by
     `datasetId`, or by `s3BaseUrl` for curated rows), **or**
   - the caller owns it (`ownerId`), or it is `adminOwned` and the caller is
     an admin — **and**, when the grant is limited to projects, it is a
     member of one of them.

   `lib/agent/visibility.ts` is that rule; nothing else decides it.
3. **Projects**: the caller's own (`ownerId`), intersected with the grant's
   projects. The site lets any admin open any project; an agent does not
   get that power.

Something not visible is answered exactly as something that does not
exist (404), so ids cannot be probed.

## 3. Tool catalogue

Names carry the `meyes_` prefix (short for MERFISHEYES; the server's display name stays "MERFISHEYES"). Every list takes `limit` (default
20, max 100) and `offset` and returns `total`, `count`, `offset`,
`has_more`, `next_offset`. Every tool takes `response_format`
(`markdown` default, `json`). A result stays under 25,000 characters by
returning fewer items, never by cutting text.

**Tier** is §4 of the plan: 1 view, 2 built-in analysis, 3 agent
compute. **Via**: `api` — the MCP process calls §4; `tab` — a command
over the bridge (§7). **Ann.**: `R` read-only, `W` writes, `D`
destructive (needs `confirm: true`), `I` idempotent.

### 3.1 Catalog and account (step 2, read)

| Tool | Scope | Via | Ann. | Does |
|---|---|---|---|---|
| `meyes_search_datasets` | `catalog:read` | api | R I | Explore search: text, species, tissue, platform, genes (substring) or genes_exact, dataset type, category (all, featured, bil, sddb, community). Cards, no gene lists. |
| `meyes_list_facets` | `catalog:read` | api | R I | The species, tissue and platform values in use — call before filtering so the spelling matches. |
| `meyes_get_dataset` | `catalog:read` | api | R I | One catalog row or app dataset: metadata, entries (each with its viewer URL), and for a store with sidecar files the merged `catalog.json` + `stats.json` + `qc.json` (cells, genes, molecules, sections, DE columns, QC, taxonomy levels). |
| `meyes_list_genes` | `catalog:read` | api | R I | A dataset's genes, filtered by a substring, paged. |
| `meyes_whoami` | — | api | R I | The person, the scopes and the projects the token carries. |
| `meyes_list_my_datasets` | `projects:read` | api | R I | The caller's datasets (in scope), with status. |
| `meyes_get_processing_status` | `projects:read` | api | R I | Status and progress of a dataset being processed. |
| `meyes_list_projects` | `projects:read` | api | R I | The caller's projects in scope, with member datasets. |
| `meyes_get_project` | `projects:read` | api | R I | One project with its datasets. |

### 3.2 Account writes (step 7)

Each is the site's own route logic, moved into `lib/` and called with the
caller (plan §5) — the route and the tool are one implementation.

| Tool | Scope | Ann. | Site action it is |
|---|---|---|---|
| `meyes_create_project` | `projects:write` | W | `POST /api/projects`. A project created under a project-limited grant joins the grant. |
| `meyes_update_project` | `projects:write` | W I | `PATCH /api/projects/[id]` |
| `meyes_add_dataset_to_project` | `projects:write` | W I | `POST /api/projects/[id]/datasets` |
| `meyes_remove_dataset_from_project` | `projects:write` | W I | `DELETE /api/projects/[id]/datasets/[datasetId]` (reversible) |
| `meyes_delete_project` | `projects:write` | D | `DELETE /api/projects/[id]` — datasets are kept |
| `meyes_update_dataset` | `datasets:write` | W I | `PATCH /api/ingest/mine/[id]`: title, description, species, tissue, platform, tags, links, metadata, viewer defaults |
| `meyes_set_image_overlay` | `datasets:write` | W I | `PATCH /api/datasets/[id]/overlay` |
| `meyes_attach_molecules` | `datasets:write` | W | `POST /api/datasets/[id]/overlay` (groups both in a project) |
| `meyes_detach_molecules` | `datasets:write` | D | `DELETE /api/datasets/[id]/overlay` |
| `meyes_register_url` | `datasets:write` | W | `POST /api/datasets/register-s3` — add a public store to the account |
| `meyes_submit_to_explore` | `datasets:write` | W | `POST /api/community/submit` (dataset or project; goes to review) |
| `meyes_withdraw_submission` | `datasets:write` | D | `DELETE /api/community/submit` |
| `meyes_delete_dataset` | `datasets:write` | D | `DELETE /api/ingest/[id]` — permanent; refused while on Explore |
| `meyes_set_thumbnail` | `datasets:write` | W I | `POST /api/datasets/[id]/thumbnail`, the bytes captured by the tab |

**Excluded**, with the reason:

| Site action | Why not a tool |
|---|---|
| everything under `/api/admin` | Admin surfaces; an agent acts as a person, never as staff. |
| `/api/auth/*` | The session issuer. |
| upload plumbing: `*/initiate*`, `*/complete`, `files/*`, `upload-policy`, `upload-progress`, `check-duplicate`, `ingest/[id]/abort` | Moves bytes. Pinned as plan P1 (agent upload). |
| `datasets/[id]/object`, `/list`, `single-molecule/[id]/gene/*` | Byte access; the tab reads data, the agent does not. |
| `datasets/[id]/view` | Analytics; an agent would inflate view counts. |
| `ingest/[id]/callback`, `send-email*`, `uploads/failed` | Worker and mail plumbing. |

### 3.3 Viewer (step 5, tier 1, via tab, scope `viewer:control`)

| Tool | Ann. | Does |
|---|---|---|
| `meyes_viewer_status` | R | Connected tab(s), what each panel shows, the current view state. |
| `meyes_open_dataset` | W I | Open a dataset (left or right panel; right opens split screen). |
| `meyes_set_view` | W I | One partial update of the view: colour by gene / column, second gene (co-expression), shown categories, colormap, scales, point size, 2D/3D, rotation, flips; molecule genes and colours. Mirrors `CellVizUrlState` / `SMVizUrlState`. |
| `meyes_set_camera` | W I | Fit, or frame a box in microns. |
| `meyes_select` | W | Select cells by a polygon or box in microns, or by category; returns a selection id. |
| `meyes_get_view_link` | R I | The shareable URL of what is on screen. |
| `meyes_screenshot` | R | A small image of the scene (longest side ≤ 768 px). |

### 3.4 Built-in analysis (step 6, tier 2, via tab, `viewer:control`)

Every result is a bounded table, and each can also be drawn in the Agent
panel (`show: true`).

| Tool | Does |
|---|---|
| `meyes_list_columns` | obs columns: type, number of categories, whether DE stats exist for it. |
| `meyes_get_markers` | Top genes of one or more groups of a column, from the precomputed DE stats (mean, % expressing, fold change). |
| `meyes_gene_by_group` | Mean and % expressing of ≤ 20 genes across the groups of a column. |
| `meyes_group_composition` | Counts and fractions of a column's groups — overall, within a selection, or split by another column (e.g. per section). |
| `meyes_gene_distribution` | Summary and histogram of one gene, overall or per group. |
| `meyes_coexpression` | Two genes: correlation and the fraction of cells expressing each and both. |
| `meyes_summarize_selection` | Cells, composition and top genes of a selection. |

A DE column listed in `stats.json` `de_stats` is not necessarily a cell
type: the marmoset Xenium stores carry DE stats for QC columns
(`segmentation_method`, `nucleus_count`, …). `meyes_list_columns`
reports every column with its kind; the tool descriptions tell the agent
to group by annotation columns.

### 3.5 Agent compute (step 7, tier 3)

| Tool | Scope | Does |
|---|---|---|
| `meyes_export_table` | `export` | Writes a table the tab assembles (cells × chosen genes and columns, or a group summary) to a file on the person's disk, capped (default 200,000 rows); returns the path, columns and row count — never the rows. |
| `meyes_show_plot` | `viewer:control` | Draws in the Agent panel: a built-in plot by name with its parameters, or a Vega-Lite spec whose data is inline (≤ 5,000 values) or a reference to an exported table. |
| `meyes_show_table` | `viewer:control` | Shows ≤ 500 rows in the Agent panel. |

## 4. The agent API (`/api/agent/v1`)

What the MCP process calls. Versioned in the path (§8), JSON only,
stateless, every response bounded. Authentication is a bearer token
(§6); without one the caller is anonymous.

| Method and path | Tool |
|---|---|
| `GET /catalog/search` | `meyes_search_datasets` |
| `GET /catalog/facets` | `meyes_list_facets` |
| `GET /catalog/{id}` | `meyes_get_dataset` (catalog row) |
| `GET /catalog/{id}/genes` | `meyes_list_genes` |
| `GET /datasets/{id}` | `meyes_get_dataset` (app dataset) |
| `GET /datasets/{id}/genes` | `meyes_list_genes` |
| `GET /me` | `meyes_whoami` |
| `GET /me/datasets` | `meyes_list_my_datasets` |
| `GET /me/datasets/{id}/status` | `meyes_get_processing_status` |
| `GET /projects`, `GET /projects/{id}` | `meyes_list_projects`, `meyes_get_project` |

Errors are `{ "error": code, "message": text }` with the HTTP status of
§6.5. `message` says how to recover.

**Sidecar files.** For an entry whose `s3BaseUrl` is a store
(`…/<name>.zarr`), the folder holding it is read for `catalog.json`,
`stats.json` and `qc.json` ([schemas/](schemas/)); a folder URL is read
as is. Only HTTPS URLs on S3 hosts are fetched (no server-side request to
an arbitrary host); each file is ≤ 256 KB, fetched with a 5 s timeout and
cached for an hour. A missing or invalid file is reported as absent —
legacy and non-BIL folders have none.

## 5. Not in this version

Agent upload (plan P1). Admin actions. Any byte transfer through the agent.

## 6. Authorization: OAuth 2.1

MERFISHEYES is its own authorization server, following the MCP
authorization specification (draft revision read 2026-11-06:
`modelcontextprotocol.io/specification/draft/basic/authorization`). The
site's sign-in (NextAuth: Google, email code) stays the login; the
authorization server adds consent, grants and tokens.

### 6.1 Endpoints and discovery

| Endpoint | Path |
|---|---|
| Issuer | the site origin (`NEXT_PUBLIC_APP_URL`), no path |
| AS metadata (RFC 8414) | `/.well-known/oauth-authorization-server` |
| Protected resource metadata (RFC 9728) | `/.well-known/oauth-protected-resource` and `/.well-known/oauth-protected-resource/api/agent/v1` |
| Authorization | `/oauth/authorize` (a page: sign-in, then consent) |
| Token | `/api/oauth/token` |
| Revocation (RFC 7009) | `/api/oauth/revoke` |
| Dynamic registration (RFC 7591, deprecated in MCP, kept for older clients) | `/api/oauth/register` |

AS metadata advertises: `response_types_supported: ["code"]`,
`grant_types_supported: ["authorization_code", "refresh_token"]`,
`code_challenge_methods_supported: ["S256"]`,
`token_endpoint_auth_methods_supported: ["none"]`,
`client_id_metadata_document_supported: true`,
`authorization_response_iss_parameter_supported: true`,
`scopes_supported`, `registration_endpoint`, `revocation_endpoint`.

### 6.2 Clients

All clients are public (no secret). In the order the spec gives:

1. **Client ID Metadata Documents.** A `client_id` that is an HTTPS URL
   with a path is fetched (§4's limits; no redirects; private and
   loopback addresses refused), and must be JSON whose `client_id` equals
   the URL exactly and that lists `client_name` and `redirect_uris`.
   Cached by HTTP headers, at most a day.
2. **Dynamic registration**: `redirect_uris` (required), `client_name`;
   returns a `client_id`.

Redirect URIs are HTTPS or loopback (`http://127.0.0.1`, `http://[::1]`,
`http://localhost`). Matching is exact, except that a loopback URI may use
any port (RFC 8252 §7.3).

### 6.3 Scopes

| Scope | Grants |
|---|---|
| `catalog:read` | Public catalog. Held by everyone; never needs consent. |
| `projects:read` | The person's projects and datasets in the grant's reach. |
| `projects:write` | Create, edit, delete projects; change membership. |
| `datasets:write` | Edit, overlay, register, submit, withdraw, delete the person's datasets. |
| `viewer:control` | Drive a paired tab. |
| `export` | Write tables to the person's disk. |

`scopes_supported` in the resource metadata is the minimal set
(`catalog:read projects:read`); a tool that needs more answers 403 with
the scope to ask for (§6.5) and the client steps up.

### 6.4 Grants and tokens

- **Consent** shows the client's name, its `client_id` host and the
  redirect host, warns when the redirect is loopback only ("an app on your
  computer"), lists the scopes in words, and asks for the reach: **all my
  projects** or **these projects** (a picker). Both are stored on the grant.
- One **grant** per (person, client): scopes, reach, resource. Consenting
  again updates it.
- **Authorization code**: single use, 5 minutes, bound to the client,
  redirect URI, PKCE challenge (`S256` only) and `resource`. The redirect
  carries `code`, `state` and `iss`.
- **Access token**: opaque, 1 hour, stored as a SHA-256 hash, bound to one
  `resource` (RFC 8707). The API accepts only tokens whose resource is one
  of its own (`/api/agent/v1`, later `/api/mcp`).
- **Refresh token**: opaque, 30 days, rotated on every use. Presenting a
  used refresh token revokes the whole grant (theft detection).
- **Revocation**: the person revokes a connected agent in `/account`; that
  revokes the grant and every token issued under it.

### 6.5 Errors at the API

| Case | Answer |
|---|---|
| No token, or invalid / expired / revoked / wrong-resource token | 401, `WWW-Authenticate: Bearer resource_metadata="…", scope="…"` |
| Valid token, scope missing | 403, `WWW-Authenticate: Bearer error="insufficient_scope", scope="<all scopes this call needs>", resource_metadata="…"` |
| Not visible (§2) | 404 |
| Over a rate limit | 429, `Retry-After: <seconds>`, `{ "error": "rate_limited", "message": … }` |

**Rate limits** (`lib/rate-limit.ts`, counted in Postgres since serverless
instances share nothing; fixed windows): without a token 60 calls a minute
per address; with one 300 a minute per grant and 10,000 a day per person;
`GET /catalog/{id}`, `GET /datasets/{id}` and `GET /me/datasets/{id}/status`
(they read S3 or ask AWS Batch) 30 a minute per caller. On the
authorization server: registration 10 an hour per address, the token
endpoint 30 a minute per client and address, consent decisions 20 an hour
per person.

## 7. The bridge (step 5)

The MCP process listens on `127.0.0.1` (a WebSocket, random port); the tab
connects after the person clicks *Connect agent* and types the 6-digit
code the process printed. Every message is JSON with a `type`.

| Message | From | Fields |
|---|---|---|
| `hello` | tab | `protocol` (`"1.0"`), `app` (site version), `commands` (names it implements, with versions), `panels` (what each shows), `pairing` (the code) |
| `welcome` | process | `protocol`, `session` (id) — or `refuse` with a reason |
| `command` | process | `id`, `name`, `args` |
| `result` | tab | `id`, `ok`, `data` or `error: { code, message }` |
| `event` | tab | `name` (`view_changed`, `selection`, `panel_closed`, `disconnecting`), `data` |
| `ping` / `pong` | either | — |

The process checks the WebSocket `Origin` is the site, refuses a second
tab unless the person allows it, and drops the connection on a wrong
code. The tab shows an "agent connected" indicator with a disconnect
button for as long as the connection is open. JSON Schemas for these
messages are added in step 5 (`schemas/bridge.schema.json`).

## 8. Stability

- **Tools and the agent API**: within a major version, changes are
  additive only — new tools, new optional arguments, new response fields.
  A tool or field to be removed is marked deprecated in its description
  for one minor release first. A breaking change is a new path
  (`/api/agent/v2`) served beside the old.
- **Bridge**: `protocol` is `major.minor`. The tab and the process refuse
  a different major; within a major each side uses only the commands the
  other listed in `hello`.
- **Data the tools read**: `stats.json`, `qc.json`, `catalog.json` follow
  their schemas; `scripts/validate_sidecars.py` is the gate for changing
  what writes them. DE stats and expression are read through the viewer's
  own readers (`ChunkedDataAdapter`, `SpatialDataTableAdapter`,
  `lib/utils/de-stats.ts`), so a format change is made once.
