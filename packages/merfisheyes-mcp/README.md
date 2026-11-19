# @merfisheyes/mcp

An [MCP](https://modelcontextprotocol.io) server for
[MERFISHEYES](https://www.merfisheyes.com), the web viewer for spatial
transcriptomics. It lets an AI agent — Claude Code, Claude Desktop, or any
MCP client — search the public catalog of datasets and, once you sign in,
read your own projects and datasets.

It runs on your computer (stdio). MERFISHEYES does the work; the agent asks
and gets short answers back.

## Use it

You need Node.js 20.3 or newer. Check with:

```bash
node -v
```

Claude Code, from any folder:

```bash
claude mcp add --scope user merfisheyes -- npx -y @merfisheyes/mcp@latest
```

`--scope user` makes the tools available in every Claude Code session, not
only in the folder you ran it from; `@latest` picks up new versions when
they are published. Start a new Claude Code session and type `/mcp`:
`merfisheyes` is listed with its tools.

Claude Desktop, Cursor and other clients — add to the client's MCP settings:

```json
{
  "mcpServers": {
    "merfisheyes": {
      "command": "npx",
      "args": ["-y", "@merfisheyes/mcp@latest"]
    }
  }
}
```

To run your own changes instead, from a checkout of the
[MERFISHEYES repository](https://github.com/kresnajenie/merfisheyes), in
`packages/merfisheyes-mcp`:

```bash
npm ci
npm run build
claude mcp add --scope user merfisheyes -- node "$(pwd)/dist/index.js"
```

## Your own data, without uploading it

Install once, then open anything — a SpatialData store, an `.h5ad`, a
Xenium export, a MERSCOPE region, a molecule table — in the MERFISHEYES
viewer on your computer:

```bash
npm install -g @merfisheyes/mcp
merfisheyes view my_sample.h5ad
```

(Or without installing: `npx -y @merfisheyes/mcp@latest view my_sample.h5ad`.)

Or ask the agent ("show my sample in MERFISHEYES"): it calls
`meyes_open_local`. The folder is served from your computer on
`127.0.0.1` — only to the MERFISHEYES site, only under a random address,
read-only, only while the command or the agent runs — and the viewer in
your browser reads it from there. Nothing is uploaded, and the folder's
address is never sent to the site.

It opens a SpatialData store (`.zarr`) or a folder holding one, a folder
MERFISHEYES processed (`manifest.json`), or an OME-Zarr image as they are.
An `.h5ad` (with `obsm["spatial"]`), a Xenium export, a MERSCOPE region,
a Visium output or a molecule table (`.parquet` / `.csv` with gene, x, y)
is first converted once by MERFISHEYES's own pipeline into a SpatialData
store in a cache folder (`~/.cache/merfisheyes/converted`, or
`MERFISHEYES_CACHE`), and reused until the input changes; the input is
never changed. Converting runs MERFISHEYES's own Python pipeline, which
needs a one-time setup: Python 3.12 and the pipeline's pinned packages
(about 1.3 GB), downloaded into MERFISHEYES's own folder, plus
[uv](https://docs.astral.sh/uv/) (about 15 MB, from its official release,
checksum verified) when the computer has none. Nothing is installed
system-wide, and it never happens without asking: `merfisheyes view`
asks "Set it up now?", an agent asks you before calling its setup tool,
and `merfisheyes setup` does it directly. Viewing data that is already a
store needs none of it. To use a Python you already have instead, set
`MERFISHEYES_PYTHON`.

The viewer itself comes with it (`@merfisheyes/viewer`) and is served from
your computer too, so nothing goes to the internet. Without it, links open
the site's viewer instead, and in Chrome you allow www.merfisheyes.com to
reach apps on this device the first time.

## Signing in

Public catalog tools work without an account. For your own data, ask the
agent to sign in (or run `npx -y @merfisheyes/mcp@latest login` in a terminal): a
MERFISHEYES page opens, you sign in as usual, and choose what the agent may
reach — **all your projects or only the ones you pick**. Disconnect it any
time from your account page on MERFISHEYES, or with `merfisheyes-mcp logout`.

Tokens are kept in `~/.config/merfisheyes-mcp/credentials.json`
(`%APPDATA%\merfisheyes-mcp\` on Windows), readable only by you.

## Tools

| Tool | Needs sign-in | Does |
|---|---|---|
| `meyes_find_datasets` | no | Search every public dataset (~12,700, MERFISHEYES's and SpatialData-DB's) by organism, tissue or disease with everything below it, technology, age, sex, genotype, genes measured or expressed — answered as a count, a breakdown and the first rows |
| `meyes_index_values` | no | The values a field takes (tissues, technologies, diseases…), with counts, optionally among filtered datasets |
| `meyes_search_datasets` | no | Search the catalog: text, species, tissue, platform, genes measured |
| `meyes_list_facets` | no | The species, tissue and platform spellings filters use |
| `meyes_get_dataset` | no¹ | One dataset: viewer links, cells, genes, molecules, sections, cell-type levels, DE columns, QC |
| `meyes_list_genes` | no¹ | A dataset's genes, filtered by a name fragment |
| `meyes_open_view` | no¹ | A link to a dataset in the viewer already showing what you asked: genes coloured, a cell-type column, 2D or 3D, dot size, opacity, or two datasets side by side in split screen. Opens a tab only when you ask |
| `meyes_open_local` | no | A folder on your computer in the viewer, served from your machine (never uploaded) |
| `meyes_setup_converter` | no | The one-time converter setup, only after you agree |
| `meyes_view_state` | no | What the open local view shows: colouring, hidden values, columns and their values |
| `meyes_view_set` | no | Change the open local view: colour by a gene or column, show only / hide values (a bad FOV), dot size, 2D/3D |
| `meyes_view_snapshot` | no | A PNG of the open local view with its legend, saved and shown to the agent |
| `meyes_view_add_column` | no | Add the agent's own results (pseudotime, labels, scores) to the open view from a CSV keyed by cell id, without changing your data |
| `meyes_sign_in` / `meyes_sign_out` | — | Connect or disconnect your account |
| `meyes_whoami` | — | Who is signed in, and what the agent may do |
| `meyes_list_my_datasets` | yes | Your datasets and their processing status |
| `meyes_get_processing_status` | yes | Where one of your uploads is in processing |
| `meyes_list_projects` / `meyes_get_project` | yes | Your projects and their datasets |

¹ Your own datasets open too once you are signed in.

Every tool answers in markdown, or JSON with `response_format: "json"`; lists
are paged (`limit`, `offset`, `next_offset`).

## Settings

| Variable | Default | |
|---|---|---|
| `MERFISHEYES_URL` | `https://www.merfisheyes.com` | Another deployment, e.g. `https://dev.merfisheyes.com` (staging) or `http://localhost:3000` |
| `MERFISHEYES_CREDENTIALS` | see above | Where tokens are kept |
| `MERFISHEYES_NO_BROWSER` | unset | Never open a browser; only give links |
| `MERFISHEYES_LOCAL_PORT` | a free port | The port `view` and `meyes_open_local` serve folders on |
| `MERFISHEYES_PYTHON` | its own, after the one-time setup | The Python that converts `.h5ad` and raw instrument output |
| `MERFISHEYES_UV` | found, or downloaded by the setup | The uv that builds that Python |
| `MERFISHEYES_CACHE` | the platform's cache folder | Where converted stores are kept |
| `MERFISHEYES_VIEWER` | unset | `site`: open the site's viewer even when the standalone one is installed |
| `MERFISHEYES_SIGN_IN_WAIT_MS` | `120000` | How long sign-in waits for you in the browser before handing back the link |
| `MERFISHEYES_VERCEL_BYPASS` | unset | Only for a deployment behind Vercel's Deployment Protection (staging): its *Protection Bypass for Automation* secret, sent with every request. The package then registers itself instead of using the site's client document. Keep it on your machine |

## Develop

```bash
npm ci
npm run build
npm test
npx @modelcontextprotocol/inspector node dist/index.js
```

`build` compiles strict TypeScript into `dist/`; `test` runs the unit tests
(sign-in against a stand-in server, paging, viewer links); the inspector
lets you call each tool by hand.

The end-to-end test — this package over stdio against the site, signing in
through a browser — is `e2e/tests/upload/agent-mcp.spec.ts` in the
repository. The protocol it follows is `docs/specs/MCP-PROTOCOL.md`.
