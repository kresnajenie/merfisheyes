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

Claude Code:

```bash
claude mcp add merfisheyes -- npx -y @merfisheyes/mcp
```

Claude Desktop, Cursor and other clients — add to the client's MCP settings:

```json
{
  "mcpServers": {
    "merfisheyes": { "command": "npx", "args": ["-y", "@merfisheyes/mcp"] }
  }
}
```

Until the package is published to npm, run it from a checkout of the
MERFISHEYES repository. In `packages/merfisheyes-mcp`:

```bash
npm ci
npm run build
claude mcp add merfisheyes -- node "$(pwd)/dist/index.js"
```

## Signing in

Public catalog tools work without an account. For your own data, ask the
agent to sign in (or run `npx @merfisheyes/mcp login` in a terminal): a
MERFISHEYES page opens, you sign in as usual, and choose what the agent may
reach — **all your projects or only the ones you pick**. Disconnect it any
time from your account page on MERFISHEYES, or with `merfisheyes-mcp logout`.

Tokens are kept in `~/.config/merfisheyes-mcp/credentials.json`
(`%APPDATA%\merfisheyes-mcp\` on Windows), readable only by you.

## Tools

| Tool | Needs sign-in | Does |
|---|---|---|
| `meyes_search_datasets` | no | Search the catalog: text, species, tissue, platform, genes measured |
| `meyes_list_facets` | no | The species, tissue and platform spellings filters use |
| `meyes_get_dataset` | no¹ | One dataset: viewer links, cells, genes, molecules, sections, cell-type levels, DE columns, QC |
| `meyes_list_genes` | no¹ | A dataset's genes, filtered by a name fragment |
| `meyes_open_view` | no¹ | Open a dataset in your browser already showing what you asked: genes coloured, a cell-type column, 2D or 3D |
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
| `MERFISHEYES_NO_BROWSER` | unset | Never open a browser; only print the sign-in link |
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
