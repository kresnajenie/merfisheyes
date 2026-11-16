# @merfisheyes/viewer

The [MERFISHEYES](https://www.merfisheyes.com) viewer as static files — the
site's own viewer pages (cells, molecules, SpatialData stores, images),
built with `next build` (`output: "export"`) into `out/`. Nothing to run on
its own: [`@merfisheyes/mcp`](../merfisheyes-mcp) serves these files from
your computer, beside the folder you open with `merfisheyes-mcp view
<folder>`, so viewing your own data needs neither the internet nor an
account.

Build it from the repository root's dependencies:

```bash
npm ci
cd packages/merfisheyes-viewer
npx next build
```

The pages re-export the site's (`app/*/from-s3/page.tsx`, `app/open/page.tsx`)
inside a bare layout; see `docs/plans/MCP.md`, pillar D2.
