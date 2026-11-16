const path = require("node:path");

/**
 * The standalone viewer (docs/plans/MCP.md pillar D2): the site's own viewer
 * pages and components, exported as static files that @merfisheyes/mcp serves
 * from the person's computer next to their data. No server, no account, no
 * database. Built from the repository root's node_modules and sources (the
 * `@/` alias points there), so it is never a copy of the viewer.
 */
/** @type {import('next').NextConfig} */
module.exports = {
  output: "export",
  trailingSlash: true,
  images: { unoptimized: true },
  turbopack: { root: path.resolve(__dirname, "../..") },
  env: { NEXT_PUBLIC_MERFISHEYES_STANDALONE: "1" },
};
