/**
 * Vanity subdomains that land on one dataset.
 *
 * Only "/" is special-cased. Everything else stays a normal route on the same
 * host, which is what makes the viewer's in-place dataset switching keep
 * working: it rewrites the address to /lm-viewer/from-s3?url=… via pushState,
 * and that path resolves on this host exactly as it does on the main one.
 *
 * Destinations must be fully percent-encoded. Next parses ":" in a redirect
 * destination as a path parameter, so a bare "https://" would be mangled.
 */
const VANITY_LANDINGS = [
  {
    host: "spiralia.merfisheyes.com",
    // MER6-2_E3_1 with a chosen camera.
    destination:
      "/lm-viewer/from-s3" +
      "?url=https%3A%2F%2Fmerfisheyes-bil.s3.us-west-2.amazonaws.com%2Fyiqun-spiralia%2FMER6-2_E3_1_lm" +
      "&v=eyJjYW0iOlstMjI1LjA0LC04Mi44NSwtMzYyLjc2LDAsMCwwXX0",
  },
];

/** @type {import('next').NextConfig} */
const nextConfig = {
  // Allow the E2E server to use a separate build dir so it never collides with a
  // dev server already running against .next (set NEXT_DIST_DIR=.next-e2e).
  ...(process.env.NEXT_DIST_DIR ? { distDir: process.env.NEXT_DIST_DIR } : {}),
  images: {
    remotePatterns: [{ hostname: "lh3.googleusercontent.com" }],
  },
  async redirects() {
    return VANITY_LANDINGS.map(({ host, destination }) => ({
      source: "/",
      has: [{ type: "host", value: host }],
      destination,
      // 307, not 308. A permanent redirect is cached by the browser, which
      // would pin this subdomain to one embryo on every machine that had ever
      // visited it, until each user cleared their cache.
      permanent: false,
    }));
  },
};

module.exports = nextConfig;
