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
// MER6-2_E3_1 with a chosen camera.
const SPIRALIA_LANDING =
  "/lm-viewer/from-s3" +
  "?url=https%3A%2F%2Fmerfisheyes-bil.s3.us-west-2.amazonaws.com%2Fyiqun-spiralia%2FMER6-2_E3_1_lm" +
  "&v=eyJjYW0iOlstMjI1LjA0LC04Mi44NSwtMzYyLjc2LDAsMCwwXX0";

const VANITY_LANDINGS = [
  { host: "spiralia.merfisheyes.com", destination: SPIRALIA_LANDING },
  {
    host: "demo1.merfisheyes.com",
    // ACE mouse brain (bil-psc-data2/ace-low-bag) with a saved view: column,
    // overlay genes and right-panel state carried in v / ov / rv.
    destination:
      "/viewer/from-s3" +
      "?url=https%3A%2F%2Fmerfisheyes-bil.s3.us-west-2.amazonaws.com%2Fbil-psc-data2%2Face-low-bag%2Fmeyes_output" +
      "&v=eyJjIjoic3ViY2xhc3NfbmFtZSIsImdzIjpbMCwyXSwic3oiOjIuOH0" +
      "&ov=eyJnZW5lcyI6W1siSWdmYnBsMSIsIiMwMEZGRkYiLDEsdHJ1ZV0sWyJEcmQxIiwiI0ZGMDBGRiIsMSx0cnVlXSxbIlRoIiwiI0ZGRkYwMCIsMSx0cnVlXV0sImdzIjoyLjN9" +
      "&rv=eyJnIjoiQ25yMSIsImMiOiJzdWJjbGFzc19uYW1lIiwiY3QiOlsiMDM3IERHIEdsdXQiLCIwMzggREctUElSIEV4IElNTiJdLCJtIjpbImNlbGx0eXBlIiwiZ2VuZSJdLCJncyI6WzAsMy4zODk4MzMyMTE4OTg4MDM3XX0",
  },
  // Same embryo as spiralia.merfisheyes.com, landed directly so there is no
  // second hop and in-place switching keeps working on this host.
  { host: "demo2.merfisheyes.com", destination: SPIRALIA_LANDING },
  // These two live in other Vercel projects, so they are absolute.
  { host: "demo3.merfisheyes.com", destination: "https://schier.merfisheyes.com" },
  { host: "demo4.merfisheyes.com", destination: "https://heart.merfisheyes.com/3d-heart" },
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
