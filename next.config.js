/** @type {import('next').NextConfig} */
const nextConfig = {
  // Allow the E2E server to use a separate build dir so it never collides with a
  // dev server already running against .next (set NEXT_DIST_DIR=.next-e2e).
  ...(process.env.NEXT_DIST_DIR ? { distDir: process.env.NEXT_DIST_DIR } : {}),
  images: {
    remotePatterns: [
      { hostname: "lh3.googleusercontent.com" },
    ],
  },
};

module.exports = nextConfig;
