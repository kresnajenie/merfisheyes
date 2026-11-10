import { defineConfig } from "vitest/config";

// The package's own tests — without this, vitest would find the site's
// config at the repository root.
export default defineConfig({
  test: { include: ["test/**/*.test.ts"], environment: "node" },
});
