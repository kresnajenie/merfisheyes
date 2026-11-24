import { readdir, readFile } from "node:fs/promises";
import { join } from "node:path";

import { expect, it } from "vitest";

// Claude Desktop shows the extension's tools from mcpb/manifest.json, not from
// the server: a tool added to the server must be listed there too.
it("the Desktop extension's manifest lists every tool the server registers", async () => {
  const src = join(import.meta.dirname, "..", "src");
  const files = (await readdir(src, { recursive: true })).filter((f) => f.endsWith(".ts"));
  const registered = new Set<string>();

  for (const f of files) {
    const text = await readFile(join(src, f), "utf8");

    for (const m of text.matchAll(/registerTool\(\s*"(meyes_\w+)"/g)) registered.add(m[1]!);
  }
  const manifest = JSON.parse(
    await readFile(join(import.meta.dirname, "..", "mcpb", "manifest.json"), "utf8"),
  ) as { tools: { name: string }[] };

  expect(registered.size).toBeGreaterThan(20);
  expect(manifest.tools.map((t) => t.name).sort()).toEqual([...registered].sort());
});
