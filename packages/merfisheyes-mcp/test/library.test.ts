import { cp, mkdtemp, readFile, writeFile } from "node:fs/promises";
import { tmpdir } from "node:os";
import { join } from "node:path";

import { beforeEach, expect, it, vi } from "vitest";

/**
 * The library search (MCP.md pillar B) over a 63-row index cut from the
 * published one by scripts/make-index-fixture.py. Expected numbers were
 * counted from the fixture files with pandas, not with this code.
 */

const FIXTURE = join(import.meta.dirname, "fixtures", "index");

process.env.MERFISHEYES_INDEX_URL = FIXTURE;
const { libraryIndex, resetLibraryIndex } = await import(
  "../src/library/data.ts"
);
const { search, values, hit, UnknownValue } = await import(
  "../src/library/search.ts"
);

beforeEach(() => resetLibraryIndex());

const find = async (f: Parameters<typeof search>[1]) =>
  search(await libraryIndex(FIXTURE), { collapse_rereleases: false, ...f });

it("loads the index once per release", async () => {
  const a = await libraryIndex(FIXTURE);
  const b = await libraryIndex(FIXTURE);

  expect(a).toBe(b);
  expect(a.release).toBe("fixture");
  expect(a.rows).toHaveLength(63);
  expect(typeof a.rows[0].n_obs).toBe("number"); // Int64 is not left a BigInt
});

it("finds a tissue with everything below it", async () => {
  expect((await find({ organism: ["mouse"] })).total).toBe(45);
  const brain = await find({ organism: ["mouse"], tissue: ["brain"] });

  expect(brain.total).toBe(35);
  const hip = await find({ tissue: ["hippocampal formation"] });

  expect(hip.total).toBe(2);
  // By id the same as by name.
  expect(
    (await find({ organism: ["mouse"], tissue: ["UBERON:0000955"] })).total,
  ).toBe(35);
});

it("finds a disease class by its grouping", async () => {
  expect((await find({ organism: ["human"], disease: ["cancer"] })).total).toBe(
    10,
  );
});

it("reads organism aliases and technology names", async () => {
  expect((await find({ organism: ["Mus musculus"] })).total).toBe(45);
  expect((await find({ organism: ["mice"] })).total).toBe(45);
  expect((await find({ technology: ["MERSCOPE"] })).total).toBe(22);
  expect((await find({ technology: ["vizgen_merscope"] })).total).toBe(22);
});

it("names the closest values when one is unknown", async () => {
  await expect(find({ tissue: ["hipocampus"] })).rejects.toThrow(UnknownValue);
  await expect(find({ tissue: ["hipocampus"] })).rejects.toThrow(
    /Closest: .*hippocampal formation/,
  );
  await expect(find({ organism: ["mouze"] })).rejects.toThrow(/Closest: mouse/);
});

it("says how many rows record no age, and leaves them out", async () => {
  const r = await find({
    organism: ["mouse"],
    age_days_min: 56,
    age_days_max: 200,
  });

  expect(r.total).toBe(9);
  expect(r.notes.join(" ")).toMatch(/^21 datasets .* no age in days/m);
});

it("finds datasets measuring genes over every row, from the lookup", async () => {
  const one = await find({ genes: [{ gene: "Gad1", min_fraction: 0 }] });

  expect(one.total).toBe(37);
  expect(
    (
      await find({
        organism: ["mouse"],
        genes: [{ gene: "GAD1", min_fraction: 0 }],
      })
    ).total,
  ).toBe(26);
  const both = await find({
    organism: ["mouse"],
    genes: [
      { gene: "Gad1", min_fraction: 0 },
      { gene: "Slc17a7", min_fraction: 0 },
    ],
  });

  expect(both.total).toBe(25);
  expect(
    (await find({ genes: [{ gene: "notagene", min_fraction: 0 }] })).total,
  ).toBe(0);
});

it("filters on expression where counts exist, and says what it could not check", async () => {
  const r = await find({ genes: [{ gene: "Rab3b", min_fraction: 0.05 }] });

  // Ours and SpatialData-DB's, from gene_expression.parquet's row group.
  expect(r.total).toBe(36);
  expect(r.breakdown.meyes_source).toEqual([
    { value: "sddb", count: 21 },
    { value: "bil", count: 15 },
  ]);
  expect(r.notes.join(" ")).toMatch(
    /1 datasets that measure it have no per-dataset counts/,
  );
  for (const row of r.rows)
    expect(
      hit(row, r).gene_stats?.[0].fraction_detected,
    ).toBeGreaterThanOrEqual(0.05);
  expect(
    (await find({ genes: [{ gene: "Gad1", min_fraction: 0.05 }] })).total,
  ).toBe(25);
});

it("folds re-releases into the newest release", async () => {
  const index = await libraryIndex(FIXTURE);
  const all = await search(index, { collapse_rereleases: false });
  const folded = await search(index, {});

  // Groups of 2, 2 and 3 in the fixture (one more group has a single row).
  expect(all.total - folded.total).toBe(4);
  const ids = new Set(folded.rows.map((r) => r.meyes_id));

  expect(ids.has("cmviigvvb00297aoqacdvfw7u")).toBe(true); // Space Ranger 4.0.1
  expect(ids.has("cmviihdl400eo7aoqfui61p5k")).toBe(false); // 3.0.0, published later
  expect(ids.has("cmviih4b900867aoqdski27py")).toBe(true); // 1.1.0 over 1.0.0
  const h = hit(index.byMeyesId.get("cmviigxzr003r7aoq65eddfka")!, folded);

  expect(h.also_published_as?.map((o) => o.release)).toEqual([
    "Space Ranger 1.0.0",
    "Space Ranger 1.1.0",
  ]);
});

it("breaks results down and lists a field's values", async () => {
  const index = await libraryIndex(FIXTURE);
  const r = await search(index, { organism: ["mouse"] });

  expect(r.breakdown.organism).toEqual([{ value: "mouse", count: r.total }]);
  expect(values(index.rows, "meyes_sex")).toEqual([
    { value: "male", count: 18 },
    { value: "female", count: 3 },
  ]);
});

it("refuses gene bitsets when the catalog and lookup are from different builds", async () => {
  const dir = await mkdtemp(join(tmpdir(), "meyes-index-"));

  await cp(FIXTURE, dir, { recursive: true });
  const meta = JSON.parse(
    await readFile(join(dir, "gene_search_meta.json"), "utf8"),
  );

  meta.meyes_ids.reverse();
  await writeFile(join(dir, "gene_search_meta.json"), JSON.stringify(meta));
  const index = await libraryIndex(dir);

  await expect(
    search(index, { genes: [{ gene: "Gad1", min_fraction: 0 }] }),
  ).rejects.toThrow(/different builds/);
  // …and drops the index, so the next call loads one build whole.
  expect(await libraryIndex(dir)).not.toBe(index);
});

it("drops the index when a gene file no longer reads (a rebuild mid-publish)", async () => {
  const dir = await mkdtemp(join(tmpdir(), "meyes-index-"));

  await cp(FIXTURE, dir, { recursive: true });
  const index = await libraryIndex(dir);

  await writeFile(join(dir, "gene_expression.parquet"), "not parquet");
  await expect(
    search(index, { genes: [{ gene: "Rab3b", min_fraction: 0.05 }] }),
  ).rejects.toThrow(/may be being updated/);
  expect(await libraryIndex(dir)).not.toBe(index);
});

it("reloads a release rebuilt the same day, and only then", async () => {
  const dir = await mkdtemp(join(tmpdir(), "meyes-index-"));

  await cp(FIXTURE, dir, { recursive: true });
  vi.useFakeTimers({ toFake: ["Date"] });
  try {
    const first = await libraryIndex(dir);

    vi.setSystemTime(Date.now() + 11 * 60 * 1000); // past the recheck
    expect(await libraryIndex(dir)).toBe(first); // same build: kept
    const info = JSON.parse(
      await readFile(join(dir, "release_info.json"), "utf8"),
    );

    await writeFile(
      join(dir, "release_info.json"),
      JSON.stringify({ ...info, built_at: "2026-11-18T23:00:00+00:00" }),
    );
    vi.setSystemTime(Date.now() + 11 * 60 * 1000);
    const second = await libraryIndex(dir);

    expect(second).not.toBe(first); // same release name, new build: reloaded
    expect(second.builtAt).toBe("2026-11-18T23:00:00+00:00");
  } finally {
    vi.useRealTimers();
  }
});

it("answers through the MCP tools", async () => {
  const { McpServer, InMemoryTransport } = await import(
    "@modelcontextprotocol/server"
  );
  const { Client } = await import("@modelcontextprotocol/client");
  const { registerLibraryTools } = await import("../src/tools/library.ts");
  const server = new McpServer({ name: "t", version: "0" });

  registerLibraryTools(server);
  const [a, b] = InMemoryTransport.createLinkedPair();
  const client = new Client({ name: "c", version: "0" });

  await Promise.all([server.connect(a), client.connect(b)]);
  const md = await client.callTool({
    name: "meyes_find_datasets",
    arguments: {
      organism: "mouse",
      tissue: "brain",
      genes_measured: ["Gad1"],
      limit: 2,
    },
  });
  const body = (md.content as { text: string }[])[0].text;

  expect(body).toMatch(/^# \d+ datasets/);
  expect(body).toMatch(/## Breakdown/);
  expect(body).toMatch(/More: call again with offset=2/);
  const js = await client.callTool({
    name: "meyes_find_datasets",
    arguments: { technology: ["MERSCOPE"], response_format: "json", limit: 3 },
  });
  const page = js.structuredContent as {
    total: number;
    items: unknown[];
    next_offset: number;
  };

  expect(page.total).toBeGreaterThan(0);
  expect(page.items).toHaveLength(3);
  expect(page.next_offset).toBe(3);
  const bad = await client.callTool({
    name: "meyes_find_datasets",
    arguments: { tissue: "hipocampus" },
  });

  expect(bad.isError).toBe(true);
  const vals = await client.callTool({
    name: "meyes_index_values",
    arguments: { column: "technology", organism: "mouse" },
  });

  expect((vals.content as { text: string }[])[0].text).toMatch(
    /vizgen_merscope/,
  );
  await client.close();
});
