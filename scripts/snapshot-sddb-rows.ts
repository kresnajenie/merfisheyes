/**
 * Snapshot every SpatialData-DB catalog row (source = "sddb", with its
 * entries) to a JSON file before a sync touches them. Read-only.
 *
 *   npx tsx scripts/snapshot-sddb-rows.ts <out.json> --env=.env
 *   npx tsx scripts/snapshot-sddb-rows.ts <out.json> --env=/path/to/.env.prod --direct
 *
 * The env file is read literally, like scripts/sync-sddb-catalog.ts: never
 * `source` .env.prod in a shell (its password holds characters the shell
 * rewrites). --direct uses the file's DIRECT_URL. Keep production snapshots
 * out of the repo (e.g. `.cache/sddb-snapshots/`, git-ignored).
 * Spec: docs/specs/sddb-catalog.md ("Bumping to a new release").
 */
import { readFileSync, writeFileSync } from "fs";

const arg = (name: string) =>
  process.argv.find((a) => a.startsWith(`--${name}=`))?.slice(name.length + 3);
const out = process.argv[2];
const envFile = arg("env");

if (!out || out.startsWith("--")) {
  throw new Error("usage: snapshot-sddb-rows.ts <out.json> [--env=<file>] [--direct]");
}
if (envFile) {
  for (const line of readFileSync(envFile, "utf8").split("\n")) {
    const m = /^([A-Z_]+)=(.*)$/.exec(line);

    if (m) process.env[m[1]] = m[2].replace(/^"|"$/g, "");
  }
  if (process.argv.includes("--direct") && process.env.DIRECT_URL) {
    process.env.DATABASE_URL = process.env.DIRECT_URL;
  }
}

async function main() {
  // The database only now: the env is loaded, so the client sees DATABASE_URL.
  const { prisma } = await import("../lib/prisma");
  const rows = await prisma.catalogDataset.findMany({
    where: { source: "sddb" },
    include: { entries: true },
  });

  writeFileSync(out, JSON.stringify({ at: new Date().toISOString(), rows }, null, 1));
  const releases: Record<string, number> = {};

  for (const r of rows) {
    releases[r.sourceRelease ?? "?"] = (releases[r.sourceRelease ?? "?"] ?? 0) + 1;
  }
  console.log(
    `${rows.length} sddb rows (${rows.filter((r) => r.isPublished).length} published) → ${out}`,
    releases,
  );
  await prisma.$disconnect();
}

main().catch((e) => {
  console.error(e);
  process.exit(1);
});
