/**
 * Materialize non-committed test datasets (medium/large tiers) on demand.
 *
 * Reads tests/data/manifest.json and, for each requested dataset whose files are
 * missing, runs its `source`:
 *   - kind "generate": invokes scripts/testdata/generate.py with the recorded args
 *   - kind "download":  fetches source.url, verifies sha256, unpacks .tar.gz/.zip
 *
 * Usage:
 *   npm run test:data:fetch                 # materialize every non-committed dataset
 *   npm run test:data:fetch -- h5ad-medium  # materialize specific id(s)
 *   npm run test:data:fetch -- --tier medium
 */

import { spawnSync } from "node:child_process";
import { createHash } from "node:crypto";
import fs from "node:fs";
import path from "node:path";

const ROOT = path.resolve(__dirname, "..", "..");
const MANIFEST = path.join(ROOT, "tests", "data", "manifest.json");

interface Source {
  kind: "generate" | "download";
  args?: string[];
  outDir?: string;
  url?: string;
  sha256?: string;
}
interface Dataset {
  id: string;
  format: string;
  tier: string;
  path: string;
  committed: boolean;
  source: Source | null;
}

function exists(p: string): boolean {
  return fs.existsSync(path.join(ROOT, p));
}

function sha256(file: string): string {
  return createHash("sha256").update(fs.readFileSync(file)).digest("hex");
}

function generate(ds: Dataset): void {
  const src = ds.source!;
  const outDir = src.outDir ?? ds.path;
  const args = [
    path.join(ROOT, "scripts", "testdata", "generate.py"),
    ...(src.args ?? []),
    "--out",
    path.join(ROOT, outDir),
  ];
  console.log(`  generate ${ds.id}: python3 ${args.join(" ")}`);
  const r = spawnSync("python3", args, { stdio: "inherit" });
  if (r.status !== 0) throw new Error(`generate failed for ${ds.id}`);
}

async function download(ds: Dataset): Promise<void> {
  const src = ds.source!;
  if (!src.url) throw new Error(`${ds.id}: download source missing url`);
  const dest = path.join(ROOT, ds.path);
  fs.mkdirSync(path.dirname(dest), { recursive: true });
  console.log(`  download ${ds.id}: ${src.url}`);
  const res = await fetch(src.url);
  if (!res.ok) throw new Error(`${ds.id}: HTTP ${res.status} fetching ${src.url}`);
  const buf = Buffer.from(await res.arrayBuffer());
  const tmp = dest + ".part";
  fs.writeFileSync(tmp, buf);
  if (src.sha256) {
    const got = sha256(tmp);
    if (got !== src.sha256) {
      fs.unlinkSync(tmp);
      throw new Error(`${ds.id}: sha256 mismatch (expected ${src.sha256}, got ${got})`);
    }
  }
  fs.renameSync(tmp, dest);
}

async function main(): Promise<void> {
  const argv = process.argv.slice(2);
  const tierIdx = argv.indexOf("--tier");
  const tier = tierIdx >= 0 ? argv[tierIdx + 1] : null;
  const ids = argv.filter((a) => !a.startsWith("--") && a !== tier);

  const manifest = JSON.parse(fs.readFileSync(MANIFEST, "utf8"));
  const datasets: Dataset[] = manifest.datasets;

  const targets = datasets.filter((d) => {
    if (d.committed) return false;
    if (ids.length && !ids.includes(d.id)) return false;
    if (tier && d.tier !== tier) return false;
    return true;
  });

  if (!targets.length) {
    console.log("No matching non-committed datasets to materialize.");
    return;
  }

  for (const ds of targets) {
    if (exists(ds.path)) {
      console.log(`  skip ${ds.id}: already present at ${ds.path}`);
      continue;
    }
    if (!ds.source) {
      console.warn(`  WARN ${ds.id}: no source defined; cannot materialize`);
      continue;
    }
    if (ds.source.kind === "generate") generate(ds);
    else await download(ds);
  }
  console.log("Done.");
}

main().catch((err) => {
  console.error(err);
  process.exit(1);
});
