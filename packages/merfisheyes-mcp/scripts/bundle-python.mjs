// Copy MERFISHEYES's pipeline into the package (python/), in the layout the
// scripts expect (scripts/ beside docs/specs/), so `merfisheyes-mcp view` can
// convert an .h5ad or raw instrument output without a checkout (MCP.md D1b).
// Run on publish; python/ is not committed.
import { cpSync, existsSync, mkdirSync, readdirSync, rmSync } from "node:fs";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";

const repo = join(dirname(fileURLToPath(import.meta.url)), "..", "..", "..");
const out = join(dirname(fileURLToPath(import.meta.url)), "..", "python");
// prepare.py and every script it imports or runs (checked by test/convert.test.ts),
// meyes_analysis.py (the analysis tools, src/local/analysis.ts) and
// raw_merfish.py (meyes_raw_check)
export const SCRIPTS = [
  "combine_slices_v3.py", "de_stats.py", "image_to_ome_zarr.py", "legacy_to_points.py",
  "map_my_cell.py", "meyes_analysis.py", "meyes_common.py", "points_zarr.py", "prepare.py", "prepare_folder.py",
  "prepare_merscope.py", "prepare_visium.py", "prepare_xenium.py", "process_single_molecule.py",
  "raw_merfish.py",
  "process_spatial_data.py", "slices_to_points.py", "spatialdata_cluster.py",
  "spatialdata_points_to_sm.py", "transcripts_to_points.py", "validate_store.py",
  "requirements.txt",
];

if (process.argv[1] === fileURLToPath(import.meta.url)) {
  rmSync(out, { recursive: true, force: true });
  mkdirSync(join(out, "scripts"), { recursive: true });
  for (const f of SCRIPTS) cpSync(join(repo, "scripts", f), join(out, "scripts", f));
  mkdirSync(join(out, "docs", "specs"), { recursive: true });
  cpSync(join(repo, "docs", "specs", "conventions.json"), join(out, "docs", "specs", "conventions.json"));
  cpSync(join(repo, "docs", "specs", "schemas"), join(out, "docs", "specs", "schemas"), { recursive: true });
  if (!existsSync(join(out, "scripts", "prepare.py"))) process.exit(1);
  console.log(`python/: ${readdirSync(join(out, "scripts")).length} scripts + specs`);
}
