// Build the Claude Desktop extension (MCP Bundle, .mcpb) — MCP.md pillar D:
// the package, its production dependencies, the standalone viewer's files and
// the Python pipeline, so a person installs MERFISHEYES into Claude Desktop
// by opening one file. Needs `npm run build` and the viewer built
// (packages/merfisheyes-viewer: next build). Output: build/merfisheyes-<v>.mcpb
import { execFileSync } from "node:child_process";
import { cpSync, existsSync, mkdirSync, readFileSync, rmSync, writeFileSync } from "node:fs";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";

const pkgDir = join(dirname(fileURLToPath(import.meta.url)), "..");
const pkg = JSON.parse(readFileSync(join(pkgDir, "package.json"), "utf8"));
const viewer = join(pkgDir, "..", "merfisheyes-viewer", "out");
const stage = join(pkgDir, "build", "mcpb");
const npm = process.platform === "win32" ? "npm.cmd" : "npm";

if (!existsSync(join(pkgDir, "dist", "index.js"))) throw new Error("run npm run build first");
if (!existsSync(join(viewer, "open", "index.html")))
  throw new Error("build the viewer first: cd packages/merfisheyes-viewer && npx next build");

rmSync(stage, { recursive: true, force: true });
mkdirSync(join(stage, "server"), { recursive: true });
cpSync(join(pkgDir, "dist"), join(stage, "server", "dist"), { recursive: true });
cpSync(join(pkgDir, "overrides"), join(stage, "server", "overrides"), { recursive: true });
cpSync(viewer, join(stage, "server", "viewer"), { recursive: true });
execFileSync(process.execPath, [join(pkgDir, "scripts", "bundle-python.mjs")], { stdio: "inherit" });
cpSync(join(pkgDir, "python"), join(stage, "server", "python"), { recursive: true });
// production dependencies only
writeFileSync(
  join(stage, "server", "package.json"),
  JSON.stringify({ name: pkg.name, version: pkg.version, type: "module", dependencies: pkg.dependencies }, null, 2),
);
execFileSync(npm, ["install", "--omit=dev", "--no-audit", "--no-fund", "--silent"], {
  cwd: join(stage, "server"),
  stdio: "inherit",
});
const manifest = JSON.parse(readFileSync(join(pkgDir, "mcpb", "manifest.json"), "utf8"));

manifest.version = pkg.version;
writeFileSync(join(stage, "manifest.json"), JSON.stringify(manifest, null, 2));
cpSync(join(pkgDir, "mcpb", "icon.png"), join(stage, "icon.png"));
const out = join(pkgDir, "build", `merfisheyes-${pkg.version}.mcpb`);

execFileSync("npx", ["-y", "@anthropic-ai/mcpb@2", "pack", stage, out], { stdio: "inherit" });
console.log(`\n${out}`);
