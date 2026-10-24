#!/usr/bin/env node
/**
 * The CHANGELOG.md section of one version, as the notes of its GitHub
 * Release (docs/RELEASING.md).
 *
 *   node scripts/release-notes.mjs [VERSION]   # default: package.json's
 *
 * Prints the section's body (without its `## [X.Y.Z] - date` heading);
 * exits 1 when CHANGELOG.md has no section for that version, so a release
 * whose files disagree stops instead of going out without notes.
 */
import { readFileSync } from "node:fs";
import { fileURLToPath } from "node:url";

/**
 * The body of `## [version]` in a Keep a Changelog file: up to the next
 * version heading or the link references at the end. Null when absent.
 */
export function releaseNotes(changelog, version) {
  const lines = changelog.split("\n");
  const start = lines.findIndex((l) => l.startsWith(`## [${version}]`));

  if (start < 0) return null;
  const rest = lines.slice(start + 1);
  const end = rest.findIndex(
    (l) => l.startsWith("## [") || /^\[[^\]]+\]: /.test(l),
  );

  return rest
    .slice(0, end < 0 ? rest.length : end)
    .join("\n")
    .trim();
}

if (process.argv[1] === fileURLToPath(import.meta.url)) {
  const root = new URL("..", import.meta.url);
  const version =
    process.argv[2] ??
    JSON.parse(readFileSync(new URL("package.json", root), "utf8")).version;
  const notes = releaseNotes(
    readFileSync(new URL("CHANGELOG.md", root), "utf8"),
    version,
  );

  if (!notes) {
    console.error(`CHANGELOG.md has no "## [${version}]" section`);
    process.exit(1);
  }
  console.log(notes);
}
