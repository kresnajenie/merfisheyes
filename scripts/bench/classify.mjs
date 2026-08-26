// The browser's folder filter, mirrored for the benchmark harness.
//
// Kept in lockstep with lib/ingest/classify-folder.ts: a folder export is
// mostly files neither pipeline opens, and both bench paths must feed the app
// exactly what the real upload would.
import fs from "node:fs";
import path from "node:path";

const XENIUM_EXPRESSION = new Set([
  "cell_by_gene.csv", "cell_feature_matrix.h5", "features.tsv", "features.tsv.gz",
  "barcodes.tsv", "barcodes.tsv.gz", "matrix.mtx", "matrix.mtx.gz",
]);

const isTop = (key) => !key.includes("/");
const baseName = (key) => key.split("/").pop().toLowerCase();

/** Every file under `dir`, keyed by its path relative to `dir`. */
export function walk(dir, prefix = "") {
  const out = [];
  let entries;
  try {
    entries = fs.readdirSync(dir, { withFileTypes: true });
  } catch {
    return out;
  }
  for (const e of entries) {
    if (e.name.startsWith(".")) continue;
    const abs = path.join(dir, e.name);
    const key = prefix ? `${prefix}/${e.name}` : e.name;
    let st;
    try { st = fs.statSync(abs); } catch { continue; } // follows symlinks
    if (st.isDirectory()) out.push(...walk(abs, key));
    else out.push({ key, abs, size: st.size });
  }
  return out;
}

/**
 * The files the app would actually consume from a dataset directory.
 * Returns [{ key, abs, size }] with `key` relative to the dataset root, so the
 * Xenium `cell_feature_matrix/` nesting the loader looks for is preserved.
 */
export function classifyDir(dir, format) {
  const all = walk(dir);

  if (format === "xenium") {
    const cells = all.filter(
      (f) => isTop(f.key) && ["cells.csv", "cells.csv.gz"].includes(baseName(f.key)),
    );
    const expression = all.filter(
      (f) => XENIUM_EXPRESSION.has(baseName(f.key)) &&
        (isTop(f.key) || f.key.toLowerCase().startsWith("cell_feature_matrix/")),
    );
    return [...cells, ...expression];
  }

  if (format === "merscope") {
    return all.filter((f) => {
      const n = baseName(f.key);
      if (!isTop(f.key) || !n.endsWith(".csv")) return false;
      const isMetadata = n.includes("metadata") && !n.includes("gene");
      const isCellByGene = n.includes("cell") && n.includes("gene");
      return isMetadata || isCellByGene;
    });
  }

  return all;
}

/**
 * A directory containing only `files`, as symlinks, with their relative paths
 * intact — what to hand a folder file-input.
 *
 * Playwright's setInputFiles(directory) streams every file in the tree over
 * CDP; a real browser hands over lazy File handles and reads only what the
 * classifier picks. Staging keeps the harness faithful to the browser (and
 * stops a 2 GB transcripts file being transferred for a single-cell run).
 */
export function stageDir(files, dir) {
  fs.rmSync(dir, { recursive: true, force: true });
  for (const f of files) {
    const dest = path.join(dir, f.key);
    const src = fs.realpathSync(f.abs);
    fs.mkdirSync(path.dirname(dest), { recursive: true });
    // Hardlink, never symlink: Chromium's webkitdirectory enumeration stalls on
    // symlinked entries (a 168 MB MERSCOPE folder of symlinks never finished
    // uploading, while the same folder of real files went through in seconds).
    // Falls back to a copy when the stage dir is on another filesystem.
    try {
      fs.linkSync(src, dest);
    } catch {
      fs.copyFileSync(src, dest);
    }
  }
  return dir;
}
