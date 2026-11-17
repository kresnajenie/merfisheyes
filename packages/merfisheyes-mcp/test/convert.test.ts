import {
  mkdirSync,
  mkdtempSync,
  readFileSync,
  utimesSync,
  writeFileSync,
} from "node:fs";
import { tmpdir } from "node:os";
import { join } from "node:path";

import { beforeAll, describe, expect, it } from "vitest";

/**
 * D1b: what the viewer cannot read is converted once by the pipeline into
 * the cache, reused while the input is unchanged. A stand-in "Python" (node
 * running a stand-in prepare.py) keeps the test fast and Python-free.
 */
const root = mkdtempSync(join(tmpdir(), "meyes-convert-"));
const scripts = join(root, "scripts");
const runs = join(root, "runs.log");

beforeAll(() => {
  mkdirSync(scripts);
  // prepare.py INPUT OUT: writes OUT/mapping.json + a store; exits 2 for "junk"
  writeFileSync(
    join(scripts, "prepare.py"),
    `const fs = require("fs"), path = require("path");
const [input, out] = process.argv.slice(2);
fs.appendFileSync(${JSON.stringify(runs)}, input + "\\n");
if (input.endsWith("junk.txt")) { console.error("unrecognised"); process.exit(2); }
if (input.endsWith("nopkg.h5ad")) { console.error("ModuleNotFoundError: No module named 'spatialdata'"); process.exit(1); }
fs.mkdirSync(path.join(out, "s.zarr"), { recursive: true });
fs.writeFileSync(path.join(out, "s.zarr", "zarr.json"), "{}");
fs.writeFileSync(path.join(out, "mapping.json"), '{"cells":"s.zarr"}');
console.log("OK done");
`,
  );
  process.env.MERFISHEYES_SCRIPTS_DIR = scripts;
  process.env.MERFISHEYES_PYTHON = process.execPath;
  process.env.MERFISHEYES_CACHE = join(root, "cache");
});

const { convert, viewable } = await import("../src/local/convert.ts");
const { SCRIPTS } = await import("../scripts/bundle-python.mjs");
const conversions = () =>
  readFileSync(runs, "utf8").split("\n").filter(Boolean).length;

describe("opening what the viewer cannot read yet", () => {
  it("knows a store or MERFISHEYES folder from something to convert", () => {
    const store = join(root, "a.zarr");
    const holder = join(root, "holder");

    mkdirSync(store);
    writeFileSync(join(store, "zarr.json"), "{}");
    mkdirSync(join(holder, "x.zarr"), { recursive: true });
    writeFileSync(join(holder, "x.zarr", "zarr.json"), "{}");
    expect(viewable(store)).toBe(true);
    expect(viewable(holder)).toBe(true);
    const file = join(root, "cells.h5ad");

    writeFileSync(file, "h5");
    expect(viewable(file)).toBe(false);
    expect(viewable(root + "/missing")).toBe(false);
  });

  it("converts once into the cache, reuses it, and converts again after a change", async () => {
    const input = join(root, "sample.h5ad");

    writeFileSync(input, "v1");
    const first = await convert(input);

    expect(first.reused).toBe(false);
    expect(first.folder.startsWith(join(root, "cache"))).toBe(true);
    expect(viewable(first.folder)).toBe(true);
    expect(readFileSync(input, "utf8")).toBe("v1");
    expect((await convert(input)).reused).toBe(true);
    const before = conversions();

    writeFileSync(input, "v2, longer");
    utimesSync(input, new Date(), new Date(Date.now() + 5000));
    expect((await convert(input)).reused).toBe(false);
    expect(conversions()).toBe(before + 1);
  });

  it("says what is wrong instead of failing silently", async () => {
    const junk = join(root, "junk.txt");
    const nopkg = join(root, "nopkg.h5ad");

    writeFileSync(junk, "x");
    writeFileSync(nopkg, "x");
    await expect(convert(junk)).rejects.toThrow(
      /does not recognise this as spatial data/,
    );
    await expect(convert(nopkg)).rejects.toThrow(
      /lacks the pipeline's packages \(no module "spatialdata"\)/,
    );
    await expect(convert(join(root, "nope.h5ad"))).rejects.toThrow(
      /does not exist/,
    );
  });

  it("bundles every script prepare.py needs", () => {
    const repo = join(__dirname, "..", "..", "..", "scripts");
    const src = readFileSync(join(repo, "prepare.py"), "utf8");
    const imported = [...src.matchAll(/^(?:from|import) (\w+)/gm)]
      .map((m) => `${m[1]}.py`)
      .filter((f) => {
        try {
          readFileSync(join(repo, f));

          return true;
        } catch {
          return false;
        }
      });

    for (const f of imported) expect(SCRIPTS).toContain(f);
  });
});
