import {
  chmodSync,
  mkdirSync,
  mkdtempSync,
  readFileSync,
  utimesSync,
  writeFileSync,
} from "node:fs";
import { execFileSync } from "node:child_process";
import { createHash } from "node:crypto";
import { createServer } from "node:http";
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

const { SetupNeeded, convert, findUv, setupPipeline, uvAsset, viewable } =
  await import("../src/local/convert.ts");
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

describe("MERFISHEYES's own Python: set up once, only when the person agrees", () => {
  it.skipIf(process.platform === "win32")(
    "asks first, then builds it with uv, then converts",
    async () => {
      const uvLog = join(root, "uv.log");
      const uv = join(root, "fake-uv");

      // venv DIR: a "python" that is node (it runs the stand-in prepare.py)
      writeFileSync(
        uv,
        `#!/bin/sh
echo "$@" >> ${uvLog}
case "$1" in
  --version) echo "uv 0.0.0" ;;
  venv) for d; do dir="$d"; done; mkdir -p "$dir/bin"; ln -sf ${process.execPath} "$dir/bin/python" ;;
  pip) exit 0 ;;
esac
`,
      );
      chmodSync(uv, 0o755);
      const saved = process.env.MERFISHEYES_PYTHON;

      delete process.env.MERFISHEYES_PYTHON;
      process.env.MERFISHEYES_UV = uv;
      try {
        writeFileSync(
          join(scripts, "requirements.txt"),
          "spatialdata==0.0.0\n",
        );
        const input = join(root, "managed.h5ad");

        writeFileSync(input, "x");
        // not set up: no conversion, an explanation of the download instead
        const asked = await convert(input).catch((e) => e);

        expect(asked).toBeInstanceOf(SetupNeeded);
        expect(asked.message).toMatch(/one-time setup.*about 1\.3 GB/);
        // the person agreed
        const lines: string[] = [];
        const python = await setupPipeline((l) => lines.push(l));

        expect(lines.join("\n")).toContain("Setting up MERFISHEYES's Python");
        expect(readFileSync(uvLog, "utf8")).toMatch(/venv .*--python 3\.12/);
        expect(readFileSync(uvLog, "utf8")).toMatch(
          /pip install .* -r .*requirements\.txt/,
        );
        const calls = readFileSync(uvLog, "utf8").split("\n").length;

        expect(await setupPipeline()).toBe(python);
        expect(readFileSync(uvLog, "utf8").split("\n").length).toBe(calls);
        expect(viewable((await convert(input)).folder)).toBe(true);
      } finally {
        process.env.MERFISHEYES_PYTHON = saved;
        delete process.env.MERFISHEYES_UV;
      }
    },
  );

  it("knows which uv to download for each computer", () => {
    expect(uvAsset("darwin", "arm64")).toBe("uv-aarch64-apple-darwin.tar.gz");
    expect(uvAsset("darwin", "x64")).toBe("uv-x86_64-apple-darwin.tar.gz");
    expect(uvAsset("linux", "x64")).toBe("uv-x86_64-unknown-linux-gnu.tar.gz");
    expect(uvAsset("win32", "x64")).toBe("uv-x86_64-pc-windows-msvc.zip");
    expect(uvAsset("linux", "ia32")).toBeNull();
  });

  it.skipIf(process.platform === "win32")(
    "downloads uv only when its checksum matches",
    async (ctx) => {
      const asset = uvAsset()!;
      const pack = join(root, "uvpack");
      const inner = join(pack, asset.replace(/\.tar\.gz$/, ""));

      mkdirSync(inner, { recursive: true });
      writeFileSync(join(inner, "uv"), `#!/bin/sh\necho "uv 9.9.9"\n`);
      chmodSync(join(inner, "uv"), 0o755);
      const tarball = join(root, asset);

      execFileSync("tar", [
        "-czf",
        tarball,
        "-C",
        pack,
        asset.replace(/\.tar\.gz$/, ""),
      ]);
      const bytes = readFileSync(tarball);
      let sum = createHash("sha256").update(bytes).digest("hex");
      const release = createServer((req, res) => {
        if (req.url!.endsWith(".sha256")) res.end(`${sum}  ${asset}\n`);
        else res.end(bytes);
      });

      await new Promise<void>((r) => release.listen(0, "127.0.0.1", r));
      process.env.MERFISHEYES_UV_RELEASE_URL = `http://127.0.0.1:${(release.address() as { port: number }).port}`;
      const saved = { ...process.env };

      // no uv where this test looks (system folders only), so setup downloads one
      process.env.PATH = "/usr/bin:/bin";
      process.env.HOME = join(root, "home");
      // a fresh cache: nothing set up yet
      process.env.MERFISHEYES_CACHE = join(root, "fresh", "converted");
      delete process.env.MERFISHEYES_UV;
      try {
        if (await findUv()) ctx.skip(); // this computer has a uv there
        sum = "0".repeat(64);
        await expect(findUv()).resolves.toBeNull();
        await expect(setupPipeline()).rejects.toThrow(
          /did not match its published checksum/,
        );
        sum = createHash("sha256").update(bytes).digest("hex");
        await expect(setupPipeline()).rejects.not.toThrow(/checksum/);
        // the downloaded uv is now the one MERFISHEYES finds
        expect(await findUv()).toMatch(/uv-0\.\d+\.\d+\/uv$/);
      } finally {
        process.env.PATH = saved.PATH;
        process.env.HOME = saved.HOME;
        process.env.MERFISHEYES_CACHE = saved.MERFISHEYES_CACHE;
        delete process.env.MERFISHEYES_UV_RELEASE_URL;
        release.close();
      }
    },
  );
});

describe("installing on every computer", () => {
  it("adjusts the Linux-locked pins on Macs, nowhere else", async () => {
    const { overridesFile } = await import("../src/local/convert.ts");
    const intel = overridesFile("darwin", "x64")!;
    const silicon = overridesFile("darwin", "arm64")!;

    expect(readFileSync(intel, "utf8")).toMatch(/^llvmlite==0\.45\.1$/m);
    expect(readFileSync(intel, "utf8")).toMatch(/^numba==0\.62\.1$/m);
    expect(readFileSync(intel, "utf8")).toMatch(/^numpy==2\.3\.\d+$/m);
    expect(readFileSync(silicon, "utf8")).toMatch(/^pyproj==3\.6\.1$/m);
    expect(overridesFile("linux", "x64")).toBeNull();
    expect(overridesFile("win32", "x64")).toBeNull();
  });

  it("says a Mac older than macOS 13 cannot convert, before downloading", async () => {
    const { unsupportedMac } = await import("../src/local/convert.ts");

    expect(unsupportedMac("darwin", 21)).toBe(true); // macOS 12
    expect(unsupportedMac("darwin", 22)).toBe(false); // macOS 13
    expect(unsupportedMac("linux", 5)).toBe(false);
  });
});
