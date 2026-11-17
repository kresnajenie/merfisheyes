import { spawn } from "node:child_process";
import { createHash } from "node:crypto";
import {
  existsSync,
  mkdirSync,
  readFileSync,
  readdirSync,
  renameSync,
  rmSync,
  statSync,
  writeFileSync,
} from "node:fs";
import { homedir } from "node:os";
import { basename, join, resolve } from "node:path";
import { fileURLToPath } from "node:url";

/**
 * D1b, first part (MCP.md pillar D): data the viewer cannot read as it is —
 * an .h5ad, a Xenium export, a MERSCOPE region, a molecule table — is
 * converted once by MERFISHEYES's own pipeline (`scripts/prepare.py`, the
 * same the server runs) into a SpatialData store in a cache folder, and that
 * store is opened. The input is never changed. A conversion is reused until
 * the input changes.
 */

/** Can the viewer open this folder as it is? (A store, a folder holding
 * one, a MERFISHEYES folder.) */
export function viewable(path: string): boolean {
  if (!existsSync(path) || !statSync(path).isDirectory()) return false;
  const has = (f: string) => existsSync(join(path, f));

  if (
    has("zarr.json") ||
    has(".zattrs") ||
    has("manifest.json") ||
    has("manifest.json.gz") ||
    has("mapping.json")
  )
    return true;

  return readdirSync(path, { withFileTypes: true }).some(
    (d) =>
      d.isDirectory() &&
      d.name.endsWith(".zarr") &&
      existsSync(join(path, d.name, "zarr.json")),
  );
}

/** Where converted stores are kept: MERFISHEYES_CACHE, else the platform's
 * cache folder. */
export function cacheRoot(): string {
  if (process.env.MERFISHEYES_CACHE)
    return resolve(process.env.MERFISHEYES_CACHE);
  const base =
    process.platform === "win32"
      ? process.env.LOCALAPPDATA || join(homedir(), "AppData", "Local")
      : process.platform === "darwin"
        ? join(homedir(), "Library", "Caches")
        : process.env.XDG_CACHE_HOME || join(homedir(), ".cache");

  return join(base, "merfisheyes", "converted");
}

/** The pipeline's scripts: MERFISHEYES_SCRIPTS_DIR, the copy the npm
 * package carries (`python/scripts`), or the repository's `scripts/`. */
export function findPrepare(): string | null {
  const candidates = [
    process.env.MERFISHEYES_SCRIPTS_DIR &&
      resolve(process.env.MERFISHEYES_SCRIPTS_DIR),
    // src/local or dist/local → the package's python/scripts
    fileURLToPath(new URL("../../python/scripts", import.meta.url)),
    // → the repository's scripts/
    fileURLToPath(new URL("../../../../scripts", import.meta.url)),
  ].filter(Boolean) as string[];
  const dir = candidates.find((d) => existsSync(join(d, "prepare.py")));

  return dir ? join(dir, "prepare.py") : null;
}

/** What identifies the input's current contents: size and modification
 * time of the file, or of a folder and everything one level inside it. */
function fingerprint(path: string): string {
  const st = statSync(path);
  const parts = [`${st.size}:${st.mtimeMs}`];

  if (st.isDirectory())
    for (const name of readdirSync(path).sort()) {
      const s = statSync(join(path, name));

      parts.push(`${name}:${s.size}:${s.mtimeMs}`);
    }

  return createHash("sha256").update(parts.join("|")).digest("hex");
}

export interface Converted {
  /** The folder to open (it holds `mapping.json` and the store). */
  folder: string;
  /** True when an earlier conversion of the same input was reused. */
  reused: boolean;
  seconds: number;
}

const EXIT: Record<number, string> = {
  2: "MERFISHEYES does not recognise this as spatial data: it opens a SpatialData store, an .h5ad (with obsm['spatial']), a Xenium export, a MERSCOPE region, a Visium output or a molecule table (gene, x, y columns)",
  3: "the input is incomplete (a file it needs is missing)",
  4: "the converted store did not pass validation",
};

/**
 * Convert `input` into a store in the cache (or reuse the last conversion of
 * the same, unchanged input). Throws an Error saying what to do when Python
 * or its packages are missing, or the input is not spatial data.
 */
export async function convert(
  input: string,
  opts: { onLine?: (line: string) => void } = {},
): Promise<Converted> {
  const abs = resolve(input);

  if (!existsSync(abs)) throw new Error(`${abs} does not exist.`);
  const id = createHash("sha256").update(abs).digest("hex").slice(0, 10);
  const name = basename(abs).replace(/\.(h5ad|parquet|csv|tsv|zarr)$/i, "");
  const folder = join(cacheRoot(), `${name}-${id}`);
  const stamp = join(folder, ".merfisheyes-source.json");
  const print = fingerprint(abs);

  if (existsSync(stamp)) {
    try {
      const was = JSON.parse(readFileSync(stamp, "utf8"));

      if (was.input === abs && was.fingerprint === print)
        return { folder, reused: true, seconds: 0 };
    } catch {
      // convert again
    }
  }

  const prepare = findPrepare();

  if (!prepare)
    throw new Error(
      "MERFISHEYES's pipeline (prepare.py) is not available here. Set MERFISHEYES_SCRIPTS_DIR to the scripts/ folder of a MERFISHEYES checkout.",
    );
  const python = process.env.MERFISHEYES_PYTHON || "python3";
  const requirements = join(prepare, "..", "requirements.txt");
  const work = `${folder}.tmp-${process.pid}`;

  rmSync(work, { recursive: true, force: true });
  mkdirSync(cacheRoot(), { recursive: true });
  const started = Date.now();
  const output: string[] = [];
  const code = await new Promise<number>((done, fail) => {
    const child = spawn(python, [prepare, abs, work], {
      stdio: ["ignore", "pipe", "pipe"],
    });
    const take = (b: Buffer) =>
      b
        .toString()
        .split("\n")
        .filter(Boolean)
        .forEach((l) => {
          output.push(l);
          opts.onLine?.(l);
        });

    child.stdout.on("data", take);
    child.stderr.on("data", take);
    child.on("error", (e: NodeJS.ErrnoException) =>
      fail(
        e.code === "ENOENT"
          ? new Error(
              `No Python found ("${python}"). Install Python 3.12 and the pipeline's packages (pip install -r ${requirements}), or set MERFISHEYES_PYTHON to a Python that has them.`,
            )
          : e,
      ),
    );
    child.on("close", (c) => done(c ?? 1));
  });

  if (code !== 0) {
    rmSync(work, { recursive: true, force: true });
    const missing = output
      .map((l) => /No module named '([^']+)'/.exec(l)?.[1])
      .find(Boolean);

    if (missing)
      throw new Error(
        `The Python at "${python}" lacks the pipeline's packages (no module "${missing}"). Install them once: pip install -r ${requirements} — or set MERFISHEYES_PYTHON to a Python that has them.`,
      );
    throw new Error(
      `Could not convert ${abs}: ${EXIT[code] ?? `prepare.py stopped (exit ${code})`}.\n${output.slice(-8).join("\n")}`,
    );
  }
  writeFileSync(
    join(work, ".merfisheyes-source.json"),
    JSON.stringify(
      {
        input: abs,
        fingerprint: print,
        converted_at: new Date().toISOString(),
      },
      null,
      2,
    ),
  );
  rmSync(folder, { recursive: true, force: true });
  renameSync(work, folder);

  return { folder, reused: false, seconds: (Date.now() - started) / 1000 };
}
