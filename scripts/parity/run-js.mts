// Run the BROWSER pipelines headlessly in Node and write the same file tree the
// upload would put in S3 (docs/PARITY-TESTING.md).
//
//   npx tsx scripts/parity/run-js.mts <inputs-dir> <out-root> [case-name ...]
//
// Single cell goes through the exact worker parse path (workerApi.parseH5ad /
// parseXenium / parseMerscope → StandardizedDataset.fromSerializedData) and
// then GeneChunkProcessor + createManifest/prepareFilesForUpload — the same
// calls components/upload-settings-modal.tsx makes. Single molecule goes through
// SingleMoleculeDataset.fromParquet/fromCSV + SingleMoleculeProcessor.
import fs from "node:fs";
import path from "node:path";

// ---- Node shims for the browser APIs the pipeline touches -------------------
// PapaParse streams a File through FileReader, which Node does not have.
if (typeof (globalThis as any).FileReader === "undefined") {
  (globalThis as any).FileReader = class {
    result: unknown = null;
    onload: ((e: unknown) => void) | null = null;
    onerror: ((e: unknown) => void) | null = null;
    readAsText(blob: Blob) {
      blob.text().then(
        (t) => { this.result = t; this.onload?.({ target: this }); },
        (e) => this.onerror?.({ target: { error: e } }),
      );
    }
    readAsArrayBuffer(blob: Blob) {
      blob.arrayBuffer().then(
        (b) => { this.result = b; this.onload?.({ target: this }); },
        (e) => this.onerror?.({ target: { error: e } }),
      );
    }
  };
}
// The adapters/processors log a lot; keep it out of the harness output
// (PARITY_DEBUG=1 keeps it, for chasing a divergence).
const say = (s: string) => process.stdout.write(s + "\n");
if (!process.env.PARITY_DEBUG) console.log = () => {};

const [inputsDir, outRoot, ...only] = process.argv.slice(2);
if (!inputsDir || !outRoot) {
  process.stderr.write("usage: run-js.mts <inputs-dir> <out-root> [case ...]\n");
  process.exit(2);
}

type Case = {
  name: string;
  kind: "sc" | "sm";
  format: "h5ad" | "xenium" | "merscope" | "parquet" | "csv";
  input: string;
  datasetType?: "xenium" | "merscope" | "custom";
  cellIdCol?: string | null;
};
const cases: Case[] = JSON.parse(fs.readFileSync(path.join(inputsDir, "cases.json"), "utf8"));

const { workerApi } = await import("@/lib/workers/standardized-dataset.worker");
const { StandardizedDataset } = await import("@/lib/StandardizedDataset");
const { GeneChunkProcessor } = await import("@/lib/utils/GeneChunkProcessor");
const { createManifest, prepareFilesForUpload } = await import("@/lib/utils/sc-upload-files");
const { SingleMoleculeDataset } = await import("@/lib/SingleMoleculeDataset");
const { SingleMoleculeProcessor } = await import("@/lib/utils/SingleMoleculeProcessor");

function fileFrom(absPath: string, relPath?: string, name?: string): File {
  // `name` must be unique per case: H5adAdapter writes the file into h5wasm's
  // virtual FS under file.name, and every case here is called input.h5ad.
  const f = new File([fs.readFileSync(absPath)], name ?? path.basename(absPath));
  // Folder picks in the browser carry webkitRelativePath ("<folder>/<rel>");
  // the adapters match files on it.
  Object.defineProperty(f, "webkitRelativePath", { value: relPath ?? "", configurable: true });
  return f;
}

function folderFiles(dir: string): File[] {
  const base = path.basename(dir);
  const out: File[] = [];
  const walk = (d: string) => {
    for (const e of fs.readdirSync(d, { withFileTypes: true })) {
      const p = path.join(d, e.name);
      if (e.isDirectory()) walk(p);
      else out.push(fileFrom(p, path.posix.join(base, path.relative(dir, p).split(path.sep).join("/"))));
    }
  };
  walk(dir);
  return out;
}

async function writeBlob(root: string, key: string, blob: Blob) {
  const p = path.join(root, key);
  fs.mkdirSync(path.dirname(p), { recursive: true });
  fs.writeFileSync(p, Buffer.from(await blob.arrayBuffer()));
}

async function runSc(c: Case, input: string, out: string) {
  let serial;
  if (c.format === "h5ad") serial = await workerApi.parseH5ad(fileFrom(input, undefined, `${c.name}.h5ad`));
  else if (c.format === "xenium") serial = await workerApi.parseXenium(folderFiles(input));
  else serial = await workerApi.parseMerscope(folderFiles(input));
  const dataset = StandardizedDataset.fromSerializedData(serial);

  const processor = new GeneChunkProcessor();
  const { chunks, index } = await processor.processGenes(dataset);
  const coordinates = await processor.processCoordinates(dataset);
  const observations = await processor.processObservations(dataset);
  const palettes = await processor.processPalettes(dataset);
  const deStatsFiles = await processor.processDeStats(dataset);
  const manifestJson = await createManifest(
    dataset, "parity", c.name, chunks, index, coordinates, observations.metadata, Object.keys(deStatsFiles),
  );
  const files = await prepareFilesForUpload(
    chunks, index, coordinates, observations.files, observations.metadata, palettes, deStatsFiles, manifestJson,
  );
  for (const f of files) await writeBlob(out, f.key, f.blob);
  return files.length;
}

async function runSm(c: Case, input: string, out: string) {
  const file = fileFrom(input);
  const dsType = c.datasetType ?? "xenium";
  const dataset = c.format === "parquet"
    ? await SingleMoleculeDataset.fromParquet(file, dsType)
    : await SingleMoleculeDataset.fromCSV(file, dsType);
  const manifest = SingleMoleculeProcessor.createManifest(dataset, "parity", c.name);
  await writeBlob(out, "manifest.json.gz", await SingleMoleculeProcessor.createManifestBlob(manifest));
  const geneFiles = await SingleMoleculeProcessor.processGenes(dataset);
  for (const [name, blob] of Object.entries(geneFiles)) await writeBlob(out, `genes/${name}.bin.gz`, blob);
  return 1 + Object.keys(geneFiles).length;
}

let failed = 0;
for (const c of cases) {
  if (only.length && !only.includes(c.name)) continue;
  const input = path.join(inputsDir, c.input);
  const out = path.join(outRoot, c.name);
  fs.rmSync(out, { recursive: true, force: true });
  const t0 = Date.now();
  try {
    const n = c.kind === "sc" ? await runSc(c, input, out) : await runSm(c, input, out);
    say(`  js  ${c.name.padEnd(24)} ${n} files  ${((Date.now() - t0) / 1000).toFixed(1)}s`);
  } catch (e) {
    failed++;
    say(`  js  ${c.name.padEnd(24)} FAILED: ${(e as Error)?.stack ?? e}`);
  }
}
process.exit(failed ? 1 : 0);
