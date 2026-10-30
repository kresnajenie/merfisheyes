import fs from "node:fs";
import path from "node:path";

/**
 * What a person sees of a dataset, read the same way from both sides of a
 * save (docs/plans/PREVIEW-SAVE.md, D5):
 * - `previewOf`: the browser preview's own parse (the worker's parseH5ad /
 *   parseMerscope → StandardizedDataset), as a dropped file is opened
 * - `storeOf`: the store scripts/prepare.py wrote, opened the way the
 *   SpatialData viewer opens a dropped or saved one
 *   (openSpatialDataStoreFromFiles → loadSpatialDataCells)
 * Both end in a StandardizedDataset and are described by the one `describe`,
 * so a difference is in the data, never in how it was read.
 */

export interface Seen {
  cells: number;
  genes: string[];
  /** One [x, y] per cell, in order; NaN for a cell with no position. */
  positions: number[][];
  embeddings: string[];
  /** Every column offered for colouring, with its kind. */
  columns: Record<string, string>;
  /** The column the viewer opens on. */
  defaultColumn: string | null;
  /** Per categorical column, each label's colour. */
  colours: Record<string, Record<string, string>>;
  /** A gene's value in every cell, in order. */
  expression(gene: string): Promise<number[]>;
}

// The worker reads Files with FileReader (and Papa Parse reads CSVs as
// text with it), which Node has not
export function installFileReader(): void {
  if (typeof (globalThis as any).FileReader !== "undefined") return;
  (globalThis as any).FileReader = class {
    result: unknown = null;
    onload: ((e: unknown) => void) | null = null;
    onerror: ((e: unknown) => void) | null = null;
    private read(p: Promise<unknown>) {
      p.then(
        (r) => {
          this.result = r;
          this.onload?.({ target: this });
        },
        (e) => this.onerror?.({ target: { error: e } }),
      );
    }
    readAsArrayBuffer(blob: Blob) {
      this.read(blob.arrayBuffer());
    }
    readAsText(blob: Blob) {
      this.read(blob.text());
    }
  };
}

/** A file as a folder pick hands it over: its path below the folder kept. */
function fileOf(abs: string, relative = ""): File {
  const f = new File([fs.readFileSync(abs)], path.basename(abs));

  Object.defineProperty(f, "webkitRelativePath", { value: relative });

  return f;
}

function walk(dir: string, visit: (abs: string, rel: string) => void): void {
  const go = (d: string) => {
    for (const e of fs.readdirSync(d, { withFileTypes: true })) {
      const abs = path.join(d, e.name);

      if (e.isDirectory()) go(abs);
      else visit(abs, path.relative(dir, abs).split(path.sep).join("/"));
    }
  };

  go(dir);
}

async function describe(ds: any): Promise<Seen> {
  const { selectBestClusterColumn } = await import("@/lib/utils/dataset-utils");
  // The opening column is chosen on the dataset as it is first opened
  // (lib/viewer/open-dataset.ts), before any column is loaded
  const defaultColumn = selectBestClusterColumn(ds) ?? null;
  // A store loads its columns on demand; load them all to see them
  const adapter = ds.adapter;
  const clusters =
    typeof adapter?.getClusterColumnInfo === "function"
      ? await adapter.loadClusters(adapter.getClusterColumnInfo().names)
      : (ds.clusters ?? []);
  const { coordinates: c, dimensions: d } = ds.spatial;
  const cells: number = ds.getPointCount();
  const positions = Array.from({ length: cells }, (_, i) =>
    Array.isArray(c) ? c[i].slice(0, 2) : [c[i * d], c[i * d + 1]],
  );

  return {
    cells,
    genes: [...ds.genes],
    positions,
    embeddings: [
      ...(ds.allEmbeddingNames?.length
        ? ds.allEmbeddingNames
        : Object.keys(ds.embeddings ?? {})),
    ].sort(),
    columns: Object.fromEntries(clusters.map((k: any) => [k.column, k.type])),
    defaultColumn,
    colours: Object.fromEntries(
      clusters
        .filter((k: any) => k.type === "categorical" && k.palette)
        .map((k: any) => [k.column, { ...k.palette }]),
    ),
    expression: async (gene) =>
      Array.from((await ds.getGeneExpression(gene)) ?? []),
  };
}

/**
 * A source as the browser preview opens it: an .h5ad, a Xenium export or a
 * MERSCOPE region.
 */
export async function previewOf(source: string): Promise<Seen> {
  installFileReader();
  const { workerApi } = await import(
    "@/lib/workers/standardized-dataset.worker"
  );
  const { StandardizedDataset } = await import("@/lib/StandardizedDataset");
  let serial;

  if (source.endsWith(".h5ad")) {
    serial = await workerApi.parseH5ad(fileOf(source));
  } else {
    const files: File[] = [];
    const base = path.basename(source);
    const xenium = fs.existsSync(path.join(source, "experiment.xenium"));

    walk(source, (abs, rel) => {
      // The cell adapter never reads a Xenium export's images, molecules or
      // outlines; reading gigabytes of them into Files would be pointless
      if (
        xenium &&
        /morphology|transcripts|\.ome\.tif|\.zarr\.zip|boundaries|\.html$/.test(
          rel,
        )
      ) {
        return;
      }
      files.push(fileOf(abs, `${base}/${rel}`));
    });
    serial = xenium
      ? await workerApi.parseXenium(files)
      : await workerApi.parseMerscope(files);
  }

  return describe(StandardizedDataset.fromSerializedData(serial));
}

/** A store as the SpatialData viewer opens one dropped from disk. */
export async function storeOf(storeDir: string): Promise<Seen> {
  const { openSpatialDataStoreFromFiles } = await import(
    "@/lib/spatialdata/store"
  );
  const { loadSpatialDataCells } = await import(
    "@/lib/spatialdata/loadSpatialDataCells"
  );
  const files = new Map<string, File>();

  walk(storeDir, (abs, rel) => files.set(rel, fileOf(abs)));
  const sd = await openSpatialDataStoreFromFiles(files);

  return describe(
    await loadSpatialDataCells(sd, "parity", path.basename(storeDir)),
  );
}

/**
 * Keep Node's loop busy while stores are read: a read in flight holds no
 * handle there, so the loop can go idle mid-read (a browser tab never does).
 * Returns the function that lets it go.
 */
export function keepAwake(): () => void {
  const timer = setInterval(() => {}, 1000);

  return () => clearInterval(timer);
}
