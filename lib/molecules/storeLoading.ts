/**
 * The legacy molecule folder (`manifest.json.gz` and one gzipped float32
 * file per gene) read through an ElementStore — the app's bucket by
 * presigned URL or a public folder by URL, the same code either way. A
 * dataset made this way loads a gene whole with byte progress, or streams
 * it as it downloads, and both can be cancelled.
 */
import type { Backend, ElementStore } from "@/lib/storage/ElementStore";
import type {
  GeneFetchOptions,
  GeneStreamOptions,
  SingleMoleculeDataset,
} from "@/lib/SingleMoleculeDataset";

import { ungzip } from "pako";

import { hint } from "@/lib/storage/ElementStore";
import {
  streamFloat32Triples,
  supportsStreamingDecompression,
} from "@/lib/utils/stream-float32";

/** 12 bytes: three float32 coordinates. */
const BYTES_PER_MOLECULE = 12;

export const UNASSIGNED_SUFFIX = "_uuuuuuuuuu";

/** A gene's file name, as the processor wrote it. */
export function sanitizeGeneName(gene: string): string {
  return gene
    .replace(/[^a-zA-Z0-9]/g, "_")
    .replace(/_+/g, "_")
    .replace(/^_|_$/g, "");
}

/** The key of a gene's file in a public molecule folder. */
export function folderGeneKey(gene: string, unassigned = false): string {
  return `/genes/${sanitizeGeneName(gene)}${unassigned ? UNASSIGNED_SUFFIX : ""}.bin.gz`;
}

/** The key the app backend understands: the gene's own name. */
export function appGeneKey(gene: string, unassigned = false): string {
  return `/genes/${encodeURIComponent(gene)}${unassigned ? "?unassigned" : ""}`;
}

/**
 * A molecule dataset in the app's bucket: its manifest at the presigned URL
 * the metadata route gave, each gene file at the URL its own route signs.
 */
export function appMoleculeBackend(
  datasetId: string,
  manifestUrl: string,
): Backend {
  const urlOf = async (key: string): Promise<string | null> => {
    if (key === "/manifest.json.gz") return manifestUrl;
    const m = /^\/genes\/([^?]+)(\?unassigned)?$/.exec(key);

    if (!m) return null;
    const res = await fetch(
      `/api/single-molecule/${datasetId}/gene/${m[1]}${m[2] ? "?unassigned=true" : ""}`,
    );

    if (!res.ok) return null;

    return ((await res.json()) as { url: string }).url;
  };
  const fetchUrl = async (
    key: string,
    signal?: AbortSignal,
    init?: RequestInit,
  ) => {
    const url = await urlOf(key);

    if (!url) return undefined;
    const res = await fetch(url, { ...init, signal });

    if (res.status === 404 || res.status === 403) return undefined;
    if (!res.ok) throw new Error(`GET ${res.status} for ${key}`);

    return res;
  };

  return {
    get: async (key, signal) =>
      fetchUrl(key, signal).then(
        async (res) => res && new Uint8Array(await res.arrayBuffer()),
      ),
    getRange: async (key, range, signal) => {
      const res = await fetchUrl(key, signal, {
        headers: {
          Range:
            "suffixLength" in range
              ? `bytes=-${range.suffixLength}`
              : `bytes=${range.offset}-${range.offset + range.length - 1}`,
        },
      });

      if (!res) return undefined;
      const bytes = new Uint8Array(await res.arrayBuffer());

      if (res.status === 206) return bytes;

      return "suffixLength" in range
        ? bytes.slice(-range.suffixLength)
        : bytes.slice(range.offset, range.offset + range.length);
    },
    open: (key, signal) => fetchUrl(key, signal),
  };
}

/** The manifest of a molecule folder, or null when the store has none. */
export async function readMoleculeManifest(
  store: ElementStore,
): Promise<Record<string, any> | null> {
  const bytes = await store.get("/manifest.json.gz");

  if (!bytes) return null;

  return JSON.parse(ungzip(bytes, { to: "string" }));
}

const geneSignal = (signal?: AbortSignal) =>
  hint(signal ?? new AbortController().signal, { tag: "genes", priority: 0 });

/** A response body, whole, reporting bytes as they arrive. */
async function readBody(
  res: Response,
  opts?: GeneFetchOptions,
): Promise<Uint8Array> {
  if (!res.body || !opts?.onProgress) {
    return new Uint8Array(await res.arrayBuffer());
  }
  const totalHeader = res.headers.get("content-length");
  const total = totalHeader ? parseInt(totalHeader, 10) : null;
  const reader = res.body.getReader();
  const chunks: Uint8Array[] = [];
  let loaded = 0;

  for (;;) {
    const { done, value } = await reader.read();

    if (done) break;
    chunks.push(value);
    loaded += value.byteLength;
    opts.onProgress(loaded, total);
  }
  const data = new Uint8Array(loaded);
  let offset = 0;

  for (const chunk of chunks) {
    data.set(chunk, offset);
    offset += chunk.byteLength;
  }

  return data;
}

/**
 * Give a dataset its genes from a store: whole loads, streamed loads and
 * counts, assigned and (when the folder has them) unassigned, each cached
 * once read. `denormFactor` is the scale of a folder written with
 * normalized coordinates (1 for raw ones).
 */
export function installStoreLoading(
  dataset: SingleMoleculeDataset,
  store: ElementStore,
  options: {
    geneKey: (gene: string, unassigned: boolean) => string;
    geneIndex: Map<string, Float32Array>;
    uniqueGenes: string[];
    denormFactor: number;
    hasUnassigned: boolean;
  },
): void {
  const { geneKey, geneIndex, uniqueGenes, denormFactor, hasUnassigned } =
    options;
  const unassignedIndex = new Map<string, Float32Array>();
  const scale = (coords: Float32Array): Float32Array => {
    if (denormFactor === 1) return coords;
    const out = new Float32Array(coords.length);

    for (let i = 0; i < coords.length; i++) out[i] = coords[i] * denormFactor;

    return out;
  };
  const decode = (bytes: Uint8Array) => {
    const raw = ungzip(bytes);

    return scale(new Float32Array(raw.buffer, raw.byteOffset, raw.length / 4));
  };

  /** A gene's file whole; undefined when the folder has none. */
  const readWhole = async (
    key: string,
    opts?: GeneFetchOptions,
  ): Promise<Float32Array | undefined> => {
    const res = await store.open(key, { signal: geneSignal(opts?.signal) });

    return res && decode(await readBody(res, opts));
  };

  /** A gene's file as it downloads; undefined when the folder has none. */
  const readStreamed = async (
    key: string,
    opts: GeneStreamOptions,
  ): Promise<Float32Array | undefined> => {
    const res = await store.open(key, { signal: geneSignal(opts.signal) });

    if (!res) return undefined;
    if (!res.body) throw new Error("Streaming is not supported here");

    return streamFloat32Triples(
      res.body.pipeThrough(new DecompressionStream("gzip")),
      opts.onBatch,
      { scale: denormFactor, signal: opts.signal },
    );
  };

  dataset.getCoordinatesByGene = async function (
    gene: string,
    opts?: GeneFetchOptions,
  ): Promise<Float32Array> {
    const cached = geneIndex.get(gene);

    if (cached) return cached;
    if (!uniqueGenes.includes(gene)) {
      throw new Error(
        `Gene '${gene}' not found. Available genes: ${uniqueGenes.length}`,
      );
    }
    const coords = await readWhole(geneKey(gene, false), opts);

    if (!coords) throw new Error(`No molecule file for gene '${gene}'`);
    geneIndex.set(gene, coords);

    return coords;
  } as unknown as typeof dataset.getCoordinatesByGene;

  dataset.streamCoordinatesByGene = async function (
    gene: string,
    opts: GeneStreamOptions,
  ): Promise<Float32Array> {
    const cached = geneIndex.get(gene);

    if (cached) {
      await opts.onBatch({
        coords: cached,
        loadedMolecules: cached.length / 3,
      });

      return cached;
    }
    if (!supportsStreamingDecompression()) {
      const coords = await dataset.getCoordinatesByGene(gene, {
        signal: opts.signal,
      });

      await opts.onBatch({ coords, loadedMolecules: coords.length / 3 });

      return coords;
    }
    const coords = await readStreamed(geneKey(gene, false), opts);

    if (!coords) throw new Error(`No molecule file for gene '${gene}'`);
    geneIndex.set(gene, coords);

    return coords;
  };

  // Counts for folders whose manifest predates per-gene counts: the gzip
  // trailer holds the uncompressed size, four bytes from the end
  const counts = new Map<string, number | null>();

  dataset.resolveMoleculeCount = async function (
    gene: string,
    kind: "assigned" | "unassigned" = "assigned",
  ): Promise<number | null> {
    const known = dataset.moleculeCountFor(gene, kind);

    if (known !== null) return known;
    if (!uniqueGenes.includes(gene)) return null;
    const id = `${kind}:${gene}`;

    if (counts.has(id)) return counts.get(id)!;
    let count: number | null = null;

    try {
      const tail = await store.getRange(
        geneKey(gene, kind === "unassigned"),
        { suffixLength: 4 },
        { signal: geneSignal() },
      );

      if (tail?.byteLength === 4) {
        const size = new DataView(tail.buffer, tail.byteOffset, 4).getUint32(
          0,
          true,
        );

        if (size > 0 && size % BYTES_PER_MOLECULE === 0) {
          count = size / BYTES_PER_MOLECULE;
        }
      }
    } catch {
      count = null;
    }
    counts.set(id, count);

    return count;
  };

  if (!hasUnassigned) return;

  const empty = new Float32Array(0);
  const isAbort = (e: unknown) => (e as Error)?.name === "AbortError";

  (dataset as any).getUnassignedCoordinatesByGene = async function (
    gene: string,
    opts?: GeneFetchOptions,
  ): Promise<Float32Array> {
    const cached = unassignedIndex.get(gene);

    if (cached) return cached;
    if (!uniqueGenes.includes(gene)) return empty;
    try {
      // A gene with no unassigned molecules has no file: not an error
      const coords = (await readWhole(geneKey(gene, true), opts)) ?? empty;

      unassignedIndex.set(gene, coords);

      return coords;
    } catch (error) {
      if (isAbort(error)) throw error;
      console.warn(`Failed to load unassigned '${gene}':`, error);
      unassignedIndex.set(gene, empty);

      return empty;
    }
  };

  (dataset as any).streamUnassignedCoordinatesByGene = async function (
    gene: string,
    opts: GeneStreamOptions,
  ): Promise<Float32Array> {
    const cached = unassignedIndex.get(gene);

    if (cached) {
      await opts.onBatch({
        coords: cached,
        loadedMolecules: cached.length / 3,
      });

      return cached;
    }
    if (!supportsStreamingDecompression() || !uniqueGenes.includes(gene)) {
      const coords = await dataset.getUnassignedCoordinatesByGene(gene, {
        signal: opts.signal,
      });

      await opts.onBatch({ coords, loadedMolecules: coords.length / 3 });

      return coords;
    }
    try {
      const coords = (await readStreamed(geneKey(gene, true), opts)) ?? empty;

      unassignedIndex.set(gene, coords);
      if (coords === empty) await opts.onBatch({ coords, loadedMolecules: 0 });

      return coords;
    } catch (error) {
      if (isAbort(error)) throw error;
      unassignedIndex.set(gene, empty);
      await opts.onBatch({ coords: empty, loadedMolecules: 0 });

      return empty;
    }
  };
}
