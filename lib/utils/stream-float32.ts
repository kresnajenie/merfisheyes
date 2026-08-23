/**
 * Read a gzipped stream of float32 triples (`[x,y,z, x,y,z, …]`) in batches,
 * so a huge gene file can be rendered while it downloads instead of after.
 *
 * The single-molecule gene files are exactly this: one gzip member wrapping a
 * flat Float32Array. Nothing about the format has to change to consume it
 * progressively — a molecule is complete as soon as its 12th byte arrives.
 */

/** 12 bytes: three float32 coordinates. */
const BYTES_PER_MOLECULE = 12;

export interface StreamFloat32Options {
  /** Emit a batch once this many molecules are pending (default 250 000). */
  batchMolecules?: number;
  /** …or once this long has passed since the last batch (ms, default 120). */
  batchIntervalMs?: number;
  /** Multiply every coordinate (legacy normalized datasets). */
  scale?: number;
  signal?: AbortSignal;
}

export interface StreamFloat32Batch {
  /** Coordinates of this batch only: `[x,y,z, …]`, length divisible by 3. */
  coords: Float32Array;
  /** Molecules emitted so far, including this batch. */
  loadedMolecules: number;
}

/**
 * Consume a byte stream as float32 triples, calling `onBatch` as data arrives.
 * Returns every coordinate read, concatenated (so callers can cache the gene).
 *
 * Bytes that don't complete a triple are carried into the next chunk; a
 * trailing partial molecule at end-of-stream is dropped (truncated file) and
 * reported via the console rather than throwing, so a partly-rendered gene
 * stays usable.
 */
export async function streamFloat32Triples(
  stream: ReadableStream<Uint8Array>,
  onBatch: (batch: StreamFloat32Batch) => void | Promise<void>,
  options: StreamFloat32Options = {},
): Promise<Float32Array> {
  const batchMolecules = options.batchMolecules ?? 250_000;
  const batchIntervalMs = options.batchIntervalMs ?? 120;
  const scale = options.scale ?? 1;
  const batchBytesTarget = batchMolecules * BYTES_PER_MOLECULE;

  const reader = stream.getReader();
  const collected: Float32Array[] = [];

  // Bytes that arrived but don't yet complete a molecule.
  let carry: Uint8Array = new Uint8Array(0);
  // Complete-molecule byte runs waiting to be emitted as one batch.
  let pending: Uint8Array[] = [];
  let pendingBytes = 0;
  let loadedMolecules = 0;
  let lastFlush = Date.now();

  const flush = async () => {
    if (pendingBytes === 0) return;

    // One fresh allocation per batch: a Float32Array view needs its own
    // 4-byte-aligned buffer, which subarrays of network chunks are not.
    const bytes = new Uint8Array(pendingBytes);
    let offset = 0;

    for (const part of pending) {
      bytes.set(part, offset);
      offset += part.length;
    }
    pending = [];
    pendingBytes = 0;

    const coords = new Float32Array(bytes.buffer, 0, bytes.length / 4);

    if (scale !== 1) {
      for (let i = 0; i < coords.length; i++) coords[i] *= scale;
    }

    collected.push(coords);
    loadedMolecules += coords.length / 3;
    lastFlush = Date.now();
    await onBatch({ coords, loadedMolecules });
  };

  try {
    for (;;) {
      if (options.signal?.aborted) {
        throw new DOMException("Aborted", "AbortError");
      }

      const { done, value } = await reader.read();

      if (done) break;
      if (!value || value.length === 0) continue;

      let chunk: Uint8Array;

      if (carry.length > 0) {
        chunk = new Uint8Array(carry.length + value.length);
        chunk.set(carry, 0);
        chunk.set(value, carry.length);
      } else {
        chunk = value;
      }

      const usable = Math.floor(chunk.length / BYTES_PER_MOLECULE) * BYTES_PER_MOLECULE;

      if (usable > 0) {
        pending.push(chunk.subarray(0, usable));
        pendingBytes += usable;
      }
      // slice (not subarray): the leftover must survive the chunk it came from.
      carry = chunk.slice(usable);

      if (
        pendingBytes >= batchBytesTarget ||
        (pendingBytes > 0 && Date.now() - lastFlush >= batchIntervalMs)
      ) {
        await flush();
      }
    }

    await flush();

    if (carry.length > 0) {
      console.warn(
        `[streamFloat32Triples] ${carry.length} trailing byte(s) did not complete a molecule — file may be truncated`,
      );
    }
  } finally {
    reader.releaseLock();
  }

  return concatFloat32(collected);
}

/** Join batch arrays into the single array a gene cache expects. */
export function concatFloat32(parts: Float32Array[]): Float32Array {
  if (parts.length === 1) return parts[0];
  let total = 0;

  for (const p of parts) total += p.length;
  const out = new Float32Array(total);
  let offset = 0;

  for (const p of parts) {
    out.set(p, offset);
    offset += p.length;
  }

  return out;
}

/**
 * Fetch a gzipped gene file and stream it as float32 triples.
 *
 * `DecompressionStream` keeps the gunzip incremental — the whole point, since
 * the alternative (`arrayBuffer()` then `ungzip`) is what makes a 50M-molecule
 * gene take minutes before anything appears.
 */
export async function streamGzippedFloat32(
  url: string,
  onBatch: (batch: StreamFloat32Batch) => void | Promise<void>,
  options: StreamFloat32Options = {},
): Promise<Float32Array> {
  const response = await fetch(url, { signal: options.signal });

  if (!response.ok) {
    throw new Error(`Failed to download ${url}: ${response.status}`);
  }
  if (!response.body) {
    throw new Error("Streaming is not supported for this response");
  }

  const stream = response.body.pipeThrough(new DecompressionStream("gzip"));

  return streamFloat32Triples(stream, onBatch, options);
}

/** True when this runtime can stream + gunzip incrementally. */
export function supportsStreamingDecompression(): boolean {
  return (
    typeof DecompressionStream !== "undefined" &&
    typeof ReadableStream !== "undefined"
  );
}

/**
 * Molecules in a gzipped gene file, without downloading it.
 *
 * A gzip member ends with ISIZE: the uncompressed size, little-endian uint32.
 * One 4-byte range request therefore gives the exact molecule count
 * (`ISIZE / 12`) — which is what sizes the streaming buffer for datasets whose
 * manifest predates per-gene `molecule_counts`.
 *
 * Returns null when the count can't be trusted: no range support, a size that
 * isn't a whole number of molecules, or a file at/over 4 GB (ISIZE is stored
 * mod 2^32 and would wrap — detected by the uncompressed size coming out
 * smaller than the compressed one).
 */
export async function probeGzipMoleculeCount(
  url: string,
  signal?: AbortSignal,
): Promise<number | null> {
  try {
    const response = await fetch(url, {
      headers: { Range: "bytes=-4" },
      signal,
    });

    if (response.status !== 206) return null;
    const buffer = await response.arrayBuffer();

    if (buffer.byteLength !== 4) return null;

    const uncompressedBytes = new DataView(buffer).getUint32(0, true);

    if (uncompressedBytes === 0 || uncompressedBytes % BYTES_PER_MOLECULE !== 0) {
      return null;
    }

    // content-range: "bytes <start>-<end>/<total>" — total is the compressed size.
    const total = Number(
      response.headers.get("content-range")?.split("/")[1] ?? NaN,
    );

    if (Number.isFinite(total) && total > uncompressedBytes) return null;

    return uncompressedBytes / BYTES_PER_MOLECULE;
  } catch {
    return null;
  }
}
