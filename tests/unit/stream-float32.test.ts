import { describe, it, expect, vi } from "vitest";
import { gzipSync } from "zlib";

import {
  concatFloat32,
  streamFloat32Triples,
} from "@/lib/utils/stream-float32";

/** A ReadableStream that emits `bytes` in fixed-size pieces. */
function chunkedStream(bytes: Uint8Array, chunkSize: number): ReadableStream<Uint8Array> {
  let offset = 0;

  return new ReadableStream<Uint8Array>({
    pull(controller) {
      if (offset >= bytes.length) {
        controller.close();

        return;
      }
      controller.enqueue(bytes.slice(offset, offset + chunkSize));
      offset += chunkSize;
    },
  });
}

function coordsToBytes(values: number[]): Uint8Array {
  const f32 = Float32Array.from(values);

  return new Uint8Array(f32.buffer);
}

describe("streamFloat32Triples", () => {
  const molecules = 1000;
  const values = Array.from({ length: molecules * 3 }, (_, i) => i * 0.5);

  it("reassembles every molecule when chunks split mid-triple", async () => {
    // 7 is coprime with 12, so nearly every chunk boundary lands inside a
    // molecule — the case a naive reader corrupts.
    const batches: number[][] = [];
    const all = await streamFloat32Triples(
      chunkedStream(coordsToBytes(values), 7),
      ({ coords }) => {
        batches.push(Array.from(coords));
      },
      { batchMolecules: 100, batchIntervalMs: 1_000_000 },
    );

    expect(all.length).toBe(values.length);
    expect(Array.from(all)).toEqual(values);
    // Every batch is whole molecules
    for (const b of batches) expect(b.length % 3).toBe(0);
    expect(batches.flat()).toEqual(values);
  });

  it("emits batches of the requested size and reports cumulative progress", async () => {
    const progress: number[] = [];

    await streamFloat32Triples(
      chunkedStream(coordsToBytes(values), 1200),
      ({ loadedMolecules }) => {
        progress.push(loadedMolecules);
      },
      { batchMolecules: 250, batchIntervalMs: 1_000_000 },
    );

    expect(progress[progress.length - 1]).toBe(molecules);
    // Monotonic, and each step is at most one batch's worth of molecules.
    for (let i = 1; i < progress.length; i++) {
      expect(progress[i]).toBeGreaterThan(progress[i - 1]);
    }
  });

  it("applies the legacy scale factor", async () => {
    const all = await streamFloat32Triples(
      chunkedStream(coordsToBytes([1, 2, 3]), 12),
      () => {},
      { scale: 10 },
    );

    expect(Array.from(all)).toEqual([10, 20, 30]);
  });

  it("drops a trailing partial molecule instead of throwing", async () => {
    const warn = vi.spyOn(console, "warn").mockImplementation(() => {});
    const bytes = coordsToBytes([1, 2, 3, 4, 5, 6]);
    const truncated = bytes.slice(0, bytes.length - 5); // 1.5 molecules

    const all = await streamFloat32Triples(chunkedStream(truncated, 5), () => {}, {});

    expect(Array.from(all)).toEqual([1, 2, 3]);
    expect(warn).toHaveBeenCalled();
    warn.mockRestore();
  });

  it("stops when the signal is aborted", async () => {
    const controller = new AbortController();

    await expect(
      streamFloat32Triples(
        chunkedStream(coordsToBytes(values), 12),
        ({ loadedMolecules }) => {
          if (loadedMolecules >= 10) controller.abort();
        },
        { batchMolecules: 5, batchIntervalMs: 1_000_000, signal: controller.signal },
      ),
    ).rejects.toThrow(/abort/i);
  });

  it("reads real gzip output through DecompressionStream", async () => {
    const gz = gzipSync(Buffer.from(coordsToBytes(values)));
    const stream = chunkedStream(new Uint8Array(gz), 64).pipeThrough(
      new DecompressionStream("gzip"),
    );
    const all = await streamFloat32Triples(stream, () => {}, {
      batchMolecules: 128,
      batchIntervalMs: 1_000_000,
    });

    expect(Array.from(all)).toEqual(values);
  });
});

describe("concatFloat32", () => {
  it("joins batches in order", () => {
    const out = concatFloat32([
      Float32Array.from([1, 2]),
      Float32Array.from([3]),
      Float32Array.from([4, 5, 6]),
    ]);

    expect(Array.from(out)).toEqual([1, 2, 3, 4, 5, 6]);
  });

  it("returns the single batch untouched", () => {
    const only = Float32Array.from([1, 2, 3]);

    expect(concatFloat32([only])).toBe(only);
  });
});
