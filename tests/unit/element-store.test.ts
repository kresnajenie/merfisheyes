/**
 * The one store every reader reads through, against a backend that counts
 * what it is asked: reads are cached and shared, cancelled all the way to
 * the backend, merged when ranged and neighbouring, ordered by priority
 * and bounded in flight, and answered from consolidated metadata.
 */
import type { Backend } from "@/lib/storage/ElementStore";

import { describe, expect, it, vi } from "vitest";

import { ElementStore, hint, hinted } from "@/lib/storage/ElementStore";

const abortErr = () => new DOMException("aborted", "AbortError");

/** A backend over an object map that records every call and can be slowed. */
function fake(objects: Record<string, Uint8Array> = {}) {
  const calls: { kind: "get" | "range"; key: string; range?: unknown }[] = [];
  let gate: Promise<void> | null = null;
  let open = () => {};
  const waitForGate = () => gate ?? Promise.resolve();
  const abortable = async (signal?: AbortSignal) => {
    if (signal?.aborted) throw abortErr();
    const race = new Promise<never>((_, reject) =>
      signal?.addEventListener("abort", () => reject(abortErr()), {
        once: true,
      }),
    );

    await Promise.race([waitForGate(), race]);
    if (signal?.aborted) throw abortErr();
  };
  const backend: Backend = {
    async get(key, signal) {
      calls.push({ kind: "get", key });
      await abortable(signal);

      return objects[key];
    },
    async getRange(key, range, signal) {
      calls.push({ kind: "range", key, range });
      await abortable(signal);
      const o = objects[key];

      if (!o) return undefined;

      return "suffixLength" in range
        ? o.slice(-range.suffixLength)
        : o.slice(range.offset, range.offset + range.length);
    },
  };

  return {
    backend,
    calls,
    /** Hold every request until release() */
    hold() {
      gate = new Promise((r) => (open = r));
    },
    release() {
      open();
      gate = null;
    },
  };
}

const bytes = (n: number) => Uint8Array.from({ length: n }, (_, i) => i % 251);
const tick = () => new Promise((r) => setTimeout(r, 5));

describe("ElementStore", () => {
  it("reads an object once and answers repeats from the cache", async () => {
    const f = fake({ "/a": bytes(10) });
    const store = new ElementStore(f.backend);

    expect(await store.get("/a")).toEqual(bytes(10));
    expect(await store.get("a")).toEqual(bytes(10));
    expect(await store.get("/missing")).toBeUndefined();
    expect(await store.get("/missing")).toBeUndefined();
    expect(f.calls.map((c) => c.key)).toEqual(["/a", "/missing"]);
    const s = store.stats();

    expect(s.requests).toBe(2);
    expect(s.bytes).toBe(10);
    expect(s.cacheHits).toBe(2);
    expect(s.cacheBytes).toBe(10);
  });

  it("shares a read in flight and keeps the cache within its bytes", async () => {
    const f = fake({ "/a": bytes(10), "/b": bytes(10), "/c": bytes(10) });
    const store = new ElementStore(f.backend, { cacheBytes: 25 });

    f.hold();
    const [x, y] = [store.get("/a"), store.get("/a")];

    f.release();
    expect(await x).toEqual(await y);
    expect(f.calls.length).toBe(1);
    await store.get("/b");
    await store.get("/c");
    // a was the least recently used
    expect(store.stats().cacheBytes).toBe(20);
    await store.get("/a");
    expect(f.calls.length).toBe(4);
  });

  it("aborts the backend read when its last waiter gives up, not before", async () => {
    const f = fake({ "/a": bytes(10) });
    const store = new ElementStore(f.backend);
    const one = new AbortController();
    const two = new AbortController();

    f.hold();
    const first = store.get("/a", { signal: one.signal });
    const second = store.get("/a", { signal: two.signal });

    one.abort();
    await expect(first).rejects.toMatchObject({ name: "AbortError" });
    f.release();
    expect(await second).toEqual(bytes(10));
    expect(store.stats().aborted).toBe(0);

    f.hold();
    const third = store.get("/z", { signal: two.signal });

    two.abort();
    await expect(third).rejects.toMatchObject({ name: "AbortError" });
    f.release();
    await tick();
    expect(store.stats().aborted).toBe(1);
    expect(store.stats().inFlight).toBe(0);
  });

  it("merges neighbouring ranges of one object into one request", async () => {
    const f = fake({ "/shard": bytes(1000) });
    const store = new ElementStore(f.backend, { mergeGap: 10 });
    const parts = await Promise.all([
      store.getRange("/shard", { offset: 0, length: 100 }),
      store.getRange("/shard", { offset: 105, length: 100 }),
      store.getRange("/shard", { offset: 500, length: 50 }),
    ]);

    expect(parts[0]).toEqual(bytes(1000).slice(0, 100));
    expect(parts[1]).toEqual(bytes(1000).slice(105, 205));
    expect(parts[2]).toEqual(bytes(1000).slice(500, 550));
    expect(f.calls.map((c) => c.range)).toEqual([
      { offset: 0, length: 205 },
      { offset: 500, length: 50 },
    ]);
    expect(store.stats().merged).toBe(1);
    // Each member is cached on its own
    await store.getRange("/shard", { offset: 105, length: 100 });
    expect(f.calls.length).toBe(2);
    // A suffix range (a shard index) is its own read, cached
    await store.getRange("/shard", { suffixLength: 8 });
    await store.getRange("/shard", { suffixLength: 8 });
    expect(f.calls.length).toBe(3);
  });

  it("runs the highest priority first and bounds what is in flight", async () => {
    const f = fake({ "/a": bytes(1), "/b": bytes(1), "/c": bytes(1), "/d": bytes(1) });
    const store = new ElementStore(f.backend, { concurrency: 1 });

    f.hold();
    const a = store.get("/a", { signal: hinted({ priority: 0 }) });
    const b = store.get("/b", { signal: hinted({ priority: 1 }) });
    const c = store.get("/c", { signal: hinted({ priority: 5 }) });
    const d = store.get("/d"); // no signal: before everything queued

    await Promise.resolve();
    expect(store.stats().inFlight).toBe(1);
    expect(store.stats().queued).toBe(3);
    f.release();
    await Promise.all([a, b, c, d]);
    expect(f.calls.map((x) => x.key)).toEqual(["/a", "/d", "/c", "/b"]);
  });

  it("drops a queued read that is aborted before it starts", async () => {
    const f = fake({ "/a": bytes(1), "/b": bytes(1) });
    const store = new ElementStore(f.backend, { concurrency: 1 });
    const ctl = new AbortController();

    f.hold();
    const a = store.get("/a");
    const b = store.get("/b", { signal: ctl.signal });

    ctl.abort();
    await expect(b).rejects.toMatchObject({ name: "AbortError" });
    f.release();
    await a;
    expect(f.calls.map((x) => x.key)).toEqual(["/a"]);
    expect(store.stats().aborted).toBe(1);
  });

  it("counts reads per tag", async () => {
    const f = fake({ "/a": bytes(10), "/b": bytes(5) });
    const store = new ElementStore(f.backend);
    const points = hinted({ tag: "points" });

    await store.get("/a", { signal: points });
    await store.get("/b", { signal: hint(new AbortController().signal, { tag: "image" }) });
    expect(store.stats().byTag).toEqual({
      points: { requests: 1, bytes: 10 },
      image: { requests: 1, bytes: 5 },
    });
  });

  it("answers every node's zarr.json from the consolidated metadata", async () => {
    const f = fake({});
    const store = new ElementStore(f.backend);
    const root = {
      zarr_format: 3,
      node_type: "group",
      consolidated_metadata: {
        metadata: {
          tables: { zarr_format: 3, node_type: "group" },
          "tables/table/X": { zarr_format: 3, node_type: "array", shape: [2] },
        },
      },
    };

    expect(store.seedConsolidated(root)).toEqual(["tables", "tables/table/X"]);
    const meta = JSON.parse(
      new TextDecoder().decode(await store.get("/tables/table/X/zarr.json")),
    );

    expect(meta.shape).toEqual([2]);
    expect(f.calls).toEqual([]);
    expect(store.seedConsolidated({ zarr_format: 3 })).toBeNull();
  });

  it("refuses to stream without a streaming backend and streams with one", async () => {
    const f = fake({});
    const store = new ElementStore(f.backend);

    await expect(store.open("/a")).rejects.toThrow(/cannot stream/);
    const streaming = new ElementStore({
      ...f.backend,
      open: async () =>
        new Response(bytes(4), { headers: { "content-length": "4" } }),
    });
    const res = await streaming.open("/a", { signal: hinted({ tag: "genes" }) });

    expect(new Uint8Array(await res!.arrayBuffer())).toEqual(bytes(4));
    expect(streaming.stats().byTag.genes).toEqual({ requests: 1, bytes: 4 });
    expect(streaming.stats().inFlight).toBe(0);
  });

  it("fails every waiter when the backend fails, and forgets the read", async () => {
    const backend: Backend = {
      get: vi.fn(async () => {
        throw new Error("boom");
      }),
      getRange: async () => undefined,
    };
    const store = new ElementStore(backend);

    await expect(
      Promise.all([store.get("/a"), store.get("/a")]),
    ).rejects.toThrow("boom");
    await expect(store.get("/a")).rejects.toThrow("boom");
    expect(backend.get).toHaveBeenCalledTimes(2);
  });
});
