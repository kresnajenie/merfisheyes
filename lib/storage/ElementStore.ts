/**
 * The one store every reader reads through (docs/specs/LOADING.md): a
 * SpatialData table, a points element, polygons, an image, a legacy
 * molecule folder. Whatever is behind it — a public URL, a dataset in the
 * app's bucket, files dropped in the browser — a read is the same call,
 * and the store is where loading is governed:
 *
 * - one byte cache (least recently used out first, bounded), with reads in
 *   flight shared, so a chunk that two tiles need is fetched once;
 * - every read takes an `AbortSignal` that reaches `fetch`; a shared read is
 *   aborted when its last waiter has gone;
 * - ranged reads of one object issued close together are merged into one
 *   request (sharded arrays: many small chunks in one object);
 * - a priority queue bounds the requests in flight per store, so an image
 *   tile cannot starve the molecules and a coarse level goes first;
 * - the consolidated metadata of a Zarr v3 store is read once and answers
 *   every `zarr.json` from memory;
 * - what was read is counted, per reader when the reader says who it is.
 *
 * zarrita forwards only the signal of a read to the store, so a reader's
 * hints (priority, who it is) ride on the signal: `hint(signal, {...})`.
 */
import type { AbsolutePath, AsyncReadable, RangeQuery } from "@zarrita/storage";

export interface ReadOptions {
  signal?: AbortSignal;
}

/** What a reader can say about its reads, attached to their signal. */
export interface ReadHints {
  /** Higher runs first; reads without a signal run before everything. */
  priority?: number;
  /** Who is reading (an element kind), for the per-reader counts. */
  tag?: string;
}

export interface StoreStats {
  /** Requests sent to the backend. */
  requests: number;
  /** Bytes the backend answered with. */
  bytes: number;
  inFlight: number;
  queued: number;
  /** Reads answered from the cache or from a read already in flight. */
  cacheHits: number;
  /** Ranged reads folded into another request. */
  merged: number;
  /** Requests cancelled before or while running. */
  aborted: number;
  /** Bytes held in the cache. */
  cacheBytes: number;
  byTag: Record<string, { requests: number; bytes: number }>;
}

/** What a store reads from. */
export interface Backend {
  /** A whole object; undefined when there is none. */
  get(key: string, signal?: AbortSignal): Promise<Uint8Array | undefined>;
  /** Part of an object; undefined when there is none. */
  getRange(
    key: string,
    range: RangeQuery,
    signal?: AbortSignal,
  ): Promise<Uint8Array | undefined>;
  /** An object as a stream, for files consumed as they arrive. */
  open?(key: string, signal?: AbortSignal): Promise<Response | undefined>;
}

export interface ElementStoreOptions {
  /** Bytes the cache may hold (default 256 MB). */
  cacheBytes?: number;
  /** Requests in flight at once (default 32). */
  concurrency?: number;
  /** Ranged reads this close are merged (default 64 KB). */
  mergeGap?: number;
  /** A merged request reads at most this much (default 8 MB). */
  mergeMax?: number;
}

const hints = new WeakMap<AbortSignal, ReadHints>();

/** Attach hints to the signal a reader passes with its reads. */
export function hint(signal: AbortSignal, h: ReadHints): AbortSignal {
  hints.set(signal, { ...hints.get(signal), ...h });

  return signal;
}

/** A signal carrying hints, for readers that have nothing to cancel. */
export function hinted(h: ReadHints): AbortSignal {
  return hint(new AbortController().signal, h);
}

const abortError = () =>
  typeof DOMException !== "undefined"
    ? new DOMException("The read was aborted", "AbortError")
    : Object.assign(new Error("The read was aborted"), { name: "AbortError" });

/** Reject when the signal aborts, whatever the promise does afterwards. */
function untilAborted<T>(p: Promise<T>, signal?: AbortSignal): Promise<T> {
  if (!signal) return p;
  if (signal.aborted) return Promise.reject(abortError());

  return new Promise<T>((resolve, reject) => {
    const onAbort = () => reject(abortError());

    signal.addEventListener("abort", onAbort, { once: true });
    p.then(resolve, reject).finally(() =>
      signal.removeEventListener("abort", onAbort),
    );
  });
}

/** One request to the backend: shared by every read that wants it. */
interface Flight {
  promise: Promise<Uint8Array | undefined>;
  controller: AbortController;
  waiters: number;
  started: boolean;
  priority: number;
}

interface PendingRange {
  offset: number;
  length: number;
  signal?: AbortSignal;
  resolve: (b: Uint8Array | undefined) => void;
  reject: (e: unknown) => void;
}

export class ElementStore implements AsyncReadable {
  private cache = new Map<string, Uint8Array | undefined>();
  private cacheHeld = 0;
  private flights = new Map<string, Flight>();
  private queue: { priority: number; seq: number; run: () => void }[] = [];
  private seq = 0;
  private running = 0;
  private pendingRanges = new Map<string, PendingRange[]>();
  private flushPending = false;
  private _stats: StoreStats = {
    requests: 0,
    bytes: 0,
    inFlight: 0,
    queued: 0,
    cacheHits: 0,
    merged: 0,
    aborted: 0,
    cacheBytes: 0,
    byTag: {},
  };
  readonly options: Required<ElementStoreOptions>;

  constructor(
    private backend: Backend,
    options: ElementStoreOptions = {},
  ) {
    this.options = {
      cacheBytes: options.cacheBytes ?? 256 * 1024 * 1024,
      concurrency: options.concurrency ?? 32,
      mergeGap: options.mergeGap ?? 64 * 1024,
      mergeMax: options.mergeMax ?? 8 * 1024 * 1024,
    };
  }

  stats(): StoreStats {
    return {
      ...this._stats,
      queued: this.queue.length,
      cacheBytes: this.cacheHeld,
      byTag: Object.fromEntries(
        Object.entries(this._stats.byTag).map(([k, v]) => [k, { ...v }]),
      ),
    };
  }

  /** Change how much the cache may hold; what no longer fits goes. */
  setCacheBytes(bytes: number): void {
    this.options.cacheBytes = bytes;
    this.trim();
  }

  /** Put an object read elsewhere (the root `zarr.json`) into the cache. */
  seed(key: string, bytes: Uint8Array): void {
    this.remember(norm(key), bytes);
  }

  /**
   * Answer every node's `zarr.json` from the root's consolidated metadata,
   * so opening an array or group costs no request. Returns the node paths
   * listed, or null when the root has no consolidated metadata.
   */
  seedConsolidated(root: Record<string, any>): string[] | null {
    const nodes = root.consolidated_metadata?.metadata;

    if (!nodes || typeof nodes !== "object") return null;
    const encoder = new TextEncoder();

    for (const [path, meta] of Object.entries(nodes)) {
      this.remember(
        `/${path.replace(/^\/+/, "")}/zarr.json`,
        encoder.encode(JSON.stringify(meta)),
      );
    }

    return Object.keys(nodes);
  }

  get(key: AbsolutePath | string, opts?: ReadOptions) {
    const k = norm(key);

    if (this.cache.has(k)) return this.hit(k, opts?.signal);

    return this.share(
      k,
      opts?.signal,
      (signal) => this.backend.get(k, signal),
      true,
    );
  }

  getRange(key: AbsolutePath | string, range: RangeQuery, opts?: ReadOptions) {
    const k = norm(key);

    if ("suffixLength" in range) {
      const ck = `${k}#-${range.suffixLength}`;

      if (this.cache.has(ck)) return this.hit(ck, opts?.signal);

      return this.share(
        ck,
        opts?.signal,
        (signal) => this.backend.getRange(k, range, signal),
        true,
      );
    }
    const ck = `${k}#${range.offset}+${range.length}`;

    if (this.cache.has(ck)) return this.hit(ck, opts?.signal);
    const flight = this.flights.get(ck);

    if (flight) return this.join(flight, opts?.signal);

    // Collected for a moment, so neighbouring ranges go as one request
    return new Promise<Uint8Array | undefined>((resolve, reject) => {
      const list = this.pendingRanges.get(k) ?? [];

      list.push({ ...range, signal: opts?.signal, resolve, reject });
      this.pendingRanges.set(k, list);
      if (!this.flushPending) {
        this.flushPending = true;
        nextTask(() => this.flushRanges());
      }
    });
  }

  /** The object as a stream: not cached, not merged, but queued and counted. */
  async open(
    key: AbsolutePath | string,
    opts?: ReadOptions,
  ): Promise<Response | undefined> {
    if (!this.backend.open) throw new Error("This store cannot stream");
    const k = norm(key);
    const signal = opts?.signal;
    const h = signal ? hints.get(signal) : undefined;

    await this.slot(h?.priority ?? Infinity, signal);
    this.count(h?.tag, 0);
    this._stats.inFlight++;
    try {
      const res = await this.backend.open(k, signal);
      const length = Number(res?.headers.get("content-length") ?? 0);

      if (length) this.count(h?.tag, length, false);

      return res;
    } catch (e) {
      if ((e as Error)?.name === "AbortError") this._stats.aborted++;
      throw e;
    } finally {
      this._stats.inFlight--;
      this.release();
    }
  }

  // ── cache ────────────────────────────────────────────────────────────────

  private hit(ck: string, signal?: AbortSignal) {
    this._stats.cacheHits++;
    const bytes = this.cache.get(ck);

    // Most recently used last
    this.cache.delete(ck);
    this.cache.set(ck, bytes);

    return untilAborted(Promise.resolve(bytes), signal);
  }

  private remember(ck: string, bytes: Uint8Array | undefined) {
    const old = this.cache.get(ck);

    if (old) this.cacheHeld -= old.byteLength;
    this.cache.delete(ck);
    this.cache.set(ck, bytes);
    this.cacheHeld += bytes?.byteLength ?? 0;
    this.trim();
  }

  private trim() {
    for (const [k, v] of this.cache) {
      if (this.cacheHeld <= this.options.cacheBytes) break;
      this.cache.delete(k);
      this.cacheHeld -= v?.byteLength ?? 0;
    }
  }

  // ── shared requests ──────────────────────────────────────────────────────

  /** One request for everyone who asks for `ck` while it is in flight. */
  private share(
    ck: string,
    signal: AbortSignal | undefined,
    load: (signal: AbortSignal) => Promise<Uint8Array | undefined>,
    cache: boolean,
  ): Promise<Uint8Array | undefined> {
    let flight = this.flights.get(ck);

    if (!flight) {
      const controller = new AbortController();
      const h = signal ? hints.get(signal) : undefined;
      const f: Flight = {
        controller,
        waiters: 0,
        started: false,
        priority: h?.priority ?? Infinity,
        promise: this.slot(h?.priority ?? Infinity, controller.signal)
          .then(() => {
            f.started = true;
            this._stats.inFlight++;
            this.count(h?.tag, 0);

            return load(controller.signal);
          })
          .then(
            (bytes) => {
              this.count(h?.tag, bytes?.byteLength ?? 0, false);
              if (cache) this.remember(ck, bytes);

              return bytes;
            },
            (e) => {
              // Cancelled while queued was counted by the queue
              if (f.started && (e as Error)?.name === "AbortError") {
                this._stats.aborted++;
              }
              throw e;
            },
          )
          .finally(() => {
            if (f.started) {
              this._stats.inFlight--;
              this.release();
            }
            this.flights.delete(ck);
          }),
      };

      flight = f;
      this.flights.set(ck, f);
      // Nobody listens: a failure is reported to each waiter
      f.promise.catch(() => {});
    } else {
      this._stats.cacheHits++;
    }

    return this.join(flight, signal);
  }

  /** Wait on a flight; when its last waiter aborts, so does the request. */
  private join(flight: Flight, signal?: AbortSignal) {
    flight.waiters++;
    let left = false;
    const leave = () => {
      if (left) return;
      left = true;
      if (--flight.waiters === 0) flight.controller.abort();
    };

    if (signal) {
      if (signal.aborted) {
        leave();

        return Promise.reject(abortError());
      }
      signal.addEventListener("abort", leave, { once: true });
    }

    return untilAborted(flight.promise, signal).finally(() => {
      signal?.removeEventListener("abort", leave);
      if (!left) {
        left = true;
        flight.waiters--;
      }
    });
  }

  // ── merged ranges ────────────────────────────────────────────────────────

  private flushRanges() {
    this.flushPending = false;
    const all = this.pendingRanges;

    this.pendingRanges = new Map();
    for (const [k, list] of all) {
      list.sort((a, b) => a.offset - b.offset);
      let group: PendingRange[] = [];
      let end = -Infinity;

      const flush = () => {
        if (group.length) this.requestGroup(k, group);
        group = [];
      };

      for (const r of list) {
        const start = group.length ? group[0].offset : r.offset;
        const fits =
          group.length === 0 ||
          (r.offset <= end + this.options.mergeGap &&
            Math.max(end, r.offset + r.length) - start <=
              this.options.mergeMax);

        if (!fits) flush();
        group.push(r);
        end = Math.max(
          group.length === 1 ? -Infinity : end,
          r.offset + r.length,
        );
      }
      flush();
    }
  }

  /** One request for a group of neighbouring ranges of one object. */
  private requestGroup(k: string, group: PendingRange[]) {
    const live = group.filter((r) => !r.signal?.aborted);

    for (const r of group) if (r.signal?.aborted) r.reject(abortError());
    if (!live.length) return;
    const offset = live[0].offset;
    const end = Math.max(...live.map((r) => r.offset + r.length));
    const length = end - offset;

    this._stats.merged += live.length - 1;
    // The members share one flight through one combined signal, each still
    // able to leave on its own
    const controller = new AbortController();
    let waiting = live.length;
    const hs = live.map((r) => (r.signal ? hints.get(r.signal) : undefined));
    const priority = Math.max(...hs.map((h) => h?.priority ?? Infinity));

    hint(controller.signal, {
      priority,
      tag: hs.find((h) => h?.tag)?.tag,
    });
    const ck = `${k}#${offset}+${length}`;
    const whole = this.share(
      ck,
      controller.signal,
      (signal) => this.backend.getRange(k, { offset, length }, signal),
      live.length === 1,
    );

    for (const r of live) {
      const own = untilAborted(whole, r.signal).then((bytes) => {
        if (!bytes) return undefined;
        const part = bytes.slice(
          r.offset - offset,
          r.offset - offset + r.length,
        );

        if (live.length > 1)
          this.remember(`${k}#${r.offset}+${r.length}`, part);

        return part;
      });

      own.then(r.resolve, r.reject).finally(() => {
        if (--waiting === 0) controller.abort();
      });
    }
  }

  // ── the queue ────────────────────────────────────────────────────────────

  /** Resolves when a request may start; rejects if the signal aborts first. */
  private slot(priority: number, signal?: AbortSignal): Promise<void> {
    if (this.running < this.options.concurrency && !this.queue.length) {
      this.running++;

      return Promise.resolve();
    }

    return new Promise<void>((resolve, reject) => {
      const entry = {
        priority,
        seq: this.seq++,
        run: () => {
          signal?.removeEventListener("abort", cancel);
          this.running++;
          resolve();
        },
      };
      const cancel = () => {
        const i = this.queue.indexOf(entry);

        if (i >= 0) this.queue.splice(i, 1);
        this._stats.aborted++;
        reject(abortError());
      };

      signal?.addEventListener("abort", cancel, { once: true });
      this.queue.push(entry);
      this.queue.sort((a, b) => b.priority - a.priority || a.seq - b.seq);
    });
  }

  private release() {
    this.running--;
    while (this.running < this.options.concurrency && this.queue.length) {
      this.queue.shift()!.run();
    }
  }

  private count(tag: string | undefined, bytes: number, request = true) {
    if (request) this._stats.requests++;
    this._stats.bytes += bytes;
    if (!tag) return;
    const t = (this._stats.byTag[tag] ??= { requests: 0, bytes: 0 });

    if (request) t.requests++;
    t.bytes += bytes;
  }
}

/**
 * Run after the current task: the window in which neighbouring ranged reads
 * of one object are collected. A message channel, not setTimeout, because
 * browsers clamp nested timers to 4 ms and a sharded array is read as a
 * long chain of small ranged reads.
 */
const nextTask: (fn: () => void) => void = (() => {
  if (typeof MessageChannel === "undefined") {
    return (fn) => setTimeout(fn, 0);
  }
  const channel = new MessageChannel();
  const queue: (() => void)[] = [];

  channel.port1.onmessage = () => {
    for (const fn of queue.splice(0)) fn();
  };
  // Node: a port with a listener would keep the process alive
  (channel.port1 as { unref?: () => void }).unref?.();

  return (fn) => {
    queue.push(fn);
    channel.port2.postMessage(0);
  };
})();

/** Keys are absolute paths from the store's root. */
function norm(key: string): string {
  return key.startsWith("/") ? key : `/${key}`;
}
