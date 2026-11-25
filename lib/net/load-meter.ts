/**
 * What the page is downloading right now, in bytes: one meter for every data
 * read (a store's chunks, a gene file, image tiles, cell outlines), so one bar
 * can say "18.2 MB / 80.4 MB" whatever is loading.
 *
 * Reads come in bursts: a burst starts with the first request in flight and
 * ends a moment after the last one finishes, and the counts are of the burst.
 * The total is the sum of the sizes the servers have announced so far
 * (Content-Length), so it is exact for the requests already answered and
 * grows as a view asks for more; a response with no announced size makes the
 * total unknown and the bar falls back to requests finished of requests made.
 */

export interface LoadSnapshot {
  /** A burst is under way. */
  active: boolean;
  /** Bytes received in this burst. */
  received: number;
  /** Bytes announced for the requests of this burst. */
  total: number;
  /** Every answered request announced its size. */
  totalKnown: boolean;
  /** Requests finished and made in this burst. */
  done: number;
  started: number;
  /** performance.now() when the burst began. */
  since: number;
}

/** How long with nothing in flight before a burst is over. */
const IDLE_MS = 400;
/** Listeners hear of byte progress at most this often. */
const NOTIFY_MS = 100;

const IDLE: LoadSnapshot = {
  active: false,
  received: 0,
  total: 0,
  totalKnown: true,
  done: 0,
  started: 0,
  since: 0,
};

let state: LoadSnapshot = IDLE;
let inFlight = 0;
let idleTimer: ReturnType<typeof setTimeout> | null = null;
let notifyTimer: ReturnType<typeof setTimeout> | null = null;
const listeners = new Set<() => void>();

function notify(now = false) {
  if (now) {
    if (notifyTimer) clearTimeout(notifyTimer);
    notifyTimer = null;
    listeners.forEach((l) => l());

    return;
  }
  if (notifyTimer) return;
  notifyTimer = setTimeout(() => {
    notifyTimer = null;
    listeners.forEach((l) => l());
  }, NOTIFY_MS);
}

/** One request of a burst; call `finish` exactly when its body is done with. */
export interface LoadHandle {
  /** The size the server announced, or null when it did not. */
  announce(bytes: number | null): void;
  add(bytes: number): void;
  finish(): void;
}

/** Set in a web worker: its reads are counted on the page that owns it. */
const inWorker =
  typeof document === "undefined" &&
  typeof (globalThis as { importScripts?: unknown }).importScripts ===
    "function";
const MESSAGE = "__loadMeter";

type WorkerEvent =
  | ["begin", number]
  | ["announce", number, number | null]
  | ["add", number, number]
  | ["finish", number];

let nextId = 1;

/** A handle in a worker: the same calls, sent to the page. */
function forwardedLoad(): LoadHandle {
  const id = nextId++;
  const send = (event: WorkerEvent) =>
    (self as unknown as Worker).postMessage({ [MESSAGE]: event });
  let pending = 0;
  let timer: ReturnType<typeof setTimeout> | null = null;
  const flush = () => {
    if (timer) clearTimeout(timer);
    timer = null;
    if (pending) send(["add", id, pending]);
    pending = 0;
  };

  send(["begin", id]);

  return {
    announce: (bytes) => send(["announce", id, bytes]),
    add(bytes) {
      pending += bytes;
      timer ??= setTimeout(flush, NOTIFY_MS);
    },
    finish() {
      flush();
      send(["finish", id]);
    },
  };
}

/**
 * Count a worker's reads on this page's meter. Call it on every worker that
 * fetches data; its messages carry their own key, which Comlink ignores.
 */
export function attachLoadMeter(worker: Worker): void {
  const loads = new Map<number, LoadHandle>();

  worker.addEventListener("message", (e: MessageEvent) => {
    const event = e.data?.[MESSAGE] as WorkerEvent | undefined;

    if (!event) return;
    const [kind, id] = event;

    if (kind === "begin") loads.set(id, beginLoad());
    else if (kind === "announce") loads.get(id)?.announce(event[2]);
    else if (kind === "add") loads.get(id)?.add(event[2]);
    else {
      loads.get(id)?.finish();
      loads.delete(id);
    }
  });
  // A worker that dies takes its unfinished reads with it
  worker.addEventListener("error", () => {
    loads.forEach((load) => load.finish());
    loads.clear();
  });
}

export function beginLoad(): LoadHandle {
  if (inWorker) return forwardedLoad();
  if (idleTimer) clearTimeout(idleTimer);
  idleTimer = null;
  state = state.active
    ? { ...state, started: state.started + 1 }
    : { ...IDLE, active: true, started: 1, since: performance.now() };
  inFlight++;
  notify(state.started === 1);
  let open = true;

  return {
    announce(bytes) {
      if (!open) return;
      state =
        bytes === null
          ? { ...state, totalKnown: false }
          : { ...state, total: state.total + bytes };
      notify();
    },
    add(bytes) {
      if (!open) return;
      state = { ...state, received: state.received + bytes };
      notify();
    },
    finish() {
      if (!open) return;
      open = false;
      inFlight--;
      state = { ...state, done: state.done + 1 };
      notify();
      if (inFlight > 0) return;
      idleTimer = setTimeout(() => {
        idleTimer = null;
        state = { ...state, active: false };
        notify(true);
      }, IDLE_MS);
    },
  };
}

export function loadSnapshot(): LoadSnapshot {
  return state;
}

export function subscribeLoad(listener: () => void): () => void {
  listeners.add(listener);

  return () => listeners.delete(listener);
}

/**
 * `fetch`, with the response's bytes counted on the meter as they arrive.
 * The caller reads the body as it would any other.
 */
export async function trackedFetch(
  input: RequestInfo | URL,
  init?: RequestInit,
): Promise<Response> {
  const load = beginLoad();
  let res: Response;

  try {
    res = await fetch(input, init);
  } catch (e) {
    load.finish();
    throw e;
  }
  // Nothing worth counting, and a body the caller may never read
  if (!res.ok || !res.body) {
    load.finish();

    return res;
  }
  const length = Number(res.headers.get("content-length"));

  // A compressed transfer announces the compressed size but delivers more
  load.announce(
    length > 0 && !res.headers.get("content-encoding") ? length : null,
  );
  const reader = res.body.getReader();
  const body = new ReadableStream<Uint8Array>({
    async pull(controller) {
      try {
        const { done, value } = await reader.read();

        if (done) {
          load.finish();
          controller.close();
        } else {
          load.add(value.byteLength);
          controller.enqueue(value);
        }
      } catch (e) {
        load.finish();
        controller.error(e);
      }
    },
    cancel(reason) {
      load.finish();

      return reader.cancel(reason);
    },
  });

  return new Response(body, {
    status: res.status,
    statusText: res.statusText,
    headers: res.headers,
  });
}

/** "18.2 MB", "1.24 GB", "640 KB". */
export function formatBytes(n: number): string {
  if (n >= 1e9) return `${(n / 1e9).toFixed(2)} GB`;
  if (n >= 1e6) return `${(n / 1e6).toFixed(1)} MB`;

  return `${Math.max(0, Math.round(n / 1e3))} KB`;
}

/** The bar's fill, 0 to 1, and its text. */
export function describeLoad(s: LoadSnapshot): {
  fraction: number;
  text: string;
} {
  const files = `${s.done} of ${s.started} ${s.started === 1 ? "file" : "files"}`;

  if (s.totalKnown && s.total > 0) {
    return {
      fraction: Math.min(1, s.received / s.total),
      text: `${formatBytes(s.received)} / ${formatBytes(s.total)} · ${files}`,
    };
  }

  return {
    fraction: s.started ? s.done / s.started : 0,
    text: `${formatBytes(s.received)} · ${files}`,
  };
}
