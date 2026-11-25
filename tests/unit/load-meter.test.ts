import { afterEach, beforeEach, describe, expect, it, vi } from "vitest";

import {
  beginLoad,
  describeLoad,
  formatBytes,
  type LoadHandle,
  loadSnapshot,
  subscribeLoad,
  trackedFetch,
} from "@/lib/net/load-meter";

// Every read a test opens is closed after it, so the next starts idle
let opened: LoadHandle[] = [];
const begin = () => {
  const load = beginLoad();

  opened.push(load);

  return load;
};

beforeEach(() => vi.useFakeTimers());
afterEach(() => {
  vi.useFakeTimers();
  opened.forEach((load) => load.finish());
  opened = [];
  vi.runAllTimers();
  vi.useRealTimers();
  vi.unstubAllGlobals();
});

describe("formatBytes", () => {
  it("picks the unit by size", () => {
    expect(formatBytes(640_000)).toBe("640 KB");
    expect(formatBytes(18_240_000)).toBe("18.2 MB");
    expect(formatBytes(1_240_000_000)).toBe("1.24 GB");
  });
});

describe("a burst of reads", () => {
  it("adds up what was announced and received, and ends when nothing is in flight", () => {
    const a = begin();
    const b = begin();

    a.announce(1_000_000);
    b.announce(3_000_000);
    a.add(1_000_000);
    b.add(500_000);
    a.finish();

    expect(loadSnapshot()).toMatchObject({
      active: true,
      received: 1_500_000,
      total: 4_000_000,
      totalKnown: true,
      done: 1,
      started: 2,
    });
    expect(describeLoad(loadSnapshot())).toEqual({
      fraction: 0.375,
      text: "1.5 MB / 4.0 MB · 1 of 2 files",
    });

    b.finish();
    expect(loadSnapshot().active).toBe(true);
    vi.advanceTimersByTime(500);
    expect(loadSnapshot().active).toBe(false);
  });

  it("starts from zero again after it has ended", () => {
    const a = begin();

    a.announce(10);
    a.add(10);
    a.finish();
    vi.advanceTimersByTime(500);

    begin();
    expect(loadSnapshot()).toMatchObject({ received: 0, total: 0, started: 1 });
  });

  it("a read that follows at once belongs to the same burst", () => {
    begin().finish();
    vi.advanceTimersByTime(100);
    begin();
    expect(loadSnapshot()).toMatchObject({ active: true, done: 1, started: 2 });
  });

  it("without an announced size, goes by files", () => {
    const a = begin();
    const b = begin();

    a.announce(null);
    a.add(2_000_000);
    a.finish();
    b.announce(1_000_000);

    expect(describeLoad(loadSnapshot())).toEqual({
      fraction: 0.5,
      text: "2.0 MB · 1 of 2 files",
    });
  });

  it("finishing twice counts once", () => {
    const a = begin();

    a.finish();
    a.finish();
    a.add(5);
    expect(loadSnapshot()).toMatchObject({ done: 1, received: 0 });
  });

  it("tells listeners when it starts and when it ends", () => {
    const heard = vi.fn();
    const stop = subscribeLoad(heard);

    const a = begin();

    expect(heard).toHaveBeenCalledTimes(1);
    a.finish();
    vi.runAllTimers();
    expect(loadSnapshot().active).toBe(false);
    expect(heard.mock.calls.length).toBeGreaterThanOrEqual(2);
    stop();
  });
});

describe("trackedFetch", () => {
  const respond = (body: Uint8Array, headers: Record<string, string>) =>
    vi.stubGlobal(
      "fetch",
      vi.fn(async () => new Response(body, { status: 200, headers })),
    );

  it("counts the body as it is read and hands the same bytes on", async () => {
    vi.useRealTimers();
    const body = new Uint8Array(5000).fill(7);

    respond(body, { "content-length": "5000" });
    const res = await trackedFetch("https://example.org/x");

    expect(res.status).toBe(200);
    expect(new Uint8Array(await res.arrayBuffer())).toEqual(body);
    expect(loadSnapshot()).toMatchObject({
      received: 5000,
      total: 5000,
      done: 1,
      started: 1,
    });
  });

  it("does not trust the size of a compressed transfer", async () => {
    vi.useRealTimers();
    respond(new Uint8Array(100), {
      "content-length": "40",
      "content-encoding": "gzip",
    });
    await (await trackedFetch("https://example.org/x")).arrayBuffer();
    expect(loadSnapshot().totalKnown).toBe(false);
  });

  it("is done with a response that is not ok without its body being read", async () => {
    vi.useRealTimers();
    vi.stubGlobal(
      "fetch",
      vi.fn(async () => new Response("no", { status: 403 })),
    );
    const res = await trackedFetch("https://example.org/x");

    expect(res.status).toBe(403);
    expect(loadSnapshot().done).toBe(loadSnapshot().started);
  });

  it("is done when the request fails", async () => {
    vi.useRealTimers();
    vi.stubGlobal(
      "fetch",
      vi.fn(async () => {
        throw new Error("offline");
      }),
    );
    await expect(trackedFetch("https://example.org/x")).rejects.toThrow(
      "offline",
    );
    expect(loadSnapshot().done).toBe(loadSnapshot().started);
  });
});
