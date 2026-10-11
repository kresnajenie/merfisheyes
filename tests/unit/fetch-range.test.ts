import { afterEach, describe, expect, it, vi } from "vitest";

import { fetchRange } from "@/lib/storage/fetchRange";

const BODY = Uint8Array.from({ length: 10 }, (_, i) => i);

/** A server holding one 10-byte object. `honoursRange` false ignores Range. */
function serve(honoursRange: boolean, status?: number) {
  const seen: string[] = [];

  vi.stubGlobal("fetch", async (_url: string, init: RequestInit) => {
    const range = (init.headers as Record<string, string>).Range;

    seen.push(range);
    if (status) return new Response(null, { status });
    if (!honoursRange) return new Response(BODY);
    const suffix = /^bytes=-(\d+)$/.exec(range);
    const span = /^bytes=(\d+)-(\d+)$/.exec(range);
    const part = suffix
      ? BODY.slice(-Number(suffix[1]))
      : BODY.slice(Number(span![1]), Number(span![2]) + 1);

    return new Response(part, { status: 206 });
  });

  return seen;
}

afterEach(() => vi.unstubAllGlobals());

describe("fetchRange", () => {
  it("asks for an offset range and a suffix range", async () => {
    const seen = serve(true);

    expect(await fetchRange("u", { offset: 2, length: 3 })).toEqual(
      Uint8Array.of(2, 3, 4),
    );
    expect(await fetchRange("u", { suffixLength: 4 })).toEqual(
      Uint8Array.of(6, 7, 8, 9),
    );
    expect(seen).toEqual(["bytes=2-4", "bytes=-4"]);
  });

  it("cuts the range out when the server sends the whole object", async () => {
    serve(false);

    expect(await fetchRange("u", { offset: 2, length: 3 })).toEqual(
      Uint8Array.of(2, 3, 4),
    );
    expect(await fetchRange("u", { suffixLength: 4 })).toEqual(
      Uint8Array.of(6, 7, 8, 9),
    );
  });

  it("treats 404 and 403 as a missing object, other failures as errors", async () => {
    serve(true, 404);
    expect(await fetchRange("u", { suffixLength: 1 })).toBeUndefined();
    serve(true, 403);
    expect(await fetchRange("u", { suffixLength: 1 })).toBeUndefined();
    serve(true, 500);
    await expect(fetchRange("u", { suffixLength: 1 })).rejects.toThrow(/500/);
  });
});
