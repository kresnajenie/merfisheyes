import { describe, expect, it } from "vitest";

import { PART_CONCURRENCY, uploadParts } from "@/lib/upload/uploadRawToS3";

/**
 * A big raw file's parts go up several at once (a lab's single connection to
 * S3 carried 3.3 MB/s; many carry the link's 26–29 MB/s).
 */

const MB = 1024 * 1024;

function target(parts: number, partSize = MB) {
  return {
    partSize,
    parts: Array.from({ length: parts }, (_, i) => ({
      partNumber: i + 1,
      url: `https://s3.test/part-${i + 1}`,
    })),
  };
}

/** A PUT that takes a while, reports its bytes, and records concurrency. */
function fakePut(opts: { failOn?: string } = {}) {
  let inFlight = 0;
  const seen = { peak: 0, aborted: 0, started: [] as string[] };
  const put = (
    url: string,
    body: Blob,
    _type: string | null,
    onProgress: (loaded: number) => void,
    signal?: AbortSignal,
  ) =>
    new Promise<{ etag: string | null }>((resolve, reject) => {
      seen.started.push(url);
      inFlight++;
      seen.peak = Math.max(seen.peak, inFlight);
      // Stopped as XMLHttpRequest.abort() stops it: at once
      const onAbort = () => {
        clearTimeout(timer);
        inFlight--;
        seen.aborted++;
        reject(new DOMException("Aborted", "AbortError"));
      };
      const finish = () => {
        signal?.removeEventListener("abort", onAbort);
        inFlight--;
        if (url === opts.failOn) {
          reject(new Error("Upload failed with status 403"));
        } else {
          onProgress(body.size);
          // Later parts finish first, so order must come from part numbers
          resolve({ etag: `"etag-${url.split("-").pop()}"` });
        }
      };

      onProgress(body.size / 2);
      const timer = setTimeout(finish, 20 - (seen.started.length % 5) * 3);

      signal?.addEventListener("abort", onAbort);
    });

  return { put, seen };
}

describe("uploadParts", () => {
  it(`sends ${PART_CONCURRENCY} parts at a time and gives the ETags in part order`, async () => {
    const { put, seen } = fakePut();
    const file = new Blob([new Uint8Array(20 * MB + 123)]);
    let done = 0;

    const etags = await uploadParts(
      target(21),
      file,
      "big.parquet",
      () => {},
      (bytes) => {
        done += bytes;
      },
      undefined,
      put,
    );

    expect(seen.peak).toBe(PART_CONCURRENCY);
    expect(etags.map((e) => e.partNumber)).toEqual(
      Array.from({ length: 21 }, (_, i) => i + 1),
    );
    expect(etags[0].etag).toBe('"etag-1"');
    // Every byte counted once, the short last part included
    expect(done).toBe(file.size);
  });

  it("reports the bytes of the parts in flight together, never more than the file", async () => {
    const { put } = fakePut();
    const totals: number[] = [];
    let done = 0;

    await uploadParts(
      target(3),
      new Blob([new Uint8Array(3 * MB)]),
      "t.parquet",
      (inFlight) => totals.push(done + inFlight),
      (bytes) => {
        done += bytes;
      },
      undefined,
      put,
    );

    // Three half-sent parts at once count together…
    expect(totals).toContain(1.5 * MB);
    // …and finished parts are not counted twice
    expect(Math.max(...totals)).toBe(3 * MB);
  });

  it("stops the other parts when one fails, with that part's reason", async () => {
    const { put, seen } = fakePut({ failOn: "https://s3.test/part-2" });

    await expect(
      uploadParts(
        target(12),
        new Blob([new Uint8Array(12 * MB)]),
        "t.parquet",
        () => {},
        () => {},
        undefined,
        put,
      ),
    ).rejects.toThrow("Upload failed with status 403");
    // The parts already going were stopped, and no new ones started
    expect(seen.aborted).toBeGreaterThan(0);
    expect(seen.started.length).toBeLessThan(12);
  });
});
