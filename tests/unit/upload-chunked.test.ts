import { describe, it, expect, vi, afterEach } from "vitest";

import {
  uploadChunkedToS3,
  type ChunkedUploadFile,
} from "@/lib/upload/uploadChunkedToS3";

/**
 * uploadChunkedToS3 against a fake server: one POST policy for the whole
 * prefix, batched UploadFile completion, manifest last, policy refresh.
 */

const S3_URL = "https://s3.test/bucket";

function file(key: string): ChunkedUploadFile {
  const blob = new Blob([key]);

  return { key, blob, size: blob.size, contentType: "application/gzip" };
}

function policy(expiresInMs: number) {
  return {
    url: S3_URL,
    fields: { key: "datasets/ds_1/${filename}", Policy: "p" },
    keyPrefix: "datasets/ds_1/",
    expiresAt: new Date(Date.now() + expiresInMs).toISOString(),
  };
}

const json = (body: unknown, status = 200) =>
  new Response(JSON.stringify(body), { status });

interface FakeServer {
  s3Keys: string[];
  completeBatches: string[][];
  refreshes: number;
  initBody: any;
}

/** Stub fetch; `s3Status(n)` decides the status of the n-th S3 POST (0-based). */
function stubServer(opts: {
  initialExpiryMs?: number;
  s3Status?: (n: number) => number;
}): FakeServer {
  const server: FakeServer = {
    s3Keys: [],
    completeBatches: [],
    refreshes: 0,
    initBody: null,
  };
  let s3Calls = 0;

  vi.stubGlobal(
    "fetch",
    vi.fn(async (url: string, init: RequestInit) => {
      if (url === "/api/datasets/initiate") {
        server.initBody = JSON.parse(init.body as string);

        return json({
          datasetId: "ds_1",
          uploadId: "up_1",
          policy: policy(opts.initialExpiryMs ?? 3600_000),
        });
      }
      if (url === "/api/datasets/ds_1/upload-policy") {
        server.refreshes++;

        return json({ policy: policy(3600_000) });
      }
      if (url === "/api/datasets/ds_1/files/complete") {
        server.completeBatches.push(JSON.parse(init.body as string).fileKeys);

        return json({ success: true });
      }
      if (url === S3_URL) {
        const status = opts.s3Status?.(s3Calls++) ?? 204;
        const fd = init.body as FormData;

        expect(fd.getAll("key")).toHaveLength(1);
        if (status < 300) server.s3Keys.push(fd.get("key") as string);

        return new Response(status < 300 ? null : "AccessDenied", { status });
      }
      throw new Error(`unexpected fetch ${url}`);
    }),
  );

  return server;
}

const META = { title: "t", numCells: 10, numGenes: 3 };

afterEach(() => {
  vi.unstubAllGlobals();
});

describe("uploadChunkedToS3", () => {
  it("uploads every file under the prefix, manifest last, completes in batches", async () => {
    const server = stubServer({});
    const files = [
      file("manifest.json"),
      ...Array.from({ length: 250 }, (_, i) => file(`expr/chunk_${i}.bin.gz`)),
    ];

    const res = await uploadChunkedToS3({
      fingerprint: "fp",
      metadata: META,
      files,
    });

    expect(res).toEqual({ datasetId: "ds_1", uploadId: "up_1" });
    expect(server.initBody.uploadMode).toBe("post-policy");
    expect(server.initBody.files).toHaveLength(251);

    expect(server.s3Keys).toHaveLength(251);
    expect(server.s3Keys.at(-1)).toBe("datasets/ds_1/manifest.json");
    expect(new Set(server.s3Keys).size).toBe(251);

    const completed = server.completeBatches.flat();

    expect(completed).toHaveLength(251);
    expect(completed.at(-1)).toBe("manifest.json");
    expect(Math.max(...server.completeBatches.map((b) => b.length))).toBe(100);
    expect(server.refreshes).toBe(0);
  });

  it("re-signs the policy on a 403 and retries the file", async () => {
    const server = stubServer({ s3Status: (n) => (n === 0 ? 403 : 204) });

    await uploadChunkedToS3({
      fingerprint: "fp",
      metadata: META,
      files: [file("manifest.json"), file("coords/spatial.bin.gz")],
    });

    expect(server.refreshes).toBe(1);
    expect(server.s3Keys).toEqual([
      "datasets/ds_1/coords/spatial.bin.gz",
      "datasets/ds_1/manifest.json",
    ]);
  });

  it("refreshes a policy that is about to expire before uploading", async () => {
    const server = stubServer({ initialExpiryMs: 60_000 });

    await uploadChunkedToS3({
      fingerprint: "fp",
      metadata: META,
      files: [
        file("manifest.json"),
        ...Array.from({ length: 20 }, (_, i) => file(`expr/c${i}.bin.gz`)),
      ],
    });

    // Eight concurrent workers share one refresh.
    expect(server.refreshes).toBe(1);
    expect(server.s3Keys).toHaveLength(21);
  });

  it("rejects an upload without a manifest before initiating", async () => {
    const server = stubServer({});

    await expect(
      uploadChunkedToS3({
        fingerprint: "fp",
        metadata: META,
        files: [file("expr/chunk_0.bin.gz")],
      }),
    ).rejects.toThrow("manifest.json");
    expect(server.initBody).toBeNull();
  });
});
