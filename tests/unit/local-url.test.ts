import { describe, expect, it } from "vitest";

import {
  isLocalUrl,
  keepPrivate,
  withoutLocalUrls,
} from "@/lib/utils/local-url";

describe("data on this computer", () => {
  it("knows a local address from a public one", () => {
    for (const u of [
      "http://127.0.0.1:52341/k3y/sample.zarr",
      "http://localhost:3000/data",
      "http://[::1]:8080/x",
      "http://127.1.2.3/x",
      "http://app.localhost/x",
    ])
      expect(isLocalUrl(u), u).toBe(true);
    for (const u of [
      "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/a",
      "https://www.merfisheyes.com/viewer",
      "http://127.example.com/x",
      "not a url",
      "",
      null,
    ])
      expect(isLocalUrl(u), String(u)).toBe(false);
  });

  it("drops local addresses from a page address before analytics sees it", () => {
    const local = "http://127.0.0.1:52341/k3y/my-sample.zarr";
    const page = `https://www.merfisheyes.com/sm-viewer/from-s3?url=${encodeURIComponent(local)}&splitS3Url=${encodeURIComponent(encodeURIComponent(local))}&splitType=sm&v=abc`;
    const sent = withoutLocalUrls(page);

    expect(sent).not.toContain("127.0.0.1");
    expect(sent).not.toContain("my-sample");
    expect(new URL(sent).searchParams.get("v")).toBe("abc");
    expect(new URL(sent).searchParams.get("splitType")).toBe("sm");

    const s3 = "https://www.merfisheyes.com/viewer/from-s3?url=https%3A%2F%2Fb.s3.amazonaws.com%2Fa";

    expect(withoutLocalUrls(s3)).toBe(s3);
  });
});

describe("keeping a local address from the site", () => {
  const local = "http://127.0.0.1:52341/k3y/sample.zarr";
  const s3 = "https://b.s3.amazonaws.com/a";

  it("keeps it when the page came from the internet", () => {
    expect(keepPrivate(local, "https://www.merfisheyes.com/open?url=x")).toBe(true);
    expect(keepPrivate(s3, "https://www.merfisheyes.com/viewer")).toBe(false);
  });

  it("treats it as any bucket when the page is on this computer too", () => {
    expect(keepPrivate(local, "http://localhost:3000/viewer")).toBe(false);
  });

  it("always keeps it in the standalone viewer", () => {
    const was = process.env.NEXT_PUBLIC_MERFISHEYES_STANDALONE;

    process.env.NEXT_PUBLIC_MERFISHEYES_STANDALONE = "1";
    try {
      expect(keepPrivate(local, "http://127.0.0.1:8765/open/")).toBe(true);
      expect(keepPrivate(s3, "http://127.0.0.1:8765/open/")).toBe(true);
    } finally {
      process.env.NEXT_PUBLIC_MERFISHEYES_STANDALONE = was;
    }
  });
});

