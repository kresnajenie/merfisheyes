import { readFileSync } from "fs";
import path from "path";

import { describe, expect, it } from "vitest";

/**
 * The agent API's pure rules (docs/specs/MCP-PROTOCOL.md §4, §6): resources,
 * scopes, the fetch guard, sidecar files, paging. The database-backed
 * checks are tests/unit/agent-oauth-db.test.ts.
 */

process.env.OAUTH_ISSUER = "https://meyes.example";

const { parseScopes, resourceMetadataUrl, resources } = await import("@/lib/oauth/config");
const { isPrivateAddress, assertPublicHttps, maxAgeSeconds } = await import("@/lib/net/limited-fetch");
const { sidecarBase, isS3Host, summarise } = await import("@/lib/agent/sidecars");
const { page, paging } = await import("@/lib/agent/respond");

describe("resources and scopes", () => {
  it("names the agent API as the resource, canonically", () => {
    expect(resources()).toEqual(["https://meyes.example/api/agent/v1"]);
    expect(resourceMetadataUrl()).toBe("https://meyes.example/.well-known/oauth-protected-resource/api/agent/v1");
  });

  it("reports unknown scopes instead of dropping them", () => {
    expect(parseScopes("projects:read  bogus projects:read")).toEqual({
      scopes: ["projects:read"],
      unknown: ["bogus"],
    });
  });
});

describe("the fetch guard", () => {
  it("knows private addresses", () => {
    for (const ip of ["127.0.0.1", "10.1.2.3", "172.16.0.1", "192.168.1.1", "169.254.169.254", "100.64.0.1", "0.0.0.0", "::1", "fd00::1", "fe80::1", "::ffff:127.0.0.1"]) {
      expect(isPrivateAddress(ip), ip).toBe(true);
    }
    for (const ip of ["52.218.1.1", "8.8.8.8", "2600:1f18::1"]) {
      expect(isPrivateAddress(ip), ip).toBe(false);
    }
  });

  it("refuses http, credentials, private literals and disallowed hosts", async () => {
    await expect(assertPublicHttps(new URL("http://example.com/x"), {})).rejects.toThrow(/https/);
    await expect(assertPublicHttps(new URL("https://a:b@example.com/x"), {})).rejects.toThrow(/credentials/);
    await expect(assertPublicHttps(new URL("https://127.0.0.1/x"), {})).rejects.toThrow(/public/);
    await expect(assertPublicHttps(new URL("https://[::1]/x"), {})).rejects.toThrow(/public/);
    await expect(
      assertPublicHttps(new URL("https://8.8.8.8/x"), { allowHost: isS3Host }),
    ).rejects.toThrow(/not allowed/);
  });

  it("clamps cache lifetimes", () => {
    const h = (v: string) => new Headers({ "cache-control": v });

    expect(maxAgeSeconds(h("max-age=999999"), 3600, 300, 86400)).toBe(86400);
    expect(maxAgeSeconds(h("max-age=10"), 3600, 300, 86400)).toBe(300);
    expect(maxAgeSeconds(h("no-store"), 3600, 300, 86400)).toBe(300);
    expect(maxAgeSeconds(new Headers(), 3600, 300, 86400)).toBe(3600);
  });
});

describe("sidecar files", () => {
  it("looks beside the store, on S3 only", () => {
    const store = "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/bil-psc-zarr/ace-dip-vat/ace_dip_vat.zarr";

    expect(sidecarBase(store)).toBe("https://merfisheyes-bil.s3.us-west-2.amazonaws.com/bil-psc-zarr/ace-dip-vat/");
    expect(sidecarBase(`${store}/`)).toBe(sidecarBase(store));
    expect(sidecarBase("https://bucket.s3.amazonaws.com/legacy/meyes_output")).toBe(
      "https://bucket.s3.amazonaws.com/legacy/meyes_output/",
    );
    expect(isS3Host("merfisheyes-bil.s3.us-west-2.amazonaws.com")).toBe(true);
    expect(isS3Host("s3.amazonaws.com")).toBe(true);
    expect(isS3Host("evil.example")).toBe(false);
    expect(isS3Host("amazonaws.com.evil.example")).toBe(false);
  });

  it("summarises the fixtures without guessing", () => {
    const dir = path.resolve(__dirname, "../fixtures/sidecars");
    const read = (code: string, n: string) => JSON.parse(readFileSync(path.join(dir, code, `${n}.json`), "utf8"));
    const vat = summarise({ base: "", catalog: read("ace-dip-vat", "catalog"), stats: read("ace-dip-vat", "stats"), qc: read("ace-dip-vat", "qc") });

    expect(vat.cells).toBe(474986);
    expect(vat.molecules).toBe(646439087);
    expect(vat.sections).toBe(4);
    expect(vat.deColumns).toEqual(["class_name", "subclass_name"]);
    expect(vat.taxonomy?.levels.map((l) => l.level)).toEqual(["cluster", "class", "subclass", "supertype"]);
    expect(vat.taxonomy?.levels.find((l) => l.level === "class")?.groups).toBe(34);

    const owl = summarise({ base: "", catalog: read("ace-owl-cup", "catalog"), stats: read("ace-owl-cup", "stats"), qc: read("ace-owl-cup", "qc") });

    expect(owl.taxonomy).toBeNull();
    expect(owl.qc?.pctArtifact).toBeNull();
    expect(summarise({ base: "", catalog: null, stats: null, qc: null }).cells).toBeNull();
  });
});

describe("paging", () => {
  it("clamps limit and offset and says where the next page starts", () => {
    const p = paging(new URL("https://x/?limit=500&offset=-3"));

    expect(p).toEqual({ limit: 100, offset: 0 });
    expect(page([1, 2], 5, { limit: 2, offset: 2 })).toMatchObject({ total: 5, count: 2, has_more: true, next_offset: 4 });
    expect(page([5], 5, { limit: 2, offset: 4 })).toMatchObject({ has_more: false });
    expect(page([5], 5, { limit: 2, offset: 4 })).not.toHaveProperty("next_offset");
  });
});
