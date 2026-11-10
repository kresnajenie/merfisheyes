import { describe, expect, it } from "vitest";

import { CHARACTER_LIMIT } from "../src/config.js";
import { many } from "../src/format.js";
import { isDatasetId } from "../src/tools/catalog.js";

describe("pages", () => {
  const page = (items: string[], total = items.length, offset = 0) => ({
    total,
    count: items.length,
    offset,
    has_more: false,
    items,
  });
  const md = { heading: "H", item: (s: string) => s, empty: "nothing" };

  it("say where the next page starts", () => {
    const r = many("json", page(["a", "b"], 5, 2), md);

    expect(r.structuredContent).toMatchObject({ count: 2, has_more: true, next_offset: 4 });
    expect(JSON.parse(r.content[0].text)).toEqual(r.structuredContent);
    expect(many("markdown", page(["a"], 3, 2), md).content[0].text).not.toContain("offset=");
  });

  it("stay under the character limit by returning fewer items, never cut ones", () => {
    const big = Array.from({ length: 40 }, (_, i) => `${i}:${"x".repeat(2_000)}`);
    const r = many("markdown", page(big, 100), md);
    const text = r.content[0].text;

    expect(text.length).toBeLessThanOrEqual(CHARACTER_LIMIT);
    expect(text).toMatch(/More: call again with offset=(\d+)/);
    const shown = Number(/offset=(\d+)/.exec(text)![1]);

    for (let i = 0; i < shown; i++) expect(text).toContain(`${i}:${"x".repeat(2_000)}`);
  });

  it("say so when empty", () => {
    expect(many("markdown", page([]), md).content[0].text).toBe("nothing");
  });
});

it("tells app dataset ids from catalog ids", () => {
  expect(isDatasetId("ds_a1B2c3D4e5")).toBe(true);
  expect(isDatasetId("sm_x")).toBe(true);
  expect(isDatasetId("s3_x")).toBe(true);
  expect(isDatasetId("ace-dip-vat")).toBe(false);
  expect(isDatasetId("cm1xyz")).toBe(false);
});
