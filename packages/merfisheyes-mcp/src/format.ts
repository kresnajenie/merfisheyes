import * as z from "zod/v4";

import { CHARACTER_LIMIT } from "./config.js";

/**
 * How results are returned: markdown for the model to read (default) or JSON
 * (as structuredContent plus its text). A list that would pass
 * CHARACTER_LIMIT returns fewer items — never cut text — and next_offset then
 * points at the first one left out.
 */

export const responseFormat = z
  .enum(["markdown", "json"])
  .default("markdown")
  .describe(
    "'markdown' for readable text (default), 'json' for structured data",
  );

export const paging = {
  limit: z
    .number()
    .int()
    .min(1)
    .max(100)
    .default(20)
    .describe("How many to return (1–100)"),
  offset: z
    .number()
    .int()
    .min(0)
    .default(0)
    .describe("How many to skip; use next_offset from the previous page"),
};

export interface Page<T> {
  total: number;
  count: number;
  offset: number;
  has_more: boolean;
  next_offset?: number;
  items: T[];
}

type Result = {
  content: { type: "text"; text: string }[];
  structuredContent?: Record<string, unknown>;
};

export function text(t: string): Result {
  return { content: [{ type: "text", text: t }] };
}

/** One object, as JSON or as the markdown `md` builds. */
export function one(
  format: "markdown" | "json",
  data: object,
  md: () => string,
): Result {
  if (format === "json") {
    return {
      content: [{ type: "text", text: JSON.stringify(data, null, 2) }],
      structuredContent: data as Record<string, unknown>,
    };
  }

  return text(md());
}

/** A page of items, fitted under the limit. `heading` and `item` build the
 * markdown; `empty` is said when there is nothing. */
export function many<T>(
  format: "markdown" | "json",
  page: Page<T>,
  md: { heading: string; item: (t: T) => string; empty: string },
): Result {
  const render = (items: T[]) => {
    const next = page.offset + items.length;
    const hasMore = page.total > next;
    const out: Page<T> = {
      ...page,
      items,
      count: items.length,
      has_more: hasMore,
    };

    if (hasMore) out.next_offset = next;
    else delete out.next_offset;
    if (format === "json") return { out, body: JSON.stringify(out, null, 2) };
    const lines = [
      `# ${md.heading}`,
      `Showing ${items.length} of ${page.total}${page.offset ? ` from ${page.offset + 1}` : ""}.`,
      "",
    ];

    lines.push(items.map(md.item).join("\n\n"));
    if (hasMore) lines.push("", `More: call again with offset=${next}.`);

    return { out, body: lines.join("\n") };
  };
  let items = page.items;
  let r = render(items);

  while (items.length > 1 && r.body.length > CHARACTER_LIMIT) {
    items = items.slice(0, Math.ceil(items.length / 2));
    r = render(items);
  }
  if (format === "json")
    return {
      content: [{ type: "text", text: r.body }],
      structuredContent: r.out as unknown as Record<string, unknown>,
    };
  if (page.total === 0) return text(md.empty);

  return text(r.body);
}

/** "12,345" — null and undefined become "—". */
export const n = (v: number | null | undefined) =>
  v === null || v === undefined ? "—" : v.toLocaleString("en-US");

/** "a · b · c" of the parts that are present. */
export const join = (...parts: (string | null | undefined | false)[]) =>
  parts.filter(Boolean).join(" · ");
