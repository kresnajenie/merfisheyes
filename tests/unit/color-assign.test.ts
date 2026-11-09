import { describe, expect, it } from "vitest";

import { DEFAULT_COLOR_PALETTE } from "@/lib/utils/color-palette";
import { shuffleColors, spreadColors } from "@/lib/utils/color-assign";

/**
 * Giving several legend chips distinct colours again
 * (lib/utils/color-assign.ts), after they were all set to one.
 */

const palette = new Set<string>(DEFAULT_COLOR_PALETTE);

describe("spreading colours over a selection", () => {
  it("gives every item a different colour from the palette", () => {
    const out = spreadColors(["a", "b", "c"]);
    const colors = Object.values(out);

    expect(Object.keys(out)).toEqual(["a", "b", "c"]);
    expect(new Set(colors).size).toBe(3);
    for (const c of colors) expect(palette.has(c)).toBe(true);
  });

  it("keeps away from colours other items already use", () => {
    const taken = DEFAULT_COLOR_PALETTE.slice(0, 5);
    const out = spreadColors(["a", "b"], taken);

    for (const c of Object.values(out)) {
      expect(taken).not.toContain(c);
    }
  });

  it("ignores the case a colour is written in", () => {
    const first = DEFAULT_COLOR_PALETTE[0];
    const out = spreadColors(["a"], [first.toUpperCase()]);

    expect(out.a.toLowerCase()).not.toBe(first.toLowerCase());
  });

  it("ignores an eight-digit colour's alpha, and anything not a colour", () => {
    const first = DEFAULT_COLOR_PALETTE[0];
    const out = spreadColors(["a"], [`${first}80`, "", "not a colour"]);

    expect(out.a.toLowerCase()).not.toBe(first.toLowerCase());
  });

  it("reuses the palette rather than running out, past 40 items", () => {
    const keys = Array.from({ length: 45 }, (_, i) => `g${i}`);
    const out = spreadColors(keys);

    expect(Object.keys(out)).toHaveLength(45);
    for (const c of Object.values(out)) expect(palette.has(c)).toBe(true);
  });

  it("falls back to the whole palette when too much is taken", () => {
    const out = spreadColors(["a", "b"], DEFAULT_COLOR_PALETTE);

    expect(new Set(Object.values(out)).size).toBe(2);
  });
});

describe("shuffling a selection's colours", () => {
  it("gives every item a different colour from the palette", () => {
    const out = shuffleColors(["a", "b", "c", "d"]);

    expect(new Set(Object.values(out)).size).toBe(4);
    for (const c of Object.values(out)) expect(palette.has(c)).toBe(true);
  });

  it("does not always land on the same arrangement", () => {
    const keys = Array.from({ length: 6 }, (_, i) => `g${i}`);
    const runs = new Set(
      Array.from({ length: 20 }, () =>
        JSON.stringify(Object.values(shuffleColors(keys))),
      ),
    );

    expect(runs.size).toBeGreaterThan(1);
  });
});
