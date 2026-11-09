import { describe, expect, it } from "vitest";

import { nextSelection, prune } from "@/lib/utils/legend-selection";

/**
 * What a click on a legend chip does to the set picked for recolouring
 * (lib/utils/legend-selection.ts). A plain click must stay the chip's own
 * action — opening its colour picker — or the legend loses what it did
 * before multi-select existed.
 */

const KEYS = ["a", "b", "c", "d", "e"];
const PLAIN = { shiftKey: false, metaKey: false, ctrlKey: false };
const SHIFT = { ...PLAIN, shiftKey: true };
const META = { ...PLAIN, metaKey: true };
const CTRL = { ...PLAIN, ctrlKey: true };

const keysOf = (s: Set<string>) => [...s].sort();

describe("a plain click", () => {
  it("is left to the chip, and clears the selection", () => {
    const step = nextSelection(KEYS, new Set(["a", "b"]), "a", "c", PLAIN);

    expect(step.handled).toBe(false);
    expect(step.selected.size).toBe(0);
    expect(step.anchor).toBe("c");
  });
});

describe("⌘ or ctrl-click", () => {
  it("adds a chip, and takes it out again", () => {
    const added = nextSelection(KEYS, new Set(["a"]), "a", "c", META);

    expect(added.handled).toBe(true);
    expect(keysOf(added.selected)).toEqual(["a", "c"]);

    const removed = nextSelection(KEYS, added.selected, "c", "a", CTRL);

    expect(keysOf(removed.selected)).toEqual(["c"]);
  });
});

describe("shift-click", () => {
  it("takes every chip between the anchor and this one", () => {
    const step = nextSelection(KEYS, new Set(["b"]), "b", "d", SHIFT);

    expect(step.handled).toBe(true);
    expect(keysOf(step.selected)).toEqual(["b", "c", "d"]);
    // The anchor stays put, so widening the range from the same end works
    expect(step.anchor).toBe("b");
    expect(keysOf(nextSelection(KEYS, step.selected, step.anchor, "e", SHIFT).selected)).toEqual(
      ["b", "c", "d", "e"],
    );
  });

  it("works backwards, from a later chip to an earlier one", () => {
    const step = nextSelection(KEYS, new Set(["d"]), "d", "b", SHIFT);

    expect(keysOf(step.selected)).toEqual(["b", "c", "d"]);
  });

  it("with no anchor yet, picks just the chip clicked", () => {
    const step = nextSelection(KEYS, new Set(), null, "c", SHIFT);

    expect(keysOf(step.selected)).toEqual(["c"]);
  });

  it("falls back to a plain click when the anchor has gone", () => {
    const step = nextSelection(KEYS, new Set(["z"]), "z", "c", SHIFT);

    expect(step.handled).toBe(false);
    expect(step.selected.size).toBe(0);
  });
});

describe("chips that leave the legend", () => {
  it("drop out of the selection", () => {
    expect(keysOf(prune(new Set(["a", "z", "c"]), KEYS))).toEqual(["a", "c"]);
  });

  it("leave an empty set when they all go", () => {
    expect(prune(new Set(["y", "z"]), KEYS).size).toBe(0);
  });
});
