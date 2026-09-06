import { describe, expect, it } from "vitest";

import {
  buildPaletteLut,
  buildSelectionLut,
  countVisible,
  hexToRgb255,
} from "@/lib/webgl/labelled-molecule-lut";
import { SIZE_ENCODE_RANGE } from "@/lib/webgl/labelled-molecule-shaders";

describe("buildSelectionLut", () => {
  const dict = ["a", "b", "c"];

  it("an empty menu is all-pass, not none-pass", () => {
    // The shader multiplies the three masks; an unconstrained menu has to
    // contribute 1.0 or the intersection collapses to nothing.
    expect([...buildSelectionLut(dict, new Set())]).toEqual([255, 255, 255]);
  });

  it("marks only the checked values", () => {
    expect([...buildSelectionLut(dict, new Set(["a", "c"]))]).toEqual([
      255, 0, 255,
    ]);
  });

  it("hides a value on top of the selection", () => {
    expect([
      ...buildSelectionLut(dict, new Set(["a", "c"]), new Set(["c"])),
    ]).toEqual([255, 0, 0]);
  });

  it("hides even when nothing is selected", () => {
    expect([...buildSelectionLut(dict, new Set(), new Set(["b"]))]).toEqual([
      255, 0, 255,
    ]);
  });

  it("ignores values that aren't in the column", () => {
    expect([...buildSelectionLut(dict, new Set(["zzz"]))]).toEqual([0, 0, 0]);
  });
});

describe("countVisible", () => {
  //   i:  0  1  2  3  4  5
  // gene: a  a  b  b  c  c
  // cell: X  Y  X  Y  X  Y
  const gene = { dict: ["a", "b", "c"], codes: [0, 0, 1, 1, 2, 2] };
  const cell = { dict: ["X", "Y"], codes: [0, 1, 0, 1, 0, 1] };
  const N = 6;

  const count = (g: string[], c: string[]) =>
    countVisible(
      [
        { indices: gene.codes, lut: buildSelectionLut(gene.dict, new Set(g)) },
        { indices: cell.codes, lut: buildSelectionLut(cell.dict, new Set(c)) },
      ],
      N,
    );

  it("shows everything when nothing is selected", () => {
    expect(count([], [])).toBe(6);
  });

  it("ORs within one menu", () => {
    expect(count(["a", "b"], [])).toBe(4);
  });

  it("ANDs across menus", () => {
    // gene ∈ {a,b} AND cell ∈ {X} -> rows 0 and 2
    expect(count(["a", "b"], ["X"])).toBe(2);
  });

  it("treats an empty menu as no constraint", () => {
    expect(count([], ["Y"])).toBe(3);
  });

  it("can intersect to nothing", () => {
    expect(count(["a"], [])).toBe(2);
    expect(
      countVisible(
        [
          {
            indices: gene.codes,
            lut: buildSelectionLut(gene.dict, new Set(["zzz"])),
          },
        ],
        N,
      ),
    ).toBe(0);
  });
});

describe("buildPaletteLut", () => {
  it("packs colour in RGB and size in A", () => {
    const lut = buildPaletteLut(["a", "b"], (v) => (v === "a" ? "#ff8800" : "#0000ff"), (v) =>
      v === "a" ? 2 : 1,
    );

    expect([...lut.slice(0, 3)]).toEqual([255, 136, 0]);
    expect([...lut.slice(4, 7)]).toEqual([0, 0, 255]);

    // A decodes back to the size multiplier within 8-bit precision.
    expect((lut[3] / 255) * SIZE_ENCODE_RANGE).toBeCloseTo(2, 1);
    expect((lut[7] / 255) * SIZE_ENCODE_RANGE).toBeCloseTo(1, 1);
  });

  it("clamps sizes to the encodable range", () => {
    const lut = buildPaletteLut(["a", "b"], () => "#ffffff", (v) =>
      v === "a" ? 999 : -5,
    );

    expect(lut[3]).toBe(255);
    expect(lut[7]).toBe(0);
  });
});

describe("hexToRgb255", () => {
  it("handles long and short hex", () => {
    expect(hexToRgb255("#00FF88")).toEqual([0, 255, 136]);
    expect(hexToRgb255("#f80")).toEqual([255, 136, 0]);
  });
});
