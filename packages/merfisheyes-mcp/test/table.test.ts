import { describe, expect, it } from "vitest";

import { idColumn, parseTable } from "../src/local/table.ts";

describe("an agent's results file", () => {
  it("reads CSV with quotes and pandas' unnamed index, and TSV", () => {
    const csv =
      ',pseudotime,fate\ncell_1,0.25,"radial glia"\ncell_2,0.9,"immature, late"\n';
    const t = parseTable(csv);

    expect(t.header).toEqual(["", "pseudotime", "fate"]);
    expect(t.rows[1]).toEqual(["cell_2", "0.9", "immature, late"]);
    expect(idColumn(t.header)).toBe(0);
    const tsv = parseTable("cell_id\tlabel\r\na\tx\r\nb\ty\r\n");

    expect(tsv.rows).toEqual([
      ["a", "x"],
      ["b", "y"],
    ]);
    expect(idColumn(tsv.header)).toBe(0);
  });

  it("finds the id column by name, or none (values in cell order)", () => {
    expect(idColumn(["score", "obs_names"])).toBe(1);
    expect(idColumn(["score", "label"])).toBeNull();
    expect(idColumn(["a", "b"], "b")).toBe(1);
    expect(() => idColumn(["a"], "nope")).toThrow(/No column "nope"/);
  });
});
