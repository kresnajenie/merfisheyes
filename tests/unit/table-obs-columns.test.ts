import { describe, expect, it } from "vitest";

import { buildObsColumnsFromRows } from "@/lib/adapters/table-obs-columns";

/**
 * MERSCOPE's cell_metadata.csv starts with a column that has no name: the
 * file's row index. It is not a property of the cells, so the preview must
 * not offer it as a colour (docs/plans/PREVIEW-SAVE.md, step 0).
 */
describe("columns built from a row table", () => {
  const rows = [
    { "": "1001", center_x: "1", center_y: "2", leiden: "0", "Unnamed: 0": "7" },
    { "": "1002", center_x: "3", center_y: "4", leiden: "1", "Unnamed: 0": "8" },
  ];

  it("leaves out a column with no name, and pandas' Unnamed: N", () => {
    const columns = buildObsColumnsFromRows(rows, {
      coordinateKeys: ["center_x", "center_y"],
    }).map((c) => c.column);

    expect(columns).toEqual(["leiden"]);
  });
});
