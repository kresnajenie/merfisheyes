import { describe, expect, it } from "vitest";

import { matchMoleculeColumns } from "@/lib/ingest/molecule-columns";

describe("matchMoleculeColumns", () => {
  it("reads a lab's own table by its column names", () => {
    // Beaumont4GBM's merged parquet
    expect(
      matchMoleculeColumns([
        "iGene",
        "X",
        "Y",
        "Z",
        "Cell ID",
        "Gene",
        "RNA_iset",
        "iset",
        "ifov",
      ]),
    ).toEqual({ gene: "Gene", x: "X", y: "Y", z: "Z" });
  });

  it("recognises the platforms' names too", () => {
    expect(
      matchMoleculeColumns(["feature_name", "x_location", "y_location"]),
    ).toEqual({ gene: "feature_name", x: "x_location", y: "y_location" });
    expect(
      matchMoleculeColumns(["gene", "global_x", "global_y", "global_z"]),
    ).toEqual({ gene: "gene", x: "global_x", y: "global_y", z: "global_z" });
  });

  it("is not fooled by a cell table", () => {
    // MERSCOPE's cell_metadata.csv: coordinates but no gene column
    expect(
      matchMoleculeColumns(["EntityID", "center_x", "center_y", "volume"]),
    ).toBeNull();
  });
});
