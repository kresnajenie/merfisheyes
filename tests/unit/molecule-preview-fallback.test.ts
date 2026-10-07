import { describe, it, expect } from "vitest";

import { readMoleculePreview } from "@/lib/services/molecule-preview";
import { pickSchema, readCsvHeader } from "@/lib/services/molecule-file-sniffer";

const mk = (content: string, name = "test_sm.csv") =>
  new File([content], name, { type: "text/csv" });

describe("molecule preview falls back to a header-only read", () => {
  // PapaParse's File streaming is unavailable in this environment, which is
  // exactly the "no preview" condition the fallback exists for.
  it("still reports columns when the streaming preview yields nothing", async () => {
    const p = await readMoleculePreview(mk("gene,global_x,global_y,cell_id\nANJING,1,2,-1\n"));

    expect(p.columns).toEqual(["gene", "global_x", "global_y", "cell_id"]);
    expect(p.autoType).toBe("merscope");
  });

  it("handles a ragged file (extra / short rows)", async () => {
    const p = await readMoleculePreview(mk("gene,global_x,global_y,cell_id\nxxx,1,2,3,-1\nnnn\n"));

    expect(p.columns).toEqual(["gene", "global_x", "global_y", "cell_id"]);
    expect(p.autoType).toBe("merscope");
  });

  it("handles CRLF, BOM and quoted headers", async () => {
    for (const csv of [
      "gene,global_x,global_y,cell_id\r\nA,1,2,-1\r\n",
      "﻿gene,global_x,global_y,cell_id\nA,1,2,-1\n",
      '"gene","global_x","global_y","cell_id"\nA,1,2,-1\n',
    ]) {
      const p = await readMoleculePreview(mk(csv));

      expect(p.columns).toEqual(["gene", "global_x", "global_y", "cell_id"]);
      expect(p.autoType).toBe("merscope");
    }
  });

  it("detects xenium headers too", async () => {
    const p = await readMoleculePreview(
      mk("feature_name,x_location,y_location,z_location\nA,1,2,3\n"),
    );

    expect(p.autoType).toBe("xenium");
  });

  it("header-only file still yields columns", async () => {
    const p = await readMoleculePreview(mk("gene,global_x,global_y,cell_id\n"));

    expect(p.columns).toHaveLength(4);
    expect(p.autoType).toBe("merscope");
  });

  it("readCsvHeader reads only the first 16 KB", async () => {
    const big = "gene,global_x,global_y,cell_id\n" + "A,1,2,-1\n".repeat(500_000);
    const cols = await readCsvHeader(mk(big));

    expect(cols).toEqual(["gene", "global_x", "global_y", "cell_id"]);
    expect(pickSchema(cols)).toBe("merscope");
  });
});
