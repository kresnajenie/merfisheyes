import { describe, expect, it, vi } from "vitest";

vi.mock("@aws-sdk/client-batch", () => ({
  BatchClient: class {},
  SubmitJobCommand: class {},
}));

const { estimatePeakGb, tierForUpload } = await import(
  "@/lib/ingest/compute-tier"
);
const { statsFromPrepare, usesPreparePipeline } = await import(
  "@/lib/ingest/pipeline"
);

const GB = 1024 ** 3;

/**
 * The tier a prepare.py run is sent to, from what was uploaded: the runs
 * measured in docs/plans/ONE-PIPELINE.md §1 each land on a tier they peak
 * under 75% of.
 */
describe("tierForUpload", () => {
  it("keeps a Xenium export without images on the standard tier", () => {
    expect(
      tierForUpload([
        { key: "lung/transcripts.parquet", size: 20e6 },
        { key: "lung/cell_feature_matrix.h5", size: 3e6 },
        { key: "lung/cells.parquet", size: 1e6 },
      ]),
    ).toBe("standard");
  });

  it("follows the Beaumont molecule tables", () => {
    // set4 peaked at 8.5 GB, set3 at 12.8 GB (over 75% of 16), set1 with its
    // table at 15.7 GB
    expect(tierForUpload([{ key: "set4.parquet", size: 2302240139 }])).toBe(
      "standard",
    );
    expect(tierForUpload([{ key: "set3.parquet", size: 3403358700 }])).toBe(
      "large",
    );
    expect(
      tierForUpload([
        { key: "set1/set1.h5ad", size: 457e6 },
        { key: "set1/set1__single_molecules.parquet", size: 6069179850 },
      ]),
    ).toBe("large");
  });

  it("counts a MERSCOPE cell table more heavily than its transcripts CSV", () => {
    // colon: 704k cells, 20 GB of transcripts, 12.5 GB peak
    const colon = [
      { key: "colon/cell_by_gene.csv", size: 1.2e9 },
      { key: "colon/cell_metadata.csv", size: 0.25e9 },
      { key: "colon/detected_transcripts.csv", size: 20e9 },
    ];

    expect(estimatePeakGb(colon)).toBeGreaterThan(12);
    expect(tierForUpload(colon)).toBe("large");
  });

  it("never goes below the owner's own tier, and tops out at xlarge", () => {
    expect(tierForUpload([{ key: "a.h5ad", size: 1e6 }], "large")).toBe(
      "large",
    );
    expect(tierForUpload([{ key: "huge.parquet", size: 40 * GB }])).toBe(
      "xlarge",
    );
  });
});

describe("prepare.py's stats.json", () => {
  it("gives cells and genes, or molecules and their genes", () => {
    expect(statsFromPrepare({ cells: 10, genes: 5, molecules: 99 })).toEqual({
      numCells: 10,
      numGenes: 5,
    });
    expect(
      statsFromPrepare({ cells: 0, genes: 0, molecules: 99, molecule_genes: 7 }),
    ).toEqual({ numCells: 99, numGenes: 7 });
  });

  it("is what a dataset uploaded under the prepare pipeline expects", () => {
    expect(usesPreparePipeline({ kind: "single_cell", pipeline: "prepare" })).toBe(
      true,
    );
    expect(usesPreparePipeline({ kind: "single_cell" })).toBe(false);
    expect(usesPreparePipeline(null)).toBe(false);
  });
});
