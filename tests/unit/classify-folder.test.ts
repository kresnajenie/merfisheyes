import { describe, it, expect } from "vitest";

import {
  classifyFolder,
  datasetRelativeKey,
  scoreTranscriptFile,
} from "@/lib/ingest/classify-folder";

/**
 * Build a File whose webkitRelativePath is set, the way a directory picker
 * produces it. The property is read-only on a real File, so it is defined here.
 */
function pickedFile(relativePath: string, size = 1024): File {
  const file = new File(["x"], relativePath.split("/").pop() as string);

  Object.defineProperty(file, "webkitRelativePath", { value: relativePath });
  Object.defineProperty(file, "size", { value: size });

  return file;
}

const keysOf = (files?: { key: string }[]) =>
  (files ?? []).map((f) => f.key).sort();

describe("datasetRelativeKey", () => {
  it("strips the picked folder's own name", () => {
    // Without this the worker unpacks into raw/{id}/<FolderName>/… and the
    // processor's format detection, which only looks at the top level, sees
    // nothing but a directory.
    expect(datasetRelativeKey(pickedFile("MyRun/cell_metadata.csv"))).toBe(
      "cell_metadata.csv",
    );
  });

  it("keeps structure below the root", () => {
    expect(
      datasetRelativeKey(pickedFile("MyRun/cell_feature_matrix/features.tsv.gz")),
    ).toBe("cell_feature_matrix/features.tsv.gz");
  });

  it("falls back to the name for a single picked file", () => {
    expect(datasetRelativeKey(pickedFile("data.h5ad"))).toBe("data.h5ad");
  });
});

describe("scoreTranscriptFile", () => {
  it("scores a detected_transcripts file highest", () => {
    expect(scoreTranscriptFile("detected_transcripts.csv")).toBe(100);
  });

  it("accepts parquet, which the worker passes explicitly", () => {
    expect(scoreTranscriptFile("detected_transcripts.parquet")).toBe(100);
  });

  it("gives a transcripts-only name the alternative score", () => {
    expect(scoreTranscriptFile("transcripts.csv")).toBe(80);
  });

  it("does not match single-cell files", () => {
    expect(scoreTranscriptFile("cell_by_gene.csv")).toBe(0);
    expect(scoreTranscriptFile("cell_metadata.csv")).toBe(0);
    expect(scoreTranscriptFile("cells.csv")).toBe(0);
  });
});

describe("classifyFolder", () => {
  it("picks exactly the two files a MERSCOPE folder needs", () => {
    // Mirrors the real test folders, which hold precisely these two.
    const result = classifyFolder([
      pickedFile("Human_Ovarian_Cancer_150mb/cell_by_gene.csv"),
      pickedFile("Human_Ovarian_Cancer_150mb/cell_metadata.csv"),
    ]);

    expect(result.singleCell?.format).toBe("merscope");
    expect(keysOf(result.singleCell?.files)).toEqual([
      "cell_by_gene.csv",
      "cell_metadata.csv",
    ]);
    expect(result.singleMolecule).toBeUndefined();
    expect(result.ignored).toHaveLength(0);
  });

  it("drops macOS junk", () => {
    const result = classifyFolder([
      pickedFile("Run/cell_metadata.csv"),
      pickedFile("Run/cell_by_gene.csv"),
      pickedFile("Run/.DS_Store"),
      pickedFile("Run/._cell_metadata.csv"),
      pickedFile("Run/__MACOSX/whatever"),
    ]);

    expect(keysOf(result.singleCell?.files)).toEqual([
      "cell_by_gene.csv",
      "cell_metadata.csv",
    ]);
    expect(keysOf(result.ignored)).toEqual([
      ".DS_Store",
      "._cell_metadata.csv",
      "__MACOSX/whatever",
    ]);
  });

  it("returns both groups when a Xenium folder also has transcripts", () => {
    const result = classifyFolder([
      pickedFile("XeniumRun/cells.csv.gz"),
      pickedFile("XeniumRun/cell_feature_matrix/features.tsv.gz"),
      pickedFile("XeniumRun/cell_feature_matrix/barcodes.tsv.gz"),
      pickedFile("XeniumRun/cell_feature_matrix/matrix.mtx.gz"),
      pickedFile("XeniumRun/transcripts.parquet"),
      pickedFile("XeniumRun/morphology.ome.tif"),
      pickedFile("XeniumRun/analysis_summary.html"),
    ]);

    expect(result.singleCell?.format).toBe("xenium");
    expect(keysOf(result.singleCell?.files)).toEqual([
      "cell_feature_matrix/barcodes.tsv.gz",
      "cell_feature_matrix/features.tsv.gz",
      "cell_feature_matrix/matrix.mtx.gz",
      "cells.csv.gz",
    ]);
    expect(keysOf(result.singleMolecule?.files)).toEqual(["transcripts.parquet"]);

    // The bulky, unread files are exactly what this is for.
    expect(keysOf(result.ignored)).toEqual([
      "analysis_summary.html",
      "morphology.ome.tif",
    ]);
  });

  it("treats a folder of h5ad files as single cell", () => {
    const result = classifyFolder([
      pickedFile("Batch/sample_a.h5ad"),
      pickedFile("Batch/notes.txt"),
    ]);

    expect(result.singleCell?.format).toBe("h5ad");
    expect(keysOf(result.singleCell?.files)).toEqual(["sample_a.h5ad"]);
    expect(keysOf(result.ignored)).toEqual(["notes.txt"]);
  });

  it("finds only single-molecule data when there is no cell-level file", () => {
    const result = classifyFolder([
      pickedFile("Run/detected_transcripts.csv"),
      pickedFile("Run/readme.txt"),
    ]);

    expect(result.singleCell).toBeUndefined();
    expect(keysOf(result.singleMolecule?.files)).toEqual([
      "detected_transcripts.csv",
    ]);
  });

  it("prefers the better-scoring transcripts file when several match", () => {
    const result = classifyFolder([
      pickedFile("Run/transcripts.csv"),
      pickedFile("Run/detected_transcripts.csv"),
    ]);

    expect(keysOf(result.singleMolecule?.files)).toEqual([
      "detected_transcripts.csv",
    ]);
  });

  it("reports nothing usable for a folder of unrelated files", () => {
    const result = classifyFolder([
      pickedFile("Run/image.tif"),
      pickedFile("Run/notes.txt"),
    ]);

    expect(result.singleCell).toBeUndefined();
    expect(result.singleMolecule).toBeUndefined();
    expect(result.ignored).toHaveLength(2);
  });
});
