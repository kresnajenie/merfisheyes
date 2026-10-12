import { describe, expect, it } from "vitest";

import {
  classifyBundle,
  hasParts,
  partFileMap,
  partFiles,
} from "@/lib/local/classify-bundle";
import {
  fileStore,
  openSpatialDataStoreFromFiles,
} from "@/lib/spatialdata/localStore";

/** A File as a directory picker produces it, with optional contents. */
function picked(relativePath: string, contents = "x"): File {
  const file = new File([contents], relativePath.split("/").pop() as string);

  Object.defineProperty(file, "webkitRelativePath", {
    value: relativePath,
    configurable: true,
  });

  return file;
}

const image = (folder: string, axes = ["c", "y", "x"]) => [
  picked(
    `${folder}/zarr.json`,
    JSON.stringify({
      attributes: {
        ome: { multiscales: [{ axes: axes.map((name) => ({ name })) }] },
      },
    }),
  ),
  picked(`${folder}/0/zarr.json`, "{}"),
  picked(`${folder}/0/c/0/0/0`),
];
const cells = (folder: string) => [
  picked(`${folder}/manifest.json`, "{}"),
  picked(`${folder}/expr/index.json`),
  picked(`${folder}/expr/chunk_00000.bin.gz`),
  picked(`${folder}/coords/spatial.bin.gz`),
];
const molecules = (folder: string) => [
  picked(`${folder}/manifest.json.gz`),
  picked(`${folder}/genes/Gad1.bin.gz`),
];
const store = (folder: string) => [
  picked(`${folder}/zarr.json`, JSON.stringify({ attributes: {} })),
  picked(`${folder}/tables/table/zarr.json`, "{}"),
  picked(`${folder}/tables/table/obsm/spatial/zarr.json`, "{}"),
  picked(`${folder}/images/dapi/zarr.json`, "{}"),
];

describe("classifyBundle", () => {
  it("finds merfisheyes cells, molecules and an image by what is in them", async () => {
    // The folder names say nothing about what they hold
    const bundle = await classifyBundle([
      ...cells("s/one"),
      ...molecules("s/two"),
      ...image("s/three"),
      picked("s/notes.txt"),
    ]);

    expect(bundle).toMatchObject({
      name: "s",
      cells: { kind: "chunked", root: "one" },
      molecules: "two",
      image: "three",
      skipped: [],
    });
  });

  it("finds a SpatialData store, whose image is inside it", async () => {
    const bundle = await classifyBundle([
      ...store("s/skin.zarr"),
      ...molecules("s/molecules"),
    ]);

    expect(bundle.cells).toEqual({ kind: "spatialdata", root: "skin.zarr" });
    expect(bundle.molecules).toBe("molecules");
    // images/ belongs to the store; it is not a standalone image
    expect(bundle.image).toBeUndefined();
  });

  it("takes the dropped folder itself as a single part", async () => {
    expect((await classifyBundle(image("dapi.zarr"))).image).toBe("");
    expect((await classifyBundle(store("skin.zarr"))).cells).toEqual({
      kind: "spatialdata",
      root: "",
    });
    expect((await classifyBundle(molecules("m"))).molecules).toBe("");
    expect((await classifyBundle(cells("c"))).cells).toEqual({
      kind: "chunked",
      root: "",
    });
  });

  it("accepts any subset of the parts", async () => {
    const bundle = await classifyBundle([
      ...molecules("s/m"),
      ...image("s/i.zarr"),
    ]);

    expect(bundle.cells).toBeUndefined();
    expect(hasParts(bundle)).toBe(true);
  });

  it("uses the first of several candidates and reports the rest", async () => {
    const bundle = await classifyBundle([
      ...cells("s/cells"),
      ...image("s/b.zarr"),
      ...image("s/a.zarr"),
      ...molecules("s/m2"),
      ...molecules("s/m1"),
    ]);

    expect(bundle.image).toBe("a.zarr");
    expect(bundle.molecules).toBe("m1");
    expect(bundle.skipped).toEqual(["image b.zarr", "molecules m2"]);
  });

  it("does not take a zarr the viewer cannot draw for an image", async () => {
    // An AnnData zarr, and an image with a z axis
    const bundle = await classifyBundle([
      picked("s/anndata.zarr/zarr.json", JSON.stringify({ attributes: {} })),
      picked("s/anndata.zarr/obs/zarr.json", "{}"),
      ...image("s/stack.zarr", ["z", "y", "x"]),
    ]);

    expect(hasParts(bundle)).toBe(false);
  });

  it("finds nothing in raw data, and ignores junk files", async () => {
    const bundle = await classifyBundle([
      picked("s/cells.h5ad"),
      picked("s/transcripts.parquet"),
      picked("s/.DS_Store"),
      picked("s/__MACOSX/x"),
    ]);

    expect(hasParts(bundle)).toBe(false);
    expect(Array.from(bundle.files.keys()).sort()).toEqual([
      "cells.h5ad",
      "transcripts.parquet",
    ]);
  });

  it("hands a part over as a folder of its own", async () => {
    const bundle = await classifyBundle([
      ...cells("s/one"),
      ...image("s/three"),
    ]);

    // The single-folder loaders strip the first path segment
    expect(
      partFiles(bundle, "one")
        .map((f) => f.webkitRelativePath)
        .sort(),
    ).toEqual([
      "one/coords/spatial.bin.gz",
      "one/expr/chunk_00000.bin.gz",
      "one/expr/index.json",
      "one/manifest.json",
    ]);
    expect(Array.from(partFileMap(bundle, "three").keys()).sort()).toEqual([
      "0/c/0/0/0",
      "0/zarr.json",
      "zarr.json",
    ]);
  });
});

describe("fileStore", () => {
  const files = new Map([["a/b", new File(["0123456789"], "b")]]);
  const text = async (bytes?: Uint8Array) =>
    bytes && new TextDecoder().decode(bytes);

  it("reads a whole file by key, with or without a leading slash", async () => {
    expect(await text(await fileStore(files).get("/a/b"))).toBe("0123456789");
    expect(await text(await fileStore(files).get("a/b"))).toBe("0123456789");
  });

  it("answers a missing key with undefined", async () => {
    expect(await fileStore(files).get("/nope")).toBeUndefined();
    expect(
      await fileStore(files).getRange("/nope", { suffixLength: 2 }),
    ).toBeUndefined();
  });

  it("reads only the bytes asked for", async () => {
    // What a sharded array needs: a slice, and the index at the end
    const store = fileStore(files);

    expect(
      await text(await store.getRange("/a/b", { offset: 2, length: 3 })),
    ).toBe("234");
    expect(await text(await store.getRange("/a/b", { suffixLength: 4 }))).toBe(
      "6789",
    );
  });
});

describe("openSpatialDataStoreFromFiles", () => {
  it("roots the store at its folder and looks beside it for sidecars", async () => {
    const folder = new Map([
      ["skin.zarr/zarr.json", new File(["root"], "zarr.json")],
      ["skin.zarr/tables/table/zarr.json", new File(["t"], "zarr.json")],
      ["annotations/annotations.json", new File(["ann"], "annotations.json")],
    ]);
    const opened = openSpatialDataStoreFromFiles(folder, "skin.zarr");

    expect(opened.keys.sort()).toEqual(["tables/table/zarr.json", "zarr.json"]);
    expect(new TextDecoder().decode(await opened.store.get("/zarr.json"))).toBe(
      "root",
    );
    expect(
      new TextDecoder().decode(
        await opened.sidecar.get("/annotations/annotations.json"),
      ),
    ).toBe("ann");
  });
});
