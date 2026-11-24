import { mkdtempSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join } from "node:path";

import { expect, it } from "vitest";

import { jobById, readNpy2 } from "../src/local/analysis.ts";

/** A float32 C-order .npy as numpy.save writes it (version 1.0). */
function npy(rows: number[][]): Buffer {
  const header = `{'descr': '<f4', 'fortran_order': False, 'shape': (${rows.length}, 2), }`;
  const pad = 64 - ((10 + header.length + 1) % 64);
  const h = Buffer.from(`${header}${" ".repeat(pad)}\n`, "latin1");
  const pre = Buffer.alloc(10);

  pre.write("\x93NUMPY", 0, "latin1");
  pre[6] = 1;
  pre[7] = 0;
  pre.writeUInt16LE(h.length, 8);

  return Buffer.concat([
    pre,
    h,
    Buffer.from(new Float32Array(rows.flat()).buffer),
  ]);
}

it("reads an embedding numpy saved, a cell it does not cover as null", () => {
  const f = join(mkdtempSync(join(tmpdir(), "meyes-npy-")), "e.npy");

  writeFileSync(
    f,
    npy([
      [1, 2],
      [NaN, NaN],
      [3.5, -4],
    ]),
  );
  expect(readNpy2(f)).toEqual([[1, 2], null, [3.5, -4]]);
});

it("refuses an array in another layout", () => {
  const f = join(mkdtempSync(join(tmpdir(), "meyes-npy-")), "e.npy");
  const b = npy([[1, 2]]);

  b.write("<f8", b.indexOf("<f4"), "latin1");
  writeFileSync(f, b);
  expect(() => readNpy2(f)).toThrow(/Unexpected array layout/);
});

it("says so for a job it does not know", () => {
  expect(() => jobById("recluster-nope")).toThrow(/No job "recluster-nope"/);
});
