import { describe, it, expect } from "vitest";
import { readFileSync } from "fs";
import { gunzipSync } from "zlib";
import path from "path";

import {
  decodeObsColumn,
  encodeCategoricalObs,
  encodeNumericalObs,
} from "@/lib/utils/obs-binary";
import { encodeDeStatsBuffer, parseDeStatsBuffer } from "@/lib/utils/de-stats";
import { GeneChunkProcessor } from "@/lib/utils/GeneChunkProcessor";

const FIXTURES = path.join(__dirname, "..", "fixtures", "obs-binary");

function toArrayBuffer(u8: Uint8Array): ArrayBuffer {
  return u8.buffer.slice(u8.byteOffset, u8.byteOffset + u8.byteLength) as ArrayBuffer;
}

describe("obs binary format", () => {
  it("round-trips a categorical column (uint8 codes, padded dictionary)", () => {
    const dict = ["", "Astro", "Oligo", "微"]; // non-ASCII + "" missing label
    const codes = [1, 2, 2, 0, 3, 1];
    const decoded = decodeObsColumn(toArrayBuffer(encodeCategoricalObs(dict, codes)));

    expect(decoded.kind).toBe("categorical");
    if (decoded.kind !== "categorical") return;
    expect(decoded.dict).toEqual(dict);
    expect(decoded.codes).toBeInstanceOf(Uint8Array);
    expect(Array.from(decoded.codes)).toEqual(codes);
  });

  it("widens codes to uint16/uint32 with the dictionary size", () => {
    const dict16 = Array.from({ length: 300 }, (_, i) => `c${i}`);
    const d16 = decodeObsColumn(toArrayBuffer(encodeCategoricalObs(dict16, [0, 299])));

    expect(d16.kind === "categorical" && d16.codes).toBeInstanceOf(Uint16Array);

    const dict32 = Array.from({ length: 70000 }, (_, i) => `c${i}`);
    const d32 = decodeObsColumn(toArrayBuffer(encodeCategoricalObs(dict32, [69999])));

    expect(d32.kind === "categorical" && d32.codes).toBeInstanceOf(Uint32Array);
    expect(d32.kind === "categorical" && d32.codes[0]).toBe(69999);
  });

  it("round-trips a numerical column with NaN for missing", () => {
    const decoded = decodeObsColumn(
      toArrayBuffer(encodeNumericalObs([1.5, NaN, -2, 0])),
    );

    expect(decoded.kind).toBe("numerical");
    if (decoded.kind !== "numerical") return;
    expect(decoded.values[0]).toBe(1.5);
    expect(Number.isNaN(decoded.values[1])).toBe(true);
    expect(decoded.values[2]).toBe(-2);
  });

  // Cross-language: files written by scripts/process_spatial_data.py
  // (encode_obs_binary) must decode with the TS reader.
  it("decodes a categorical column written by the Python pipeline", () => {
    const expected = JSON.parse(readFileSync(path.join(FIXTURES, "expected.json"), "utf-8"));
    const buf = gunzipSync(readFileSync(path.join(FIXTURES, "py_categorical_leiden.bin.gz")));
    const decoded = decodeObsColumn(toArrayBuffer(new Uint8Array(buf)));

    expect(decoded.kind).toBe("categorical");
    if (decoded.kind !== "categorical") return;
    // Dictionary is in sorted label order
    expect(decoded.dict).toEqual([...decoded.dict].sort());
    const labels = Array.from(decoded.codes, (c) => decoded.dict[c]);

    expect(labels).toEqual(expected.categorical.labels);
  });

  it("decodes a numerical column written by the Python pipeline", () => {
    const expected = JSON.parse(readFileSync(path.join(FIXTURES, "expected.json"), "utf-8"));
    const buf = gunzipSync(readFileSync(path.join(FIXTURES, "py_numerical.bin.gz")));
    const decoded = decodeObsColumn(toArrayBuffer(new Uint8Array(buf)));

    expect(decoded.kind).toBe("numerical");
    if (decoded.kind !== "numerical") return;
    const values = Array.from(decoded.values, (v) => (Number.isNaN(v) ? null : v));

    expect(values).toEqual(expected.numerical.values);
  });
});

describe("DE stats binary", () => {
  it("encode → parse round-trips (same layout as write_de_stats_binary)", () => {
    const stats = {
      column: "leiden",
      celltypes: ["A", "B", "Cé"],
      cellCounts: [10, 0, 5],
      genes: ["g1", "g2"],
      means: new Float32Array([1, 2, 3, 4, 5, 6]),
      pctExpressing: new Float32Array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6]),
    };
    const parsed = parseDeStatsBuffer(
      toArrayBuffer(encodeDeStatsBuffer(stats)),
      "leiden",
      stats.genes,
    );

    expect(parsed.celltypes).toEqual(stats.celltypes);
    expect(parsed.cellCounts).toEqual(stats.cellCounts);
    expect(Array.from(parsed.means)).toEqual(Array.from(stats.means));
    expect(Array.from(parsed.pctExpressing).map((v) => Number(v.toFixed(6)))).toEqual(
      Array.from(stats.pctExpressing).map((v) => Number(v.toFixed(6))),
    );
  });
});

describe("expression chunking", () => {
  it("defaults to one gene per chunk (shared spec with the server pipeline)", () => {
    const p = new GeneChunkProcessor();

    expect(p.determineChunkSize(3)).toBe(1);
    expect(p.determineChunkSize(20000)).toBe(1);
    expect(new GeneChunkProcessor(50).determineChunkSize(20000)).toBe(50);
  });
});
