/**
 * Binary obs-column format (`obs/{col}.bin.gz`, version 1).
 *
 * Replaces the per-cell JSON arrays (`["Astro","Oligo",…]`), which cost tens
 * of MB and a JSON.parse that materialises one string per cell. Here a
 * categorical column is a dictionary + one small integer code per cell —
 * exactly the representation the viewer keeps in memory — and a numerical
 * column is a float32 per cell. Decoding is a typed-array view, no parse.
 *
 * Layout (little-endian), identical in scripts/process_spatial_data.py:
 *
 *   u32  version        = 1
 *   u32  numCells
 *   u8   kind           0 = categorical, 1 = numerical
 *   u8   width          bytes per cell: 1 | 2 | 4 (codes) — always 4 for f32
 *   u16  reserved       0
 *   u32  dictBytes      length of the UTF-8 JSON array of category labels
 *                       (0 for numerical)
 *   …    dict bytes, zero-padded to a 4-byte boundary
 *   …    numCells × width: uint8/16/32 codes into dict, or float32 values
 *
 * Categorical codes are indices into `dict`, which both writers emit in
 * sorted label order (the same order used for palette assignment). Missing
 * categorical values are the label "" (a normal dictionary entry); missing
 * numerical values are NaN.
 */

export const OBS_BINARY_FORMAT = "bin-v1";
export const OBS_BINARY_VERSION = 1;
const HEADER_BYTES = 16;

export type ObsColumnDecoded =
  | { kind: "categorical"; dict: string[]; codes: Uint8Array | Uint16Array | Uint32Array }
  | { kind: "numerical"; values: Float32Array };

function codeWidth(numCategories: number): 1 | 2 | 4 {
  if (numCategories <= 0xff) return 1;
  if (numCategories <= 0xffff) return 2;

  return 4;
}

function pad4(n: number): number {
  return (n + 3) & ~3;
}

/** Encode a categorical column from its dictionary + per-cell codes. */
export function encodeCategoricalObs(
  dict: string[],
  codes: ArrayLike<number>,
): Uint8Array {
  const width = codeWidth(dict.length);
  const dictBytes = new TextEncoder().encode(JSON.stringify(dict));
  const dataOffset = HEADER_BYTES + pad4(dictBytes.length);
  const out = new Uint8Array(dataOffset + codes.length * width);
  const view = new DataView(out.buffer);

  view.setUint32(0, OBS_BINARY_VERSION, true);
  view.setUint32(4, codes.length, true);
  view.setUint8(8, 0);
  view.setUint8(9, width);
  view.setUint16(10, 0, true);
  view.setUint32(12, dictBytes.length, true);
  out.set(dictBytes, HEADER_BYTES);

  const Arr = width === 1 ? Uint8Array : width === 2 ? Uint16Array : Uint32Array;
  const data = new Arr(out.buffer, dataOffset, codes.length);

  for (let i = 0; i < codes.length; i++) data[i] = codes[i];

  return out;
}

/** Encode a numerical column (NaN for missing). */
export function encodeNumericalObs(values: ArrayLike<number>): Uint8Array {
  const out = new Uint8Array(HEADER_BYTES + values.length * 4);
  const view = new DataView(out.buffer);

  view.setUint32(0, OBS_BINARY_VERSION, true);
  view.setUint32(4, values.length, true);
  view.setUint8(8, 1);
  view.setUint8(9, 4);
  view.setUint16(10, 0, true);
  view.setUint32(12, 0, true);
  new Float32Array(out.buffer, HEADER_BYTES, values.length).set(
    values instanceof Float32Array ? values : Float32Array.from(values),
  );

  return out;
}

/** Decode a (gunzipped) obs column buffer. */
export function decodeObsColumn(buffer: ArrayBuffer): ObsColumnDecoded {
  const view = new DataView(buffer);
  const version = view.getUint32(0, true);

  if (version !== OBS_BINARY_VERSION) {
    throw new Error(`Unsupported obs column format version ${version}`);
  }
  const numCells = view.getUint32(4, true);
  const kind = view.getUint8(8);
  const width = view.getUint8(9);
  const dictBytes = view.getUint32(12, true);
  const dataOffset = HEADER_BYTES + pad4(dictBytes);

  if (kind === 1) {
    return {
      kind: "numerical",
      values: new Float32Array(buffer, dataOffset, numCells),
    };
  }

  const dict: string[] = JSON.parse(
    new TextDecoder().decode(new Uint8Array(buffer, HEADER_BYTES, dictBytes)),
  );
  const Arr = width === 1 ? Uint8Array : width === 2 ? Uint16Array : Uint32Array;

  return { kind: "categorical", dict, codes: new Arr(buffer, dataOffset, numCells) };
}
