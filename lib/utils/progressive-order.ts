/**
 * Progressive molecule order for single-molecule gene files.
 *
 * A gene file is streamed into the viewer as it downloads, so the order the
 * molecules are stored in decides what a partially loaded gene looks like. In
 * source order that is a spatial sweep — one corner of the slide fills in
 * while the rest is empty. Writing the molecules in bit-reversal order instead
 * makes any prefix of the file an evenly spread sample of the whole slide, so
 * the picture starts complete-looking and simply densifies.
 *
 * The permutation is the classic bit-reversal of `0 … 2^k - 1` (k the smallest
 * power of two covering n) with out-of-range values dropped. It is a pure
 * function of n — no RNG — so the Python pipeline
 * (`progressive_order` in scripts/process_single_molecule.py) produces exactly
 * the same order and both pipelines emit byte-identical files.
 */

/** Reverse the low `bits` bits of `value`. */
function reverseBits(value: number, bits: number): number {
  let out = 0;

  for (let i = 0; i < bits; i++) {
    out = (out << 1) | ((value >>> i) & 1);
  }

  return out >>> 0;
}

/**
 * Indices `0 … n-1`, reordered so that every prefix is spread over the whole
 * range. Returns a Uint32Array of length n (identity for n < 2).
 */
export function progressiveOrder(n: number): Uint32Array {
  const out = new Uint32Array(n);

  if (n < 2) {
    for (let i = 0; i < n; i++) out[i] = i;

    return out;
  }

  let bits = 1;

  while (1 << bits < n) bits++;

  const size = 1 << bits;
  let written = 0;

  for (let i = 0; i < size; i++) {
    const value = reverseBits(i, bits);

    if (value < n) out[written++] = value;
  }

  return out;
}

/**
 * Destination slot for each source molecule — the inverse of
 * `progressiveOrder`, i.e. the array `slots` with
 * `slots[progressiveOrder(n)[i]] === i`.
 *
 * The Python pipeline reorders a gene by *gathering*
 * (`out[i] = src[order[i]]`); a browser pass that fills molecules as it reads
 * them has to *scatter* (`out[slots[j]] = src[j]`) to land on the same bytes.
 * The two are only the same permutation when n is a power of two, so this
 * distinction is what keeps the pipelines byte-identical.
 */
export function progressiveSlots(n: number): Uint32Array {
  const order = progressiveOrder(n);
  const slots = new Uint32Array(n);

  for (let i = 0; i < n; i++) slots[order[i]] = i;

  return slots;
}
