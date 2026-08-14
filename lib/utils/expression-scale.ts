/**
 * Raw-counts vs log-normalized expression, and the toggle between them.
 *
 * The stored expression matrix is whatever the uploader provided — either raw
 * integer counts or an already-normalized (usually log1p) matrix. We can't know
 * which from metadata alone, so we sniff it from the values: raw counts are
 * whole numbers, a log/normalized matrix has fractional values. The gene scale
 * bar then lets the user flip the *display* between the two spaces.
 */

export type ExpressionScaleMode = "raw" | "log";
export type DetectedExpressionKind = "counts" | "normalized";

/** Fraction of sampled nonzero values that must be fractional to call it
 * normalized. Raw counts are all integers (→ 0), a log matrix is almost all
 * fractional (→ ~1), so any modest threshold separates them cleanly. */
const FRACTIONAL_THRESHOLD = 0.05;
const MAX_SAMPLES = 5000;

/**
 * Decide whether an expression vector is raw counts or already normalized, by
 * sampling its nonzero values and asking how many are non-integers. Zeros are
 * skipped (they're integers in either space and dominate sparse data).
 */
export function detectExpressionKind(
  values: ArrayLike<number>,
): DetectedExpressionKind {
  const n = values.length;

  if (n === 0) return "counts";

  const step = Math.max(1, Math.floor(n / MAX_SAMPLES));
  let checked = 0;
  let fractional = 0;

  for (let i = 0; i < n; i += step) {
    const v = values[i];

    if (!isFinite(v) || v === 0) continue;
    checked++;
    if (Math.abs(v - Math.round(v)) > 1e-6) fractional++;
  }

  if (checked === 0) return "counts";

  return fractional / checked > FRACTIONAL_THRESHOLD ? "normalized" : "counts";
}

/** The space the data is stored in — what "show it as-is" means. */
export function nativeModeFor(
  kind: DetectedExpressionKind,
): ExpressionScaleMode {
  return kind === "normalized" ? "log" : "raw";
}

/**
 * Map one stored value into the requested display space, given the detected
 * base. Counts ⇄ log1p is exact; going "raw" on an already-log matrix uses
 * expm1 to linearize back toward counts (approximate — not the true integers).
 */
export function toDisplaySpace(
  value: number,
  base: DetectedExpressionKind,
  mode: ExpressionScaleMode,
): number {
  if (!isFinite(value)) return value;

  if (base === "counts") {
    return mode === "log" ? Math.log1p(value) : value;
  }

  // base === "normalized" (log1p space)
  return mode === "raw" ? Math.expm1(value) : value;
}

/** Whether a (base, mode) pair is a no-op — lets callers skip a copy. */
export function isNativeSpace(
  base: DetectedExpressionKind,
  mode: ExpressionScaleMode,
): boolean {
  return nativeModeFor(base) === mode;
}

/**
 * Transform a whole expression vector into the display space. Returns the input
 * untouched when it's already in that space (no allocation).
 */
export function transformExpression(
  values: number[],
  base: DetectedExpressionKind,
  mode: ExpressionScaleMode,
): number[] {
  if (isNativeSpace(base, mode)) return values;

  return values.map((v) => toDisplaySpace(v, base, mode));
}
