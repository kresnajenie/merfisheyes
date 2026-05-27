// Resolve the effective order of a secondary column's values.
//
// `naturalOrder` is the column's `uniqueValues` (or a deduped fallback).
// `customOrder` is the user's reorder; empty means "no customisation".
//
// Returns values from `customOrder` first (in that order, filtered to ones
// still present in `naturalOrder`), then any `naturalOrder` values missing
// from `customOrder` appended at the end. The "append missing" branch
// defends against columns gaining values the user hasn't touched yet.
export function getOrderedSecondaryValues(
  naturalOrder: string[],
  customOrder: string[],
): string[] {
  if (customOrder.length === 0) return naturalOrder;
  const natural = new Set(naturalOrder);
  const seen = new Set<string>();
  const out: string[] = [];
  for (const v of customOrder) {
    if (natural.has(v) && !seen.has(v)) {
      out.push(v);
      seen.add(v);
    }
  }
  for (const v of naturalOrder) {
    if (!seen.has(v)) out.push(v);
  }
  return out;
}
