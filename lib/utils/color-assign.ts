/**
 * Giving one colour to several items, and distinct colours back again.
 *
 * Shared by the two legends (cell celltypes, single-molecule genes), which
 * differ only in where the colour is written — a cluster's palette or a
 * gene's visualization entry.
 */
import { DEFAULT_COLOR_PALETTE } from "./color-palette";

const PALETTE: readonly string[] = DEFAULT_COLOR_PALETTE;

/** `#RRGGBB` lowercased, for comparing colours written in either case. */
function sameColorKey(color: unknown): string | null {
  if (typeof color !== "string") return null;
  const hex = color.trim().slice(0, 7).toLowerCase();

  return /^#[0-9a-f]{6}$/.test(hex) ? hex : null;
}

/**
 * Distinct colours for `keys`, preferring palette entries that nothing else
 * is already using, so a group that shared one colour becomes
 * distinguishable without disturbing the items around it.
 *
 * `taken` is every colour in use, including the selected items' own; theirs
 * are freed first, since they are about to be replaced. Once the palette is
 * exhausted it is reused in order — 40 colours, so this only bites on a
 * selection larger than any legend shows at once.
 */
export function spreadColors(
  keys: string[],
  taken: Iterable<string> = [],
): Record<string, string> {
  const selected = new Set(keys);
  const used = new Set(
    [...taken].map(sameColorKey).filter((c): c is string => c !== null),
  );

  // A selected item's current colour is not a conflict: it is being replaced
  const free = PALETTE.filter((c) => !used.has(sameColorKey(c)!));
  const pool = free.length >= selected.size ? free : PALETTE;
  const out: Record<string, string> = {};

  keys.forEach((key, i) => {
    out[key] = pool[i % pool.length];
  });

  return out;
}

/**
 * Distinct colours for `keys` in a random order, for when a spread lands on
 * two colours that look alike. Draws from the whole palette (the point is a
 * different arrangement, not avoiding neighbours) and never repeats within
 * the selection while the palette lasts.
 */
export function shuffleColors(keys: string[]): Record<string, string> {
  const pool = [...PALETTE];

  // Fisher-Yates, so every arrangement is equally likely
  for (let i = pool.length - 1; i > 0; i--) {
    const j = Math.floor(Math.random() * (i + 1));

    [pool[i], pool[j]] = [pool[j], pool[i]];
  }

  const out: Record<string, string> = {};

  keys.forEach((key, i) => {
    out[key] = pool[i % pool.length];
  });

  return out;
}
