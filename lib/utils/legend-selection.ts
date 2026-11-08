/**
 * What a click on a legend chip does to the set picked for recolouring.
 *
 * Pure, so it is tested on its own; `lib/hooks/useLegendSelection.ts` holds
 * the state and calls this. The gestures the chips already use are left
 * alone — a plain click still opens that chip's own colour picker — so only
 * shift and ⌘/ctrl pick several.
 */

export interface ClickModifiers {
  shiftKey: boolean;
  metaKey: boolean;
  ctrlKey: boolean;
}

export interface SelectionStep {
  selected: Set<string>;
  /** Where the next shift range starts. */
  anchor: string;
  /**
   * True when the click was a selection gesture: the caller does nothing
   * else with it. False for a plain click, which it handles as before.
   */
  handled: boolean;
}

export function nextSelection(
  /** The chips in the order they are drawn, which is what a range means. */
  keys: string[],
  selected: Set<string>,
  anchor: string | null,
  key: string,
  mods: ClickModifiers,
): SelectionStep {
  if (mods.shiftKey) {
    const from = anchor ?? key;
    const a = keys.indexOf(from);
    const b = keys.indexOf(key);

    // A chip that has since gone: treat it as a plain click rather than
    // guessing at a range
    if (a === -1 || b === -1) {
      return { selected: new Set(), anchor: key, handled: false };
    }

    return {
      selected: new Set(keys.slice(Math.min(a, b), Math.max(a, b) + 1)),
      anchor: from,
      handled: true,
    };
  }
  if (mods.metaKey || mods.ctrlKey) {
    const next = new Set(selected);

    if (next.has(key)) next.delete(key);
    else next.add(key);

    return { selected: next, anchor: key, handled: true };
  }

  return { selected: new Set(), anchor: key, handled: false };
}

/** The picked chips still on screen, so a stale key can't be recoloured. */
export function prune(selected: Set<string>, keys: string[]): Set<string> {
  const live = new Set(keys);

  return new Set([...selected].filter((k) => live.has(k)));
}
