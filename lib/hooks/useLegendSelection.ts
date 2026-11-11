"use client";

import { useCallback, useEffect, useRef, useState } from "react";

import {
  type ClickModifiers,
  nextSelection,
  prune,
} from "@/lib/utils/legend-selection";

/**
 * Picking several legend chips to recolour at once, the same way in both
 * legends (cell celltypes, single-molecule genes). What a click means lives
 * in `lib/utils/legend-selection.ts`; this holds the state.
 *
 *   - shift-click: every chip between the last one picked and this one
 *   - ⌘/ctrl-click: add or remove this chip
 *   - a plain click: the chip's own action, and the selection goes away
 *
 * `keys` is the chips in the order they are drawn. Chips that disappear (a
 * celltype deselected, a gene removed) drop out of the selection.
 */
export interface LegendSelection {
  selected: Set<string>;
  /**
   * Handle a click on a chip. Returns true when the click was a selection
   * gesture and the caller should do nothing else; false for a plain click,
   * which the caller handles as before.
   */
  onChipClick: (key: string, e: ClickModifiers) => boolean;
  clear: () => void;
}

export function useLegendSelection(keys: string[]): LegendSelection {
  const [selected, setSelected] = useState<Set<string>>(new Set());
  // Where a shift range starts: the last chip picked on its own
  const anchor = useRef<string | null>(null);
  const keysRef = useRef(keys);
  // `onChipClick` goes to every chip, so it must not be rebuilt whenever the
  // selection changes: it reads the current keys and set through these refs
  const selectedRef = useRef(selected);

  keysRef.current = keys;
  selectedRef.current = selected;

  // Keep only chips still on screen
  const present = keys.join("\u0000");

  useEffect(() => {
    setSelected((prev) => {
      if (prev.size === 0) return prev;
      const next = prune(prev, keysRef.current);

      return next.size === prev.size ? prev : next;
    });
  }, [present]);

  const onChipClick = useCallback((key: string, e: ClickModifiers) => {
    const step = nextSelection(
      keysRef.current,
      selectedRef.current,
      anchor.current,
      key,
      e,
    );

    anchor.current = step.anchor;
    setSelected((prev) =>
      prev.size === 0 && step.selected.size === 0 ? prev : step.selected,
    );

    return step.handled;
  }, []);

  const clear = useCallback(() => {
    setSelected((prev) => (prev.size === 0 ? prev : new Set()));
  }, []);

  return { selected, onChipClick, clear };
}
