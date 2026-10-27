"use client";

import type { RefObject } from "react";

import { useEffect, useRef } from "react";

/**
 * Call `onOutside` when the pointer goes down outside every element of
 * `refs` — a flyout panel and the button that opens it — while `enabled`.
 * A toast is never "outside" (closing an error must not close the panel the
 * user is working in).
 */
export function useClickOutside(
  refs: RefObject<HTMLElement | null>[],
  onOutside: () => void,
  enabled = true,
): void {
  const latest = useRef(onOutside);

  latest.current = onOutside;

  useEffect(() => {
    if (!enabled) return;
    const handle = (event: PointerEvent) => {
      const target = event.target as Element | null;

      if (!target) return;
      if (target.closest?.("[data-react-aria-top-layer]")) return;
      if (refs.some((r) => r.current?.contains(target))) return;
      latest.current();
    };

    document.addEventListener("pointerdown", handle);

    return () => document.removeEventListener("pointerdown", handle);
    // The refs are stable objects; their contents are read at event time
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [enabled]);
}
