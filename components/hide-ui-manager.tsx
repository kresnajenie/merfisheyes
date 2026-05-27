"use client";

import { useEffect } from "react";

import { useSplitScreenStore } from "@/lib/stores/splitScreenStore";

/**
 * Mounts once at the app root. When `hideUi` is true:
 *   - Adds `.ui-hidden` to <body> so the global CSS rule hides every
 *     [data-ui-overlay] element.
 *   - Renders a small floating chip that brings the UI back.
 *   - Listens for the `H` keyboard shortcut to toggle.
 *
 * The chip itself has no `data-ui-overlay` attribute so it stays visible.
 */
export function HideUiManager() {
  const hideUi = useSplitScreenStore((s) => s.hideUi);
  const toggleHideUi = useSplitScreenStore((s) => s.toggleHideUi);
  const setHideUi = useSplitScreenStore((s) => s.setHideUi);

  // Sync body class.
  useEffect(() => {
    if (typeof document === "undefined") return;
    if (hideUi) {
      document.body.classList.add("ui-hidden");
    } else {
      document.body.classList.remove("ui-hidden");
    }
    return () => {
      document.body.classList.remove("ui-hidden");
    };
  }, [hideUi]);

  // Keyboard shortcut: H toggles. Ignored while typing in inputs/textareas
  // or when a modifier key is held.
  useEffect(() => {
    if (typeof window === "undefined") return;
    const handler = (e: KeyboardEvent) => {
      if (e.key !== "h" && e.key !== "H") return;
      if (e.ctrlKey || e.metaKey || e.altKey) return;
      const target = e.target as HTMLElement | null;
      if (target) {
        const tag = target.tagName;
        if (tag === "INPUT" || tag === "TEXTAREA" || target.isContentEditable)
          return;
      }
      e.preventDefault();
      toggleHideUi();
    };
    window.addEventListener("keydown", handler);
    return () => window.removeEventListener("keydown", handler);
  }, [toggleHideUi]);

  if (!hideUi) return null;

  return (
    <>
      <style>{`
        body.ui-hidden [data-ui-overlay],
        body.ui-hidden .Toastify {
          display: none !important;
        }
      `}</style>
      <button
        className="fixed top-4 left-4 z-[1000] flex items-center gap-2 px-3 py-1.5 rounded-full bg-black/60 hover:bg-black/80 backdrop-blur-md border border-white/20 text-white text-xs font-medium shadow-lg transition-colors"
        title="Show UI (H)"
        type="button"
        onClick={() => setHideUi(false)}
      >
      <svg
        className="w-3.5 h-3.5"
        fill="none"
        stroke="currentColor"
        strokeWidth={2}
        viewBox="0 0 24 24"
      >
        <path
          d="M1 12s4-8 11-8 11 8 11 8-4 8-11 8-11-8-11-8z"
          strokeLinecap="round"
          strokeLinejoin="round"
        />
        <circle cx={12} cy={12} r={3} strokeLinecap="round" strokeLinejoin="round" />
      </svg>
      Show UI
      </button>
    </>
  );
}
