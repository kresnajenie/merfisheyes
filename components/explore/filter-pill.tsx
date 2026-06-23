"use client";

import { useEffect, useMemo, useRef, useState } from "react";
import { createPortal } from "react-dom";

interface FilterPillProps {
  /** Label shown when no value is selected, e.g. "Species". */
  label: string;
  /** Available values to choose from. */
  options: string[];
  /** Currently-selected value, or empty string for none. */
  value: string;
  /** Called with "" to clear, or with a value to commit. */
  onChange: (next: string) => void;
  /** Optional placeholder for the in-popover search input. Defaults to "Search". */
  searchPlaceholder?: string;
}

const POPOVER_WIDTH = 240;
const POPOVER_MAX_HEIGHT = 320;

/**
 * LinkedIn-style filter chip.
 *
 *   - inactive: outlined pill with the category label + chevron
 *   - active:   filled pill showing "Label: Value" + an inline × to clear
 *
 * Clicking the body opens a portal-mounted popover with a search box and a
 * scrollable list of options. Picking an option commits immediately
 * (single-select) and closes. Outside-click and Esc close without committing.
 */
export function FilterPill({
  label,
  options,
  value,
  onChange,
  searchPlaceholder,
}: FilterPillProps) {
  const buttonRef = useRef<HTMLButtonElement>(null);
  const popoverRef = useRef<HTMLDivElement>(null);
  const searchRef = useRef<HTMLInputElement>(null);
  const [isOpen, setIsOpen] = useState(false);
  const [query, setQuery] = useState("");
  const [pos, setPos] = useState<{ top: number; left: number } | null>(null);

  const isActive = value.length > 0;

  // Compute popover position from the pill's rect; clamp to viewport.
  const reposition = () => {
    const rect = buttonRef.current?.getBoundingClientRect();

    if (!rect) return;
    const left = Math.min(
      Math.max(8, rect.left),
      window.innerWidth - POPOVER_WIDTH - 8,
    );
    const top = Math.min(rect.bottom + 6, window.innerHeight - 8);

    setPos({ top, left });
  };

  // Open / close lifecycle: position, outside-click, Esc, scroll reposition.
  useEffect(() => {
    if (!isOpen) return;
    reposition();

    const onDown = (e: MouseEvent) => {
      const target = e.target as Node;

      if (buttonRef.current?.contains(target)) return;
      if (popoverRef.current?.contains(target)) return;
      setIsOpen(false);
    };
    const onKey = (e: KeyboardEvent) => {
      if (e.key === "Escape") setIsOpen(false);
    };
    const onScroll = () => reposition();

    document.addEventListener("mousedown", onDown);
    document.addEventListener("keydown", onKey);
    window.addEventListener("scroll", onScroll, true);
    window.addEventListener("resize", onScroll);

    return () => {
      document.removeEventListener("mousedown", onDown);
      document.removeEventListener("keydown", onKey);
      window.removeEventListener("scroll", onScroll, true);
      window.removeEventListener("resize", onScroll);
    };
  }, [isOpen]);

  // Reset the in-popover search every time we open + focus it for keyboard-first
  // typing without the jsx-a11y/no-autofocus prop warning.
  useEffect(() => {
    if (!isOpen) return;
    setQuery("");
    searchRef.current?.focus();
  }, [isOpen]);

  const filtered = useMemo(() => {
    if (!query.trim()) return options;
    const q = query.toLowerCase();

    return options.filter((o) => o.toLowerCase().includes(q));
  }, [options, query]);

  const pillClass = isActive
    ? "bg-foreground text-background border-foreground"
    : "bg-transparent text-default-700 border-default-300 hover:bg-default-100";

  return (
    <>
      <button
        ref={buttonRef}
        className={`inline-flex items-center gap-1.5 h-8 px-3 rounded-full border text-sm transition-colors ${pillClass}`}
        type="button"
        onClick={() => setIsOpen((p) => !p)}
      >
        <span className="whitespace-nowrap">
          {isActive ? (
            <>
              <span className="opacity-70">{label}:</span> {value}
            </>
          ) : (
            label
          )}
        </span>
        {isActive ? (
          <span
            aria-label={`Clear ${label}`}
            className="inline-flex items-center justify-center w-4 h-4 rounded-full hover:bg-background/20"
            role="button"
            tabIndex={0}
            onClick={(e) => {
              e.stopPropagation();
              onChange("");
              setIsOpen(false);
            }}
            onKeyDown={(e) => {
              if (e.key === "Enter" || e.key === " ") {
                e.preventDefault();
                e.stopPropagation();
                onChange("");
                setIsOpen(false);
              }
            }}
          >
            <svg
              className="w-3 h-3"
              fill="none"
              stroke="currentColor"
              strokeWidth={2.5}
              viewBox="0 0 24 24"
            >
              <path
                d="M6 18L18 6M6 6l12 12"
                strokeLinecap="round"
                strokeLinejoin="round"
              />
            </svg>
          </span>
        ) : (
          <svg
            className="w-3.5 h-3.5 opacity-70"
            fill="none"
            stroke="currentColor"
            strokeWidth={2}
            viewBox="0 0 24 24"
          >
            <path
              d="M6 9l6 6 6-6"
              strokeLinecap="round"
              strokeLinejoin="round"
            />
          </svg>
        )}
      </button>

      {isOpen && pos && typeof document !== "undefined"
        ? createPortal(
            <div
              ref={popoverRef}
              className="fixed z-[200] flex flex-col rounded-xl shadow-lg border border-default-200 bg-content1 overflow-hidden"
              style={{
                top: pos.top,
                left: pos.left,
                width: POPOVER_WIDTH,
                maxHeight: POPOVER_MAX_HEIGHT,
              }}
            >
              <div className="p-2 border-b border-default-200">
                <input
                  ref={searchRef}
                  className="w-full bg-default-100 rounded-md px-2 py-1.5 text-sm outline-none focus:ring-1 focus:ring-primary/50"
                  placeholder={
                    searchPlaceholder ?? `Search ${label.toLowerCase()}`
                  }
                  type="text"
                  value={query}
                  onChange={(e) => setQuery(e.target.value)}
                />
              </div>
              <div className="overflow-y-auto py-1">
                {filtered.length === 0 ? (
                  <p className="text-xs text-default-400 text-center py-4">
                    No matches
                  </p>
                ) : (
                  filtered.map((opt) => {
                    const isSelected = opt === value;

                    return (
                      <button
                        key={opt}
                        className={`w-full text-left text-sm px-3 py-1.5 hover:bg-default-100 transition-colors ${
                          isSelected ? "bg-primary/10 text-primary" : ""
                        }`}
                        type="button"
                        onClick={() => {
                          onChange(opt);
                          setIsOpen(false);
                        }}
                      >
                        {opt}
                      </button>
                    );
                  })
                )}
              </div>
            </div>,
            document.body,
          )
        : null}
    </>
  );
}
