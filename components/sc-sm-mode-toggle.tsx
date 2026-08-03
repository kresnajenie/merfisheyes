"use client";

import { Hexagon, Atom } from "lucide-react";
import { Tooltip } from "@heroui/tooltip";

import { glassButton } from "@/components/primitives";

export type ScSmMode = "sc" | "sm";

interface ScSmModeToggleProps {
  mode: ScSmMode;
  onChange: (next: ScSmMode) => void;
  // "compact" = horizontal w-14, icons only.
  // "wide"    = horizontal w-32, icons + SC/SM labels (panel header).
  // "vertical"= vertical w-7 h-32, icons stacked (next to rail slider).
  variant?: "compact" | "wide" | "vertical";
}

// Segmented SC/SM mode toggle. Colors mirror the homepage:
// SC = blue-500, SM = purple-500. Icons: Hexagon (cell) / Atom (molecule).
export function ScSmModeToggle({
  mode,
  onChange,
  variant = "compact",
}: ScSmModeToggleProps) {
  const isVertical = variant === "vertical";
  const containerSize = (() => {
    if (variant === "wide") return "w-32 h-7";
    if (variant === "vertical") return "w-7 h-32";

    return "w-14 h-7";
  })();
  const flexDir = isVertical ? "flex-col" : "flex-row";

  return (
    <div
      aria-label="SC/SM mode"
      className={`relative flex ${flexDir} items-center ${containerSize} rounded-full border-2 border-default-200 overflow-hidden ${glassButton()}`}
      role="radiogroup"
    >
      <Tooltip content="Single cell" placement={isVertical ? "left" : "top"}>
        <button
          aria-checked={mode === "sc"}
          aria-label="Single cell"
          className={`flex-1 w-full h-full flex items-center justify-center gap-1 transition-colors ${
            mode === "sc"
              ? "bg-blue-500 text-white shadow-[0_4px_12px_rgba(59,130,246,0.45)]"
              : "text-blue-300/70 hover:bg-blue-500/15 hover:text-blue-200"
          }`}
          role="radio"
          type="button"
          onClick={() => onChange("sc")}
        >
          <Hexagon className="w-3.5 h-3.5" />
          {variant === "wide" && (
            <span className="text-[10px] font-semibold tracking-wider uppercase">
              SC
            </span>
          )}
        </button>
      </Tooltip>
      <Tooltip content="Single molecule" placement={isVertical ? "left" : "top"}>
        <button
          aria-checked={mode === "sm"}
          aria-label="Single molecule"
          className={`flex-1 w-full h-full flex items-center justify-center gap-1 transition-colors ${
            mode === "sm"
              ? "bg-purple-500 text-white shadow-[0_4px_12px_rgba(168,85,247,0.45)]"
              : "text-purple-300/70 hover:bg-purple-500/15 hover:text-purple-200"
          }`}
          role="radio"
          type="button"
          onClick={() => onChange("sm")}
        >
          <Atom className="w-3.5 h-3.5" />
          {variant === "wide" && (
            <span className="text-[10px] font-semibold tracking-wider uppercase">
              SM
            </span>
          )}
        </button>
      </Tooltip>
    </div>
  );
}
