"use client";

import { useEffect, useRef } from "react";
import { createPortal } from "react-dom";

import {
  COLORMAPS,
  DIVERGING_COLORMAPS,
  SEQUENTIAL_COLORMAPS,
  colormapCss,
  type ColormapName,
} from "@/lib/utils/colormaps";

interface ColormapPickerProps {
  x: number;
  y: number;
  current: string;
  onSelect: (name: ColormapName) => void;
  onClose: () => void;
}

const PICKER_WIDTH = 184;

export function ColormapPicker({
  x,
  y,
  current,
  onSelect,
  onClose,
}: ColormapPickerProps) {
  const ref = useRef<HTMLDivElement>(null);

  useEffect(() => {
    const handleDown = (e: MouseEvent) => {
      if (ref.current && !ref.current.contains(e.target as Node)) onClose();
    };
    const handleKey = (e: KeyboardEvent) => {
      if (e.key === "Escape") onClose();
    };
    document.addEventListener("mousedown", handleDown);
    document.addEventListener("keydown", handleKey);
    return () => {
      document.removeEventListener("mousedown", handleDown);
      document.removeEventListener("keydown", handleKey);
    };
  }, [onClose]);

  const items: ColormapName[] = [...DIVERGING_COLORMAPS, ...SEQUENTIAL_COLORMAPS];

  // Clamp to viewport
  const left = Math.min(Math.max(8, x), window.innerWidth - PICKER_WIDTH - 8);
  const approxHeight = 32 + items.length * 28;
  const top = Math.min(Math.max(8, y), window.innerHeight - approxHeight - 8);

  return createPortal(
    <div
      ref={ref}
      className="fixed z-[100] flex flex-col gap-1 p-2 rounded-lg shadow-xl backdrop-blur-md bg-background/85 border border-default-300/40 max-h-[80vh] overflow-y-auto"
      style={{ left, top, width: PICKER_WIDTH }}
    >
      <div className="text-[10px] uppercase tracking-wider text-default-500 px-1 pb-0.5">
        Colormap
      </div>
      {items.map((name) => {
        const def = COLORMAPS[name];
        const isCurrent = name === current;
        return (
          <button
            key={name}
            className={`flex items-center gap-2 px-1.5 py-1 rounded-md transition-colors text-left
              ${isCurrent ? "bg-primary/20 ring-1 ring-primary/60" : "hover:bg-white/10"}`}
            type="button"
            onClick={() => onSelect(name)}
          >
            <span
              className="block h-3 w-20 rounded-sm shadow-inner border border-black/30"
              style={{ background: colormapCss(def.name, "to right") }}
            />
            <span className="text-xs text-default-800">{def.label}</span>
          </button>
        );
      })}
    </div>,
    document.body,
  );
}
