"use client";

import { useEffect, useRef, useState } from "react";
import { createPortal } from "react-dom";

import { usePanelVisualizationStore } from "@/lib/hooks/usePanelStores";

type RangeOverride = { min: number; max: number };

/**
 * Hook that supplies right-click range editing for size sliders.
 *
 * - Effective range: explicit override if set (literal), otherwise the
 *   config default expanded to fit the current value (so URL-restored
 *   out-of-default values stay reachable).
 * - Wrap the slider in a div spreading `wrapperProps`; render `popover`
 *   alongside (it portals to body itself).
 * - Validation: setSliderRange in the store rejects min >= max.
 */
export function useSliderRange(
  rangeKey: string,
  defaultMin: number,
  defaultMax: number,
  currentValue: number,
) {
  const override = usePanelVisualizationStore(
    (s) => s.sliderRanges[rangeKey] as RangeOverride | undefined,
  );
  const setSliderRange = usePanelVisualizationStore((s) => s.setSliderRange);

  const min = override ? override.min : Math.min(defaultMin, currentValue);
  const max = override ? override.max : Math.max(defaultMax, currentValue);

  const [popoverPos, setPopoverPos] = useState<{ x: number; y: number } | null>(
    null,
  );

  const onContextMenu = (e: React.MouseEvent) => {
    e.preventDefault();
    e.stopPropagation();
    setPopoverPos({ x: e.clientX, y: e.clientY });
  };

  const popover =
    popoverPos !== null ? (
      <RangePopover
        max={max}
        min={min}
        x={popoverPos.x}
        y={popoverPos.y}
        onChange={(newMin, newMax) =>
          setSliderRange(rangeKey, newMin, newMax)
        }
        onClose={() => setPopoverPos(null)}
      />
    ) : null;

  return { min, max, onContextMenu, popover };
}

interface RangePopoverProps {
  x: number;
  y: number;
  min: number;
  max: number;
  onChange: (min: number, max: number) => void;
  onClose: () => void;
}

function RangePopover({ x, y, min, max, onChange, onClose }: RangePopoverProps) {
  const ref = useRef<HTMLDivElement>(null);
  const [minStr, setMinStr] = useState(String(min));
  const [maxStr, setMaxStr] = useState(String(max));

  useEffect(() => setMinStr(String(min)), [min]);
  useEffect(() => setMaxStr(String(max)), [max]);

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

  const tryApply = (rawMin: string, rawMax: string) => {
    const m = parseFloat(rawMin);
    const M = parseFloat(rawMax);
    if (!Number.isFinite(m) || !Number.isFinite(M)) return;
    if (!(m < M)) return;
    onChange(m, M);
  };

  const inputCls =
    "w-14 px-1 py-0.5 text-xs rounded bg-default-100/60 border border-default-300/40 text-default-900 outline-none focus:border-primary";

  // Clamp to viewport
  const left = Math.min(Math.max(8, x), window.innerWidth - 160);
  const top = Math.min(Math.max(8, y), window.innerHeight - 60);
  const invalid =
    !(parseFloat(minStr) < parseFloat(maxStr)) ||
    !Number.isFinite(parseFloat(minStr)) ||
    !Number.isFinite(parseFloat(maxStr));

  return createPortal(
    <div
      ref={ref}
      className="fixed z-[100] flex items-center gap-1 px-2 py-1.5 rounded-lg shadow-lg backdrop-blur-md bg-background/80 border border-default-300/40"
      style={{ left, top }}
    >
      <input
        aria-label="min"
        className={`${inputCls} ${invalid ? "border-danger/60" : ""}`}
        type="number"
        value={minStr}
        onChange={(e) => {
          setMinStr(e.target.value);
          tryApply(e.target.value, maxStr);
        }}
      />
      <span className="text-[10px] text-default-500">→</span>
      <input
        aria-label="max"
        className={`${inputCls} ${invalid ? "border-danger/60" : ""}`}
        type="number"
        value={maxStr}
        onChange={(e) => {
          setMaxStr(e.target.value);
          tryApply(minStr, e.target.value);
        }}
      />
    </div>,
    document.body,
  );
}
