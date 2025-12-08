"use client";

import React, { useRef, useCallback, useState } from "react";

import { cn } from "@/lib/utils";

interface NumberScrubberProps {
  value: number;
  onChange?: (value: number) => void;
  min?: number;
  max?: number;
  step?: number;
  decimals?: number;
  className?: string;
  disabled?: boolean;
}

export const NumberScrubber = React.forwardRef<
  HTMLDivElement,
  NumberScrubberProps
>(
  (
    {
      value,
      onChange,
      min = -Infinity,
      max = Infinity,
      step = 0.1,
      decimals = 2,
      className,
      disabled = false,
    },
    ref,
  ) => {
    // Use refs to avoid stale closure issues
    const isDraggingRef = useRef(false);
    const startXRef = useRef(0);
    const startValueRef = useRef(0);

    // Track the live dragging value for smooth display updates
    const [draggingValue, setDraggingValue] = useState<number | null>(null);

    // Debounce timer for onChange callback
    const debounceTimerRef = useRef<NodeJS.Timeout | null>(null);

    const clampValue = useCallback(
      (val: number) => {
        return Math.max(min, Math.min(max, val));
      },
      [min, max],
    );

    const formatValue = useCallback(
      (val: number) => {
        return val.toFixed(decimals);
      },
      [decimals],
    );

    const handleMouseDown = useCallback(
      (e: React.MouseEvent) => {
        if (disabled) return;
        e.preventDefault();
        isDraggingRef.current = true;
        startXRef.current = e.clientY; // Changed from clientX to clientY
        startValueRef.current = value;
        setDraggingValue(value);
        document.body.style.cursor = "ns-resize"; // Changed from ew-resize to ns-resize
      },
      [disabled, value],
    );

    const handleMouseMove = useCallback(
      (e: MouseEvent) => {
        if (!isDraggingRef.current) return;
        const deltaY = e.clientY - startXRef.current; // Changed from clientX to clientY
        const delta = -deltaY * step; // Negative so dragging up increases value
        const newValue = clampValue(startValueRef.current + delta);

        // Update the display immediately (no debounce)
        setDraggingValue(newValue);

        // Debounce the onChange callback (100ms after last change)
        if (debounceTimerRef.current) {
          clearTimeout(debounceTimerRef.current);
        }
        debounceTimerRef.current = setTimeout(() => {
          onChange?.(newValue);
        }, 100);
      },
      [step, clampValue, onChange],
    );

    const handleMouseUp = useCallback(() => {
      isDraggingRef.current = false;
      document.body.style.cursor = "default";

      // Clear debounce timer and immediately call onChange with final value
      if (debounceTimerRef.current) {
        clearTimeout(debounceTimerRef.current);
        debounceTimerRef.current = null;
      }

      // Immediately update store with final value when releasing
      if (draggingValue !== null) {
        onChange?.(draggingValue);
      }

      setDraggingValue(null);
    }, [draggingValue, onChange]);

    React.useEffect(() => {
      window.addEventListener("mousemove", handleMouseMove);
      window.addEventListener("mouseup", handleMouseUp);

      return () => {
        window.removeEventListener("mousemove", handleMouseMove);
        window.removeEventListener("mouseup", handleMouseUp);
        // Clean up debounce timer on unmount
        if (debounceTimerRef.current) {
          clearTimeout(debounceTimerRef.current);
        }
      };
    }, [handleMouseMove, handleMouseUp]);

    // Display the dragging value if dragging, otherwise the prop value
    const displayValue = draggingValue !== null ? draggingValue : value;

    return (
      <div
        ref={ref}
        className={cn(
          "inline-flex items-center justify-center cursor-ns-resize select-none",
          "font-mono text-sm font-medium",
          "px-3 py-1.5 rounded-lg",
          "bg-black/30 backdrop-blur-md",
          "border border-white/20",
          "shadow-lg",
          disabled && "opacity-50 cursor-not-allowed",
          className,
        )}
        onMouseDown={handleMouseDown}
      >
        {formatValue(displayValue)}
      </div>
    );
  },
);

NumberScrubber.displayName = "NumberScrubber";
