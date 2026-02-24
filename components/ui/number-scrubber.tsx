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
    const didDragRef = useRef(false);
    const startXRef = useRef(0);
    const startValueRef = useRef(0);

    // Track the live dragging value for smooth display updates
    const [draggingValue, setDraggingValue] = useState<number | null>(null);

    // Inline editing state
    const [isEditing, setIsEditing] = useState(false);
    const [editText, setEditText] = useState("");
    const inputRef = useRef<HTMLInputElement>(null);

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

    const commitEdit = useCallback(() => {
      const parsed = parseFloat(editText);

      if (Number.isFinite(parsed)) {
        onChange?.(clampValue(parsed));
      }
      setIsEditing(false);
    }, [editText, clampValue, onChange]);

    const startEditing = useCallback(() => {
      if (disabled) return;
      setEditText(formatValue(value));
      setIsEditing(true);
      // Focus the input after React renders it
      requestAnimationFrame(() => inputRef.current?.select());
    }, [disabled, value, formatValue]);

    const handleMouseDown = useCallback(
      (e: React.MouseEvent) => {
        if (disabled || isEditing) return;
        e.preventDefault();
        isDraggingRef.current = true;
        didDragRef.current = false;
        startXRef.current = e.clientY;
        startValueRef.current = value;
        setDraggingValue(value);
        document.body.style.cursor = "ns-resize";
      },
      [disabled, isEditing, value],
    );

    const handleMouseMove = useCallback(
      (e: MouseEvent) => {
        if (!isDraggingRef.current) return;
        const deltaY = e.clientY - startXRef.current;

        // Only count as a drag if moved more than 2px
        if (Math.abs(deltaY) > 2) {
          didDragRef.current = true;
        }

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
      const wasDragging = isDraggingRef.current;
      const didDrag = didDragRef.current;

      isDraggingRef.current = false;
      didDragRef.current = false;
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

      // If mouse went down and up without dragging, enter edit mode
      if (wasDragging && !didDrag) {
        startEditing();
      }
    }, [draggingValue, onChange, startEditing]);

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

    if (isEditing) {
      return (
        <input
          ref={inputRef}
          className={cn(
            "font-mono text-sm font-medium text-center",
            "px-1 py-1 rounded-lg w-full",
            "bg-black/50 backdrop-blur-md",
            "border border-blue-400/60",
            "shadow-lg outline-none",
            "text-white",
            className,
          )}
          type="text"
          inputMode="decimal"
          value={editText}
          onChange={(e) => setEditText(e.target.value)}
          onKeyDown={(e) => {
            if (e.key === "Enter") commitEdit();
            if (e.key === "Escape") setIsEditing(false);
          }}
          onBlur={commitEdit}
        />
      );
    }

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
