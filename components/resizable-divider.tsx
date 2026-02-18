"use client";

import { useCallback, useEffect, useRef, useState } from "react";
import { usePathname } from "next/navigation";
import { Tooltip } from "@heroui/tooltip";

import { useSplitScreenStore } from "@/lib/stores/splitScreenStore";

export function ResizableDivider() {
  const {
    dividerPosition,
    setDividerPosition,
    syncEnabled,
    toggleSync,
    rightPanelType,
  } = useSplitScreenStore();
  const pathname = usePathname();
  const [isDragging, setIsDragging] = useState(false);
  const containerRef = useRef<HTMLDivElement | null>(null);

  const handleMouseDown = useCallback((e: React.MouseEvent) => {
    e.preventDefault();
    setIsDragging(true);
  }, []);

  useEffect(() => {
    if (!isDragging) return;

    const handleMouseMove = (e: MouseEvent) => {
      const container = containerRef.current?.parentElement;

      if (!container) return;

      const rect = container.getBoundingClientRect();
      const percentage = ((e.clientX - rect.left) / rect.width) * 100;

      setDividerPosition(percentage);
    };

    const handleMouseUp = () => {
      setIsDragging(false);
    };

    document.addEventListener("mousemove", handleMouseMove);
    document.addEventListener("mouseup", handleMouseUp);

    return () => {
      document.removeEventListener("mousemove", handleMouseMove);
      document.removeEventListener("mouseup", handleMouseUp);
    };
  }, [isDragging, setDividerPosition]);

  const leftPanelType = pathname.startsWith("/sm-viewer") ? "sm" : "cell";
  const showSyncButton = rightPanelType !== null && leftPanelType === rightPanelType;

  return (
    <div
      ref={containerRef}
      className={`w-2 cursor-col-resize flex-shrink-0 relative group ${
        isDragging ? "bg-primary/40" : "bg-white/10 hover:bg-white/20"
      } transition-colors`}
      onMouseDown={handleMouseDown}
    >
      {/* Visual grip indicator */}
      <div className="absolute inset-y-0 left-1/2 -translate-x-1/2 w-0.5 bg-white/30 group-hover:bg-white/50 transition-colors" />

      {/* Sync toggle button */}
      {showSyncButton && (
        <Tooltip
          content={syncEnabled ? "Disable sync" : "Enable sync"}
          placement="right"
        >
          <button
            className={`absolute top-1/2 left-1/2 -translate-x-1/2 -translate-y-1/2 z-50 w-8 h-8 rounded-full flex items-center justify-center transition-all ${
              syncEnabled
                ? "bg-primary/80 shadow-[0_0_12px_rgba(0,112,243,0.5)] border border-primary/60"
                : "bg-black/60 border border-white/20 hover:bg-white/10"
            }`}
            onClick={(e) => {
              e.stopPropagation();
              toggleSync();
            }}
            onMouseDown={(e) => e.stopPropagation()}
          >
            <svg
              className={`w-4 h-4 ${syncEnabled ? "text-white" : "text-white/60"}`}
              fill="none"
              stroke="currentColor"
              strokeWidth={2}
              viewBox="0 0 24 24"
            >
              {/* Chain link icon */}
              <path
                d="M10 13a5 5 0 0 0 7.54.54l3-3a5 5 0 0 0-7.07-7.07l-1.72 1.71"
                strokeLinecap="round"
                strokeLinejoin="round"
              />
              <path
                d="M14 11a5 5 0 0 0-7.54-.54l-3 3a5 5 0 0 0 7.07 7.07l1.71-1.71"
                strokeLinecap="round"
                strokeLinejoin="round"
              />
            </svg>
          </button>
        </Tooltip>
      )}
    </div>
  );
}
