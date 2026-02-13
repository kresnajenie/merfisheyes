"use client";

import { useCallback, useEffect, useRef, useState } from "react";

import { useSplitScreenStore } from "@/lib/stores/splitScreenStore";

export function ResizableDivider() {
  const { dividerPosition, setDividerPosition } = useSplitScreenStore();
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
    </div>
  );
}
