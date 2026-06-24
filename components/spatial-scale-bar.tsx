"use client";

import { useEffect, useRef, useState } from "react";
import * as THREE from "three";

import { micronsPerPixel } from "@/lib/utils/world-units";

interface SpatialScaleBarProps {
  cameraRef: React.RefObject<THREE.PerspectiveCamera | null>;
  rendererRef: React.RefObject<THREE.WebGLRenderer | null>;
  controlsRef: React.RefObject<any | null>;
}

/** Pick the largest "nice" number that fits within the target value */
function niceNumber(value: number): number {
  const steps = [0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000];

  // Find the largest step that's <= value
  let result = steps[0];

  for (const step of steps) {
    if (step <= value) {
      result = step;
    } else {
      break;
    }
  }

  return result;
}

function formatLabel(microns: number): string {
  if (microns >= 1000) {
    const mm = microns / 1000;

    return mm % 1 === 0 ? `${mm} mm` : `${mm.toFixed(1)} mm`;
  }

  if (microns < 1) {
    return `${microns * 1000} nm`;
  }

  return `${microns} μm`;
}

export function SpatialScaleBar({ cameraRef, rendererRef, controlsRef }: SpatialScaleBarProps) {
  const [barWidth, setBarWidth] = useState(0);
  const [label, setLabel] = useState("");
  const rafRef = useRef<number>(0);

  // Custom drag-position (top/left in px). null = default bottom-left placement.
  const [pos, setPos] = useState<{ x: number; y: number } | null>(null);
  const dragRef = useRef<{ offsetX: number; offsetY: number } | null>(null);
  const elRef = useRef<HTMLDivElement | null>(null);

  useEffect(() => {
    let running = true;

    const update = () => {
      if (!running) return;

      const camera = cameraRef.current;
      const renderer = rendererRef.current;

      if (camera && renderer) {
        const canvasHeight = renderer.domElement.clientHeight;
        const canvasWidth = renderer.domElement.clientWidth;
        const worldPerPx = micronsPerPixel(
          camera,
          controlsRef?.current,
          canvasHeight,
        );

        // Target bar width ~120px, find nice micron value
        const targetPx = 120;
        const targetMicrons = worldPerPx * targetPx;
        const niceMicrons = niceNumber(targetMicrons);

        // Actual pixel width for the nice number
        const actualPx = niceMicrons / worldPerPx;

        // Clamp bar to reasonable screen size
        const clampedPx = Math.max(40, Math.min(actualPx, canvasWidth * 0.3));

        setBarWidth(clampedPx);
        setLabel(formatLabel(niceMicrons));
      }

      rafRef.current = requestAnimationFrame(update);
    };

    rafRef.current = requestAnimationFrame(update);

    return () => {
      running = false;
      cancelAnimationFrame(rafRef.current);
    };
  }, [cameraRef, rendererRef]);

  // Drag handlers — positions are stored relative to the offsetParent so the
  // bar stays inside its panel container (matters in split-screen mode where
  // each panel renders its own scale bar).
  useEffect(() => {
    const onMove = (e: MouseEvent) => {
      const d = dragRef.current;
      const el = elRef.current;
      if (!d || !el) return;
      const parent = el.offsetParent as HTMLElement | null;
      const parentRect = parent
        ? parent.getBoundingClientRect()
        : { left: 0, top: 0 };
      setPos({
        x: e.clientX - parentRect.left - d.offsetX,
        y: e.clientY - parentRect.top - d.offsetY,
      });
    };
    const onUp = () => {
      dragRef.current = null;
      document.body.style.userSelect = "";
    };
    window.addEventListener("mousemove", onMove);
    window.addEventListener("mouseup", onUp);
    return () => {
      window.removeEventListener("mousemove", onMove);
      window.removeEventListener("mouseup", onUp);
    };
  }, []);

  const onMouseDown = (e: React.MouseEvent) => {
    const el = elRef.current;
    if (!el) return;
    const rect = el.getBoundingClientRect();
    const parent = el.offsetParent as HTMLElement | null;
    const parentRect = parent
      ? parent.getBoundingClientRect()
      : { left: 0, top: 0 };
    dragRef.current = {
      offsetX: e.clientX - rect.left,
      offsetY: e.clientY - rect.top,
    };
    // Snap current position into px so subsequent moves work from a known origin.
    setPos({
      x: rect.left - parentRect.left,
      y: rect.top - parentRect.top,
    });
    document.body.style.userSelect = "none";
    e.preventDefault();
  };

  if (!barWidth || !label) return null;

  const positionStyle: React.CSSProperties = pos
    ? { top: `${pos.y}px`, left: `${pos.x}px` }
    : { bottom: "1.5rem", left: "1.5rem" };

  return (
    <div
      ref={elRef}
      className="absolute z-50 flex flex-col items-start gap-1 cursor-move select-none"
      style={positionStyle}
      title="Drag to reposition"
      onMouseDown={onMouseDown}
    >
      <div
        className="h-[2px] bg-white"
        style={{ width: `${barWidth}px` }}
      />
      <span className="text-white text-xs font-medium drop-shadow-[0_1px_2px_rgba(0,0,0,0.8)]">
        {label}
      </span>
    </div>
  );
}
