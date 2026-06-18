"use client";

import { useEffect, useRef, useState } from "react";
import * as THREE from "three";

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

  useEffect(() => {
    let running = true;

    const update = () => {
      if (!running) return;

      const camera = cameraRef.current;
      const renderer = rendererRef.current;

      if (camera && renderer) {
        const canvasHeight = renderer.domElement.clientHeight;
        const canvasWidth = renderer.domElement.clientWidth;
        // Distance from camera to look-at target (not origin — raw coords may be far from origin)
        const target = controlsRef?.current?.target;
        const cameraDistance = target
          ? camera.position.distanceTo(target)
          : camera.position.length();
        const fovRad = (camera.fov / 2) * (Math.PI / 180);

        // World units per pixel
        const worldPerPx = (2 * cameraDistance * Math.tan(fovRad)) / canvasHeight;

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

  if (!barWidth || !label) return null;

  return (
    <div className="absolute bottom-6 left-6 z-50 flex flex-col items-start gap-1">
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
