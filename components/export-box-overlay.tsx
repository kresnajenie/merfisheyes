"use client";

import { useEffect, useRef, useState } from "react";
import type * as THREE from "three";

import { micronsPerPixel } from "@/lib/utils/world-units";

interface ExportBoxState {
  enabled: boolean;
  widthMm: number;
  heightMm: number;
  centerPx: { x: number; y: number } | null;
}

interface ExportBoxOverlayProps {
  cameraRef: React.RefObject<THREE.PerspectiveCamera | null>;
  rendererRef: React.RefObject<THREE.WebGLRenderer | null>;
  // Matches SpatialScaleBar's loose controlsRef typing (Orbit/Trackball etc.).
  // We only read `.target` off it.
  controlsRef: React.RefObject<any | null>;
  state: ExportBoxState;
  setCenterPx: (pos: { x: number; y: number } | null) => void;
  // When false (e.g. 3D view, or coords not in micron-space), nothing renders.
  canShow: boolean;
}

// Drives a draggable rectangle that represents an exact spatial area in mm,
// and exports the cropped canvas region as a PNG (download or clipboard).
export function ExportBoxOverlay({
  cameraRef,
  rendererRef,
  controlsRef,
  state,
  setCenterPx,
  canShow,
}: ExportBoxOverlayProps) {
  // Pixel size of the box (driven by mm × current pixels-per-micron).
  const [pxSize, setPxSize] = useState<{ w: number; h: number } | null>(null);
  // Pixel position of the box's top-left, relative to the renderer's container.
  // Falls back to viewport center when no centerPx is set.
  const [pxPos, setPxPos] = useState<{ left: number; top: number } | null>(
    null,
  );
  const [status, setStatus] = useState<string | null>(null);

  const rafRef = useRef<number>(0);
  const dragRef = useRef<{
    offsetX: number;
    offsetY: number;
  } | null>(null);
  const boxRef = useRef<HTMLDivElement | null>(null);

  // RAF loop: recompute pixel size + position from current camera + state.
  useEffect(() => {
    if (!canShow || !state.enabled) {
      setPxSize(null);
      setPxPos(null);
      return;
    }
    let running = true;
    const tick = () => {
      if (!running) return;
      const camera = cameraRef.current;
      const renderer = rendererRef.current;
      const el = boxRef.current;
      if (camera && renderer) {
        const canvas = renderer.domElement;
        const ch = canvas.clientHeight;
        const cw = canvas.clientWidth;
        const umPerPx = micronsPerPixel(camera, controlsRef.current, ch);
        if (umPerPx > 0) {
          const wPx = (state.widthMm * 1000) / umPerPx;
          const hPx = (state.heightMm * 1000) / umPerPx;
          setPxSize({ w: wPx, h: hPx });
          // centerPx is canvas-relative. Translate to offsetParent coords so
          // the rendered absolute-positioned div lands in the right place
          // regardless of where the canvas sits inside the page layout.
          const center = state.centerPx ?? { x: cw / 2, y: ch / 2 };
          const canvasRect = canvas.getBoundingClientRect();
          const parent = (el?.offsetParent ?? canvas.offsetParent) as
            | HTMLElement
            | null;
          const parentRect = parent
            ? parent.getBoundingClientRect()
            : { left: 0, top: 0 };
          const offsetX = canvasRect.left - parentRect.left;
          const offsetY = canvasRect.top - parentRect.top;
          setPxPos({
            left: offsetX + center.x - wPx / 2,
            top: offsetY + center.y - hPx / 2,
          });
        }
      }
      rafRef.current = requestAnimationFrame(tick);
    };
    rafRef.current = requestAnimationFrame(tick);
    return () => {
      running = false;
      cancelAnimationFrame(rafRef.current);
    };
  }, [
    canShow,
    state.enabled,
    state.widthMm,
    state.heightMm,
    state.centerPx,
    cameraRef,
    rendererRef,
    controlsRef,
  ]);

  // Drag the box around. centerPx is canvas-relative so the box keeps its
  // intended on-screen position regardless of where the canvas sits in the
  // page layout.
  useEffect(() => {
    const onMove = (e: MouseEvent) => {
      const d = dragRef.current;
      const renderer = rendererRef.current;
      if (!d || !renderer || !pxSize) return;
      const canvasRect = renderer.domElement.getBoundingClientRect();
      // Mouse position translated into canvas-relative coords, anchored at
      // the drag-grab point. Then add half the box dimensions to get the
      // center.
      const canvasX = e.clientX - canvasRect.left - d.offsetX;
      const canvasY = e.clientY - canvasRect.top - d.offsetY;
      setCenterPx({
        x: canvasX + pxSize.w / 2,
        y: canvasY + pxSize.h / 2,
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
  }, [pxSize, setCenterPx, rendererRef]);

  const onBoxMouseDown = (e: React.MouseEvent) => {
    const el = boxRef.current;
    if (!el) return;
    const rect = el.getBoundingClientRect();
    dragRef.current = {
      offsetX: e.clientX - rect.left,
      offsetY: e.clientY - rect.top,
    };
    document.body.style.userSelect = "none";
    e.preventDefault();
  };

  const exportCroppedBlob = async (): Promise<Blob | null> => {
    const renderer = rendererRef.current;
    const el = boxRef.current;
    if (!renderer || !pxSize || !el) return null;
    const canvas = renderer.domElement;
    // Compute the box's bounding rect relative to the canvas (CSS px).
    const canvasRect = canvas.getBoundingClientRect();
    const boxRect = el.getBoundingClientRect();
    const cssLeft = boxRect.left - canvasRect.left;
    const cssTop = boxRect.top - canvasRect.top;
    // Renderer is constructed with preserveDrawingBuffer: true, so we can
    // read the GL buffer directly via toBlob.
    const fullBlob: Blob | null = await new Promise((resolve) =>
      canvas.toBlob((b) => resolve(b), "image/png"),
    );
    if (!fullBlob) return null;
    const bitmap = await createImageBitmap(fullBlob);
    // Canvas is rendered at devicePixelRatio × CSS pixels; clamp source rect
    // to the bitmap so we never read outside it.
    const dpr = bitmap.width / canvas.clientWidth;
    const sx = Math.max(0, Math.round(cssLeft * dpr));
    const sy = Math.max(0, Math.round(cssTop * dpr));
    const sw = Math.min(
      bitmap.width - sx,
      Math.max(1, Math.round(pxSize.w * dpr)),
    );
    const sh = Math.min(
      bitmap.height - sy,
      Math.max(1, Math.round(pxSize.h * dpr)),
    );
    const out = document.createElement("canvas");
    out.width = sw;
    out.height = sh;
    const ctx = out.getContext("2d");
    if (!ctx) return null;
    ctx.drawImage(bitmap, sx, sy, sw, sh, 0, 0, sw, sh);

    // Overlay the spatial scale bar onto the cropped image. We read the live
    // bar's DOM rect relative to the renderer's container; if the bar is
    // outside the box crop, draw it as a small overlay in the bottom-left
    // corner of the export so the figure always carries a scale reference.
    drawScaleBarOnCrop(ctx, out.width, out.height, sw, dpr, {
      umPerPx: micronsPerPixel(
        cameraRef.current!,
        controlsRef.current,
        canvas.clientHeight,
      ),
    });

    return await new Promise<Blob | null>((resolve) =>
      out.toBlob((b) => resolve(b), "image/png"),
    );
  };

  const flashStatus = (msg: string) => {
    setStatus(msg);
    setTimeout(() => setStatus(null), 1800);
  };

  const onSavePng = async () => {
    try {
      const blob = await exportCroppedBlob();
      if (!blob) return;
      const url = URL.createObjectURL(blob);
      const a = document.createElement("a");
      const stamp = new Date()
        .toISOString()
        .slice(0, 19)
        .replace(/[:T]/g, "-");
      a.href = url;
      a.download = `box-${state.widthMm}x${state.heightMm}mm-${stamp}.png`;
      document.body.appendChild(a);
      a.click();
      a.remove();
      URL.revokeObjectURL(url);
      flashStatus("Saved");
    } catch (e) {
      // eslint-disable-next-line no-console
      console.error("Export PNG failed", e);
      flashStatus("Failed");
    }
  };

  const onCopy = async () => {
    try {
      const blob = await exportCroppedBlob();
      if (!blob) return;
      if (!navigator.clipboard || !window.ClipboardItem) {
        flashStatus("Clipboard unsupported");
        return;
      }
      await navigator.clipboard.write([
        new ClipboardItem({ "image/png": blob }),
      ]);
      flashStatus("Copied");
    } catch (e) {
      // eslint-disable-next-line no-console
      console.error("Copy to clipboard failed", e);
      flashStatus("Failed");
    }
  };

  if (!canShow || !state.enabled || !pxSize || !pxPos) return null;

  return (
    <div
      ref={boxRef}
      className="absolute z-40 border-2 border-yellow-300 cursor-move"
      style={{
        left: `${pxPos.left}px`,
        top: `${pxPos.top}px`,
        width: `${pxSize.w}px`,
        height: `${pxSize.h}px`,
        boxShadow: "0 0 0 1px rgba(0,0,0,0.4)",
        // No background so the live scene shows through.
      }}
      onMouseDown={onBoxMouseDown}
    >
      {/* Floating toolbar — buttons stop propagation so they don't drag. */}
      <div
        className="absolute -top-9 right-0 flex items-center gap-1 px-1.5 py-1 rounded-md bg-black/70 backdrop-blur-sm border border-white/15"
        onMouseDown={(e) => e.stopPropagation()}
      >
        <button
          aria-label="Save PNG"
          className="text-[11px] text-white px-2 py-0.5 rounded hover:bg-white/15"
          title="Save PNG"
          type="button"
          onClick={onSavePng}
        >
          Save
        </button>
        <button
          aria-label="Copy to clipboard"
          className="text-[11px] text-white px-2 py-0.5 rounded hover:bg-white/15"
          title="Copy to clipboard"
          type="button"
          onClick={onCopy}
        >
          Copy
        </button>
        {status && (
          <span className="text-[10px] text-white/80 pl-1">{status}</span>
        )}
      </div>
      {/* Dimensions label */}
      <div
        className="absolute -bottom-5 left-0 text-[10px] text-yellow-200 bg-black/60 px-1.5 py-0.5 rounded pointer-events-none"
        style={{ whiteSpace: "nowrap" }}
      >
        {state.widthMm} × {state.heightMm} mm
      </div>
    </div>
  );
}

// Draw a small scale bar on the export crop. Mirrors spatial-scale-bar.tsx's
// "nice number" logic so the bar reads the same units as the live overlay,
// but renders independently inside the cropped image (so figures always
// carry a scale reference regardless of where the user dragged the box).
const NICE_STEPS = [
  0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000,
  20000, 50000,
];

function niceNumber(value: number): number {
  let result = NICE_STEPS[0];
  for (const step of NICE_STEPS) {
    if (step <= value) result = step;
    else break;
  }
  return result;
}

function formatLabel(microns: number): string {
  if (microns >= 1000) {
    const mm = microns / 1000;
    return mm % 1 === 0 ? `${mm} mm` : `${mm.toFixed(1)} mm`;
  }
  if (microns < 1) return `${microns * 1000} nm`;
  return `${microns} μm`;
}

function drawScaleBarOnCrop(
  ctx: CanvasRenderingContext2D,
  outW: number,
  outH: number,
  _sw: number,
  dpr: number,
  args: { umPerPx: number },
) {
  if (!args.umPerPx) return;
  // Target bar ~20% of crop width, in CSS pixels.
  const cssW = outW / dpr;
  const targetPx = Math.min(120, cssW * 0.2);
  const targetMicrons = args.umPerPx * targetPx;
  const niceMicrons = niceNumber(targetMicrons);
  const cssBarW = niceMicrons / args.umPerPx;
  const barW = cssBarW * dpr;
  const pad = 12 * dpr;
  const barH = 2 * dpr;
  const x = pad;
  const y = outH - pad;
  // Shadow then bar
  ctx.fillStyle = "rgba(0,0,0,0.55)";
  ctx.fillRect(x - 1, y - 1, barW + 2, barH + 2);
  ctx.fillStyle = "white";
  ctx.fillRect(x, y, barW, barH);
  // Label
  const fontPx = 11 * dpr;
  ctx.font = `${fontPx}px system-ui, -apple-system, sans-serif`;
  ctx.textBaseline = "bottom";
  const label = formatLabel(niceMicrons);
  const textY = y - 3 * dpr;
  ctx.fillStyle = "rgba(0,0,0,0.7)";
  ctx.fillText(label, x + 1, textY + 1);
  ctx.fillStyle = "white";
  ctx.fillText(label, x, textY);
}
