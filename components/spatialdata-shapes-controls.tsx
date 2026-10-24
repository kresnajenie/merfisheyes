"use client";

import type { ShapesLayer, ShapesMode } from "@/lib/spatialdata/ShapesLayer";

import { Button } from "@heroui/button";
import { Tooltip } from "@heroui/react";
import { useEffect, useRef, useState } from "react";
import { createPortal } from "react-dom";

import { glassButton, glassPanel } from "@/components/primitives";
import { useClickOutside } from "@/lib/hooks/useClickOutside";
import { useLayerStatus } from "@/lib/scene/useLayerStatus";
import { SHAPES_SHOW_WIDTH } from "@/lib/spatialdata/ShapesLayer";

/**
 * Rail button + panel for cell polygons: show or hide them, outline or
 * filled, and their opacity. They are drawn once the view is narrow enough
 * to tell cells apart.
 */
export function SpatialDataShapesControls({ layer }: { layer: ShapesLayer }) {
  const [open, setOpen] = useState(false);
  const anchor = useRef<HTMLDivElement>(null);
  const panel = useRef<HTMLDivElement>(null);
  // Where the panel opens: beside its button, on the page. The button sits
  // in the rail, which scrolls and would clip a panel opened inside it.
  const [at, setAt] = useState<{ left: number; top: number } | null>(null);

  useLayerStatus(layer);
  // A click anywhere else closes it, as the other flyouts do
  useClickOutside([anchor, panel], () => setOpen(false), open);

  useEffect(() => {
    if (!open) return;
    const place = () => {
      const rect = anchor.current?.getBoundingClientRect();

      if (rect) setAt({ left: rect.right + 8, top: rect.top });
    };

    place();
    // The rail scrolled or the window changed: follow the button
    window.addEventListener("resize", place);
    window.addEventListener("scroll", place, true);

    return () => {
      window.removeEventListener("resize", place);
      window.removeEventListener("scroll", place, true);
    };
  }, [open]);

  const modes: { key: ShapesMode; label: string }[] = [
    { key: "outline", label: "Outline" },
    { key: "filled", label: "Filled" },
  ];
  const note =
    layer.shapesStatus === "loading"
      ? "Loading cell shapes…"
      : layer.shapesStatus === "error"
        ? "Couldn't load the cell shapes."
        : !layer.zoomedIn
          ? `Zoom in to under ${SHAPES_SHOW_WIDTH / 1000} mm across to see them.`
          : layer.visible
            ? "Drawn in each cell's colour."
            : "";

  return (
    <div ref={anchor} className="relative">
      <Tooltip content="Cell shapes" placement="right">
        <Button
          className={`w-14 h-14 min-w-0 rounded-full font-medium text-xs ${open ? "" : glassButton()}`}
          color={open ? "primary" : "default"}
          data-testid="shapes-button"
          variant={open ? "shadow" : "light"}
          onPress={() => setOpen((o) => !o)}
        >
          Shapes
        </Button>
      </Tooltip>
      {open &&
        at &&
        createPortal(
          <div
            ref={panel}
            data-ui-overlay
            className={`fixed z-[var(--z-panel)] w-72 p-4 flex flex-col gap-3 text-sm ${glassPanel()}`}
            data-testid="shapes-panel"
            style={{ left: at.left, top: at.top }}
          >
            <label className="flex items-center gap-2">
              <input
                checked={layer.visible}
                type="checkbox"
                onChange={(e) => {
                  layer.visible = e.target.checked;
                }}
              />
              Show cell shapes
            </label>
            <div className="flex items-center gap-2">
              <span className="text-default-500">Style</span>
              <div className="ml-auto flex overflow-hidden rounded-md border border-default-300">
                {modes.map((m) => (
                  <button
                    key={m.key}
                    className={`px-2.5 py-0.5 text-xs ${
                      layer.mode === m.key
                        ? "bg-default-200 text-foreground"
                        : "text-default-500 hover:bg-default-100"
                    }`}
                    disabled={!layer.visible}
                    type="button"
                    onClick={() => {
                      layer.mode = m.key;
                    }}
                  >
                    {m.label}
                  </button>
                ))}
              </div>
            </div>
            <label className="flex items-center gap-2 text-xs text-default-500">
              <span className="w-14">Opacity</span>
              <input
                aria-label="Shape opacity"
                className="flex-1"
                disabled={!layer.visible}
                max={100}
                min={5}
                type="range"
                value={Math.round(layer.opacity * 100)}
                onChange={(e) => {
                  layer.opacity = Number(e.target.value) / 100;
                }}
              />
              <span className="w-9 text-right tabular-nums">
                {Math.round(layer.opacity * 100)}%
              </span>
            </label>
            {note && <p className="text-xs text-default-500">{note}</p>}
          </div>,
          document.body,
        )}
    </div>
  );
}
