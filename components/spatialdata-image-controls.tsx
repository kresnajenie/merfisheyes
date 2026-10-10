"use client";

import type { ImageLayer } from "@/lib/spatialdata/ImageLayer";

import { useEffect, useReducer, useState } from "react";

import { glassButton, glassPanel } from "@/components/primitives";

/**
 * Image button + channel panel for the SpatialData viewer: show/hide the
 * image, and per channel its visibility, colour and display window.
 */
export function SpatialDataImageControls({ layer }: { layer: ImageLayer }) {
  const [open, setOpen] = useState(false);
  const [, rerender] = useReducer((n: number) => n + 1, 0);

  // The layer sets each channel's window itself once it has seen the data
  useEffect(() => {
    layer.onChange = rerender;

    return () => {
      layer.onChange = undefined;
    };
  }, [layer]);

  return (
    <div className="absolute bottom-6 left-36 z-[var(--z-rail)] flex flex-col-reverse items-start gap-2">
      <button
        className={`${glassButton()} h-10 px-4 rounded-full text-sm`}
        type="button"
        onClick={() => setOpen((o) => !o)}
      >
        Image
      </button>
      {open && (
        <div className={`${glassPanel()} p-4 w-80 flex flex-col gap-3`}>
          <label className="flex items-center gap-2 text-sm">
            <input
              checked={layer.visible}
              type="checkbox"
              onChange={(e) => {
                layer.visible = e.target.checked;
                rerender();
              }}
            />
            Show image
          </label>
          {layer.channels.map((channel, i) => (
            <div key={i} className="flex flex-col gap-1">
              <div className="flex items-center gap-2 text-sm">
                <input
                  checked={channel.visible}
                  type="checkbox"
                  onChange={(e) => {
                    layer.setChannel(i, { visible: e.target.checked });
                    rerender();
                  }}
                />
                <input
                  aria-label={`${channel.label} colour`}
                  className="w-6 h-6 p-0 border-0 bg-transparent"
                  type="color"
                  value={channel.color}
                  onChange={(e) => {
                    layer.setChannel(i, { color: e.target.value });
                    rerender();
                  }}
                />
                <span className="truncate">{channel.label}</span>
                <span className="ml-auto text-xs text-default-500">
                  {channel.max}
                </span>
              </div>
              <input
                aria-label={`${channel.label} maximum`}
                disabled={!channel.visible}
                max={layer.maxValue}
                min={1}
                type="range"
                value={channel.max}
                onChange={(e) => {
                  layer.setChannel(i, { max: Number(e.target.value) });
                  rerender();
                }}
              />
            </div>
          ))}
        </div>
      )}
    </div>
  );
}
