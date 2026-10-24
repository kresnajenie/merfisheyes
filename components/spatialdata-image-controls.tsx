"use client";

import type { ImageLayer } from "@/lib/spatialdata/ImageLayer";

import { useEffect, useReducer, useRef, useState } from "react";

import { FadeScroll } from "@/components/fade-scroll";
import { glassButton, glassPanel } from "@/components/primitives";
import { useClickOutside } from "@/lib/hooks/useClickOutside";
import { usePanelId } from "@/lib/hooks/usePanelStores";
import { useLayerStatus } from "@/lib/scene/useLayerStatus";
import {
  applyImageState,
  encodeImageState,
} from "@/lib/spatialdata/imageUrlState";
import {
  readUrlVizState,
  scheduleUrlUpdate,
} from "@/lib/utils/url-state-writer";

/**
 * Image button + channel panel for the SpatialData viewer: show/hide the
 * image, and per channel its visibility, colour and display window. What is
 * set here is kept in the link (img=, or rimg= in the right split panel), so
 * a shared link looks the same.
 */
export function SpatialDataImageControls({ layer }: { layer: ImageLayer }) {
  const [open, setOpen] = useState(false);
  const [, rerender] = useReducer((n: number) => n + 1, 0);
  const root = useRef<HTMLDivElement>(null);

  // A click anywhere else closes the panel, as the other flyouts do
  useClickOutside([root], () => setOpen(false), open);

  const slot = usePanelId() === "right" ? "rightImage" : "image";

  // The layer pushes its changes (auto contrast once it has seen the data)
  useLayerStatus(layer);

  // Settings in the link are applied to a new layer; a window set there is
  // kept, any other the layer sets itself once it has seen the data.
  useEffect(() => {
    const fromUrl = readUrlVizState()[slot];

    if (fromUrl) applyImageState(layer, fromUrl);
    rerender();
  }, [layer, slot]);

  // The right panel's slot goes away with the panel.
  useEffect(() => {
    if (slot !== "rightImage") return;

    return () => scheduleUrlUpdate(slot, null);
  }, [slot]);

  // Only what someone set goes in the link, so an untouched image leaves it
  // as it was.
  const changed = () => {
    scheduleUrlUpdate(slot, encodeImageState(layer));
    rerender();
  };

  return (
    // Sits between the brand pill and the control rail's first button
    <div
      ref={root}
      className="absolute top-[4.25rem] left-4 z-[var(--z-panel)] flex flex-row items-start gap-2"
    >
      <button
        className={`${glassButton()} h-10 w-14 rounded-full text-xs`}
        type="button"
        onClick={() => setOpen((o) => !o)}
      >
        Image
      </button>
      {open && (
        // Many channels (a CODEX image has dozens) scroll inside the panel
        <div
          className={`${glassPanel()} p-4 w-80 max-h-[calc(100dvh-6rem)] flex flex-col`}
          data-testid="image-panel"
        >
          <FadeScroll className="flex flex-col gap-3 min-h-0">
            <label className="flex items-center gap-2 text-sm">
              <input
                checked={layer.visible}
                type="checkbox"
                onChange={(e) => {
                  layer.visible = e.target.checked;
                  changed();
                }}
              />
              Show image
            </label>
            {layer.channels.map((channel, i) => (
              <div key={i} className="flex flex-col gap-1">
                <div className="flex items-center gap-2 text-sm">
                  <input
                    aria-label={`Show ${channel.label}`}
                    checked={channel.visible}
                    type="checkbox"
                    onChange={(e) => {
                      layer.setChannel(i, { visible: e.target.checked });
                      changed();
                    }}
                  />
                  <input
                    aria-label={`${channel.label} colour`}
                    className="w-6 h-6 p-0 border-0 bg-transparent"
                    type="color"
                    value={channel.color}
                    onChange={(e) => {
                      layer.setChannel(i, { color: e.target.value });
                      changed();
                    }}
                  />
                  <span className="truncate">{channel.label}</span>
                  <span className="ml-auto shrink-0 whitespace-nowrap text-xs text-default-500">
                    {channel.min} – {channel.max}
                  </span>
                </div>
                {/* Display window: values at or below the lower bound draw
                  black, at or above the upper one full colour. */}
                {(["min", "max"] as const).map((bound) => (
                  <label
                    key={bound}
                    className="flex items-center gap-2 text-xs text-default-500"
                  >
                    <span className="w-7">{bound}</span>
                    <input
                      aria-label={`${channel.label} ${bound === "min" ? "minimum" : "maximum"}`}
                      className="flex-1"
                      disabled={!channel.visible}
                      max={layer.maxValue}
                      min={0}
                      type="range"
                      value={channel[bound]}
                      onChange={(e) => {
                        const value = Number(e.target.value);

                        // The bounds can't cross: the one being dragged stops
                        // one short of the other.
                        layer.setChannel(
                          i,
                          bound === "min"
                            ? { min: Math.min(value, channel.max - 1) }
                            : { max: Math.max(value, channel.min + 1) },
                        );
                        changed();
                      }}
                    />
                  </label>
                ))}
              </div>
            ))}
          </FadeScroll>
        </div>
      )}
    </div>
  );
}
