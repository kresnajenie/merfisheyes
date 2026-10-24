"use client";

import type { ImageChannel, ImageLayer } from "@/lib/spatialdata/ImageLayer";

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

const clamp = (v: number, lo: number, hi: number) =>
  Math.min(Math.max(v, lo), hi);

/**
 * A whole number typed in, applied on Enter or on leaving the box; Escape
 * puts back the value it had.
 */
function NumberBox({
  label,
  value,
  disabled,
  onCommit,
}: {
  label: string;
  value: number;
  disabled?: boolean;
  onCommit: (value: number) => void;
}) {
  const [draft, setDraft] = useState<string | null>(null);
  const commit = () => {
    const n = Number(draft);

    setDraft(null);
    if (draft !== null && draft.trim() !== "" && Number.isFinite(n)) {
      onCommit(Math.round(n));
    }
  };

  return (
    <input
      aria-label={label}
      className="w-14 shrink-0 rounded bg-default-100 px-1 py-0.5 text-right text-xs tabular-nums text-foreground disabled:opacity-50"
      disabled={disabled}
      inputMode="numeric"
      type="text"
      value={draft ?? String(value)}
      onBlur={commit}
      onChange={(e) => setDraft(e.target.value)}
      onKeyDown={(e) => {
        if (e.key === "Enter") e.currentTarget.blur();
        if (e.key === "Escape") setDraft(null);
      }}
    />
  );
}

/**
 * Image button + channel panel for the SpatialData viewer: show/hide the
 * image, and per channel its visibility, colour, display window and the
 * range its sliders cover (narrow it for finer control where the data is).
 * What is set here is kept in the link (img=, or rimg= in the right split
 * panel), so a shared link looks the same.
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

  // The sliders' span: the channel's own range, else the whole data type
  const span = (ch: ImageChannel): [number, number] =>
    ch.range ?? [0, layer.maxValue];

  // A bound dragged or typed. The bounds can't cross: the one being set
  // stops one short of the other. A typed value outside the range widens it.
  const setBound = (i: number, bound: "min" | "max", value: number) => {
    const ch = layer.channels[i];
    const v = clamp(value, 0, layer.maxValue);
    const [lo, hi] = span(ch);
    const patch: Partial<ImageChannel> =
      bound === "min"
        ? { min: Math.min(v, ch.max - 1) }
        : { max: Math.max(v, ch.min + 1) };

    if (ch.range && (v < lo || v > hi)) {
      patch.range = [Math.min(lo, v), Math.max(hi, v)];
    }
    layer.setChannel(i, patch);
    changed();
  };

  // One end of the sliders' range typed. The window is pulled inside it, so
  // the sliders show what is drawn; the whole data type is no range at all.
  const setRangeEnd = (i: number, end: 0 | 1, value: number) => {
    const ch = layer.channels[i];
    const v = clamp(value, 0, layer.maxValue);
    let [lo, hi] = span(ch);

    if (end === 0) lo = Math.min(v, hi - 1);
    else hi = Math.max(v, lo + 1);
    const min = clamp(ch.min, lo, hi - 1);
    const max = clamp(ch.max, min + 1, hi);
    const patch: Partial<ImageChannel> = {
      range: lo === 0 && hi === layer.maxValue ? undefined : [lo, hi],
    };

    if (min !== ch.min) patch.min = min;
    if (max !== ch.max) patch.max = max;
    layer.setChannel(i, patch);
    changed();
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
                </div>
                {/* Display window: values at or below the lower bound draw
                  black, at or above the upper one full colour. */}
                {(["min", "max"] as const).map((bound) => (
                  <div
                    key={bound}
                    className="flex items-center gap-2 text-xs text-default-500"
                  >
                    <span className="w-9 shrink-0">{bound}</span>
                    <input
                      aria-label={`${channel.label} ${bound === "min" ? "minimum" : "maximum"}`}
                      className="min-w-0 flex-1"
                      disabled={!channel.visible}
                      max={span(channel)[1]}
                      min={span(channel)[0]}
                      type="range"
                      value={channel[bound]}
                      onChange={(e) =>
                        setBound(i, bound, Number(e.target.value))
                      }
                    />
                    <NumberBox
                      disabled={!channel.visible}
                      label={`${channel.label} ${bound === "min" ? "minimum" : "maximum"} value`}
                      value={channel[bound]}
                      onCommit={(v) => setBound(i, bound, v)}
                    />
                  </div>
                ))}
                {/* What the sliders span: narrow it to where the data is */}
                <div className="flex items-center gap-2 text-xs text-default-500">
                  <span className="w-9 shrink-0">range</span>
                  <span className="flex-1" />
                  <NumberBox
                    disabled={!channel.visible}
                    label={`${channel.label} range start`}
                    value={span(channel)[0]}
                    onCommit={(v) => setRangeEnd(i, 0, v)}
                  />
                  <span>–</span>
                  <NumberBox
                    disabled={!channel.visible}
                    label={`${channel.label} range end`}
                    value={span(channel)[1]}
                    onCommit={(v) => setRangeEnd(i, 1, v)}
                  />
                </div>
              </div>
            ))}
          </FadeScroll>
        </div>
      )}
    </div>
  );
}
