"use client";

import type { ImageChannel, ImageLayer } from "@/lib/spatialdata/ImageLayer";

import { useEffect, useReducer, useRef, useState } from "react";
import { toast } from "react-toastify";

import { FadeScroll } from "@/components/fade-scroll";
import { glassButton, glassPanel } from "@/components/primitives";
import { useClickOutside } from "@/lib/hooks/useClickOutside";
import { usePanelId, usePanelLayersApi } from "@/lib/hooks/usePanelStores";
import { useLayerStatus } from "@/lib/scene/useLayerStatus";
import {
  applyImageState,
  encodeImageChoice,
  encodeImageState,
  imageNameInState,
  overlayNameInState,
} from "@/lib/spatialdata/imageUrlState";
import {
  readUrlVizState,
  scheduleUrlUpdate,
} from "@/lib/utils/url-state-writer";

const clamp = (v: number, lo: number, hi: number) =>
  Math.min(Math.max(v, lo), hi);

/**
 * A number typed in, applied on Enter or on leaving the box; Escape puts back
 * the value it had. A whole number, unless the image's pixels are floats.
 */
function NumberBox({
  label,
  value,
  disabled,
  float = false,
  onCommit,
}: {
  label: string;
  value: number;
  disabled?: boolean;
  float?: boolean;
  onCommit: (value: number) => void;
}) {
  const [draft, setDraft] = useState<string | null>(null);
  const commit = () => {
    const n = Number(draft);

    setDraft(null);
    if (draft !== null && draft.trim() !== "" && Number.isFinite(n)) {
      onCommit(float ? n : Math.round(n));
    }
  };

  return (
    <input
      aria-label={label}
      className="w-14 shrink-0 rounded bg-default-100 px-1 py-0.5 text-right text-xs tabular-nums text-foreground disabled:opacity-50"
      disabled={disabled}
      inputMode={float ? "decimal" : "numeric"}
      type="text"
      value={draft ?? String(float ? Number(value.toPrecision(4)) : value)}
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
 * One image's own settings: greyscale for a colour photograph, and per
 * channel its visibility, colour, display window and sliders' range.
 */
function ImageSettings({
  layer,
  changed,
  tag = "",
}: {
  layer: ImageLayer;
  changed: () => void;
  /** Starts every control's name, to tell the image on top from the other. */
  tag?: string;
}) {
  // The sliders' span: the channel's own range, else the whole data type
  const span = (ch: ImageChannel): [number, number] =>
    ch.range ?? [0, layer.maxValue];
  // The least the two bounds stay apart, and the sliders' step: one for whole
  // numbers, a thousandth of the sliders' range for floats
  const gap = (ch: ImageChannel) =>
    layer.float ? (span(ch)[1] - span(ch)[0]) / 1000 : 1;

  // A bound dragged or typed. The bounds can't cross: the one being set
  // stops one short of the other. A typed value outside the range widens it.
  const setBound = (i: number, bound: "min" | "max", value: number) => {
    const ch = layer.channels[i];
    const v = clamp(value, 0, layer.maxValue);
    const [lo, hi] = span(ch);
    const patch: Partial<ImageChannel> =
      bound === "min"
        ? { min: Math.min(v, ch.max - gap(ch)) }
        : { max: Math.max(v, ch.min + gap(ch)) };

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

    const least = gap(ch);

    if (end === 0) lo = Math.min(v, hi - least);
    else hi = Math.max(v, lo + least);
    const min = clamp(ch.min, lo, hi - least);
    const max = clamp(ch.max, min + least, hi);
    const patch: Partial<ImageChannel> = {
      range: lo === 0 && hi === layer.maxValue ? undefined : [lo, hi],
    };

    if (min !== ch.min) patch.min = min;
    if (max !== ch.max) patch.max = max;
    layer.setChannel(i, patch);
    changed();
  };

  return (
    <>
      {layer.rgb && (
        <label className="flex items-center gap-2 text-sm">
          <input
            aria-label={`${tag}Greyscale`}
            checked={layer.grayscale}
            type="checkbox"
            onChange={(e) => {
              layer.grayscale = e.target.checked;
              changed();
            }}
          />
          Greyscale
        </label>
      )}
      {layer.channels.map((channel, i) => (
        <div key={i} className="flex flex-col gap-1">
          <div className="flex items-center gap-2 text-sm">
            <input
              aria-label={`${tag}Show ${channel.label}`}
              checked={channel.visible}
              type="checkbox"
              onChange={(e) => {
                layer.setChannel(i, { visible: e.target.checked });
                changed();
              }}
            />
            <input
              aria-label={`${tag}${channel.label} colour`}
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
                aria-label={`${tag}${channel.label} ${bound === "min" ? "minimum" : "maximum"}`}
                className="min-w-0 flex-1"
                disabled={!channel.visible}
                max={span(channel)[1]}
                min={span(channel)[0]}
                step={gap(channel)}
                type="range"
                value={channel[bound]}
                onChange={(e) => setBound(i, bound, Number(e.target.value))}
              />
              <NumberBox
                float={layer.float}
                disabled={!channel.visible}
                label={`${tag}${channel.label} ${bound === "min" ? "minimum" : "maximum"} value`}
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
              float={layer.float}
              disabled={!channel.visible}
              label={`${tag}${channel.label} range start`}
              value={span(channel)[0]}
              onCommit={(v) => setRangeEnd(i, 0, v)}
            />
            <span>–</span>
            <NumberBox
              float={layer.float}
              disabled={!channel.visible}
              label={`${tag}${channel.label} range end`}
              value={span(channel)[1]}
              onCommit={(v) => setRangeEnd(i, 1, v)}
            />
          </div>
        </div>
      ))}
    </>
  );
}

/**
 * Image button + channel panel for the SpatialData viewer: show/hide the
 * image, and per channel its visibility, colour, display window and the
 * range its sliders cover (narrow it for finer control where the data is).
 * A store with several images offers which one to draw, and another to
 * draw over it, seen through.
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
  const layersApi = usePanelLayersApi();

  // Another image of the same store, drawn in place of this one
  const switchTo = async (name: string) => {
    const sibling = layer.siblings.find((s) => s.name === name);

    if (!sibling || name === layer.name) return;
    try {
      layersApi.getState().replaceElement("image", await sibling.open());
    } catch (error) {
      toast.error(`Could not open ${name}: ${(error as Error).message}`);
    }
  };

  // Another image of the same store, drawn over this one
  const putOver = async (name: string | null) => {
    const sibling = layer.siblings.find((s) => s.name === name);

    if (!sibling || name === layer.name) return layer.setOverlay(null);
    try {
      layer.setOverlay(await sibling.open());
    } catch (error) {
      toast.error(`Could not open ${name}: ${(error as Error).message}`);
    }
  };

  // The layer pushes its changes (auto contrast once it has seen the data)
  useLayerStatus(layer);

  // Settings in the link are applied to a new layer; a window set there is
  // kept, any other the layer sets itself once it has seen the data.
  useEffect(() => {
    const fromUrl = readUrlVizState()[slot];
    const named = fromUrl && imageNameInState(fromUrl);

    if (!fromUrl) return rerender();
    // Written for another image of this store: open that one first
    if (named && named !== layer.name) return void switchTo(named);
    void (async () => {
      const over = overlayNameInState(fromUrl);

      if (over && layer.overlay?.name !== over) await putOver(over);
      applyImageState(layer, fromUrl);
      rerender();
    })();
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
            {layer.siblings.length > 1 && (
              <label className="flex items-center gap-2 text-sm">
                Image
                <select
                  aria-label="Image"
                  className="min-w-0 flex-1 rounded bg-default-100 px-1 py-0.5 text-xs text-foreground"
                  value={layer.name}
                  onChange={(e) => {
                    scheduleUrlUpdate(slot, encodeImageChoice(e.target.value));
                    void switchTo(e.target.value);
                  }}
                >
                  {layer.siblings.map((s) => (
                    <option key={s.name} value={s.name}>
                      {s.name}
                    </option>
                  ))}
                </select>
              </label>
            )}
            <ImageSettings changed={changed} layer={layer} />
            {layer.siblings.length > 1 && (
              <label className="flex items-center gap-2 text-sm">
                On top
                <select
                  aria-label="Image on top"
                  className="min-w-0 flex-1 rounded bg-default-100 px-1 py-0.5 text-xs text-foreground"
                  value={layer.overlay?.name ?? ""}
                  onChange={async (e) => {
                    await putOver(e.target.value || null);
                    changed();
                  }}
                >
                  <option value="">Nothing</option>
                  {layer.siblings
                    .filter((s) => s.name !== layer.name)
                    .map((s) => (
                      <option key={s.name} value={s.name}>
                        {s.name}
                      </option>
                    ))}
                </select>
              </label>
            )}
            {layer.overlay && (
              <>
                <label className="flex items-center gap-2 text-xs text-default-500">
                  <span className="w-12 shrink-0">opacity</span>
                  <input
                    aria-label="Opacity of the image on top"
                    className="min-w-0 flex-1"
                    max={1}
                    min={0}
                    step={0.01}
                    type="range"
                    value={layer.overlay.opacity}
                    onChange={(e) => {
                      layer.overlay!.opacity = Number(e.target.value);
                      changed();
                    }}
                  />
                </label>
                <ImageSettings
                  changed={changed}
                  layer={layer.overlay}
                  tag="On top: "
                />
              </>
            )}
          </FadeScroll>
        </div>
      )}
    </div>
  );
}
