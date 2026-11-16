import type { ImageLayer } from "./ImageLayer";

import {
  fromBase64Url,
  toBase64Url,
} from "@/lib/utils/url-visualization-state";

/**
 * An image layer's display settings as they go in a link (`img=`): whether
 * the image is shown, and per channel whether it is on, its colour, its
 * display window and its sliders' range. Enough for a shared link to look
 * the same.
 */
interface ImageUrlState {
  /** Only present, as false, when the whole image is hidden. */
  s?: false;
  /** Which image of the store, where it has more than one to choose from. */
  n?: string;
  /** Only present, as 1, when a colour photograph is drawn in grey. */
  g?: 1;
  /**
   * Per channel: on, colour, and its window (lower, upper bound) once it has
   * one, then its sliders' range when narrowed. A channel never looked at
   * has no window yet, and gets it from the data.
   */
  c: (
    | [0 | 1, string]
    | [0 | 1, string, number, number]
    | [0 | 1, string, number, number, number, number]
  )[];
}

const HEX_RE = /^#[0-9a-fA-F]{6}$/;

export function encodeImageState(layer: ImageLayer): string {
  const state: ImageUrlState = {
    c: layer.channels.map((ch, i) =>
      ch.range
        ? [ch.visible ? 1 : 0, ch.color, ch.min, ch.max, ...ch.range]
        : layer.hasWindow(i)
          ? [ch.visible ? 1 : 0, ch.color, ch.min, ch.max]
          : [ch.visible ? 1 : 0, ch.color],
    ),
  };

  if (!layer.visible) state.s = false;
  if (layer.siblings.length > 1) state.n = layer.name;
  if (layer.grayscale) state.g = 1;

  return toBase64Url(state);
}

/** Settings that only name an image: what a link holds right after a switch. */
export function encodeImageChoice(name: string): string {
  return toBase64Url({ n: name, c: [] } satisfies ImageUrlState);
}

/** The image the settings in a link were written for, when they name one. */
export function imageNameInState(encoded: string): string | undefined {
  const name = fromBase64Url<ImageUrlState>(encoded)?.n;

  return typeof name === "string" ? name : undefined;
}

/**
 * Apply settings from a link. Ones written for an image with a different
 * number of channels, or that do not make sense for this one, are ignored.
 */
export function applyImageState(layer: ImageLayer, encoded: string): boolean {
  const state = fromBase64Url<ImageUrlState>(encoded);

  if (
    !Array.isArray(state?.c) ||
    state.c.length !== layer.channels.length ||
    !state.c.every(
      (ch) =>
        Array.isArray(ch) &&
        HEX_RE.test(ch[1]) &&
        (ch.length === 2 ||
          (Number.isFinite(ch[2]) &&
            Number.isFinite(ch[3]) &&
            ch[2] >= 0 &&
            ch[2] < ch[3] &&
            ch[3] <= layer.maxValue &&
            (ch.length === 4 ||
              (ch.length === 6 &&
                Number.isFinite(ch[4]) &&
                Number.isFinite(ch[5]) &&
                ch[4] >= 0 &&
                ch[4] <= ch[2] &&
                ch[3] <= ch[5] &&
                ch[5] <= layer.maxValue)))),
    )
  ) {
    return false;
  }
  layer.visible = state.s !== false;
  layer.grayscale = state.g === 1;
  state.c.forEach(([on, color, min, max, lo, hi], i) =>
    layer.setChannel(i, {
      visible: !!on,
      color,
      ...(min !== undefined && { min, max }),
      ...(lo !== undefined && { range: [lo, hi!] as [number, number] }),
    }),
  );

  return true;
}
