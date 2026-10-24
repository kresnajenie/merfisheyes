"use client";

import type { RefObject } from "react";

import { useMemo, useRef, useState } from "react";
import { toast } from "react-toastify";
import qrcode from "qrcode-generator";

import { glassPanel } from "@/components/primitives";
import { useClickOutside } from "@/lib/hooks/useClickOutside";
import { flushUrlState } from "@/lib/utils/url-state-writer";

/**
 * The Share button's state: `toggle` opens the panel on the link to this view
 * as it is now (pending URL writes flushed first) and copies it, or closes it.
 */
export function useShareLink() {
  const [url, setUrl] = useState<string | null>(null);

  const toggle = async () => {
    if (url) {
      setUrl(null);

      return;
    }
    flushUrlState();
    const href = window.location.href;

    setUrl(href);
    try {
      await navigator.clipboard.writeText(href);
    } catch {
      toast.error("Couldn't copy link — copy it from the share panel");
    }
  };

  return { url, toggle, close: () => setUrl(null) };
}

/** An SVG QR code of `text`, black on white with a quiet zone, as markup. */
export function qrSvg(text: string): string {
  // Type 0: the smallest version that holds it; M: 15% may be damaged
  const qr = qrcode(0, "M");

  qr.addData(text);
  qr.make();

  return qr.createSvgTag({ cellSize: 4, margin: 4, scalable: true });
}

/**
 * The share panel under the viewer's top-right buttons: the link to this
 * view (already copied when it opens) and its QR code, to open the same view
 * on a phone or from a slide. Closes on a click outside, like the camera
 * panel. `url` is read when it opens, so it is the view as it is now.
 */
export function SharePanel({
  url,
  onClose,
  controlsRef,
}: {
  url: string;
  onClose: () => void;
  /** The button row it opens from: a click there is not "outside". */
  controlsRef?: RefObject<HTMLElement | null>;
}) {
  const panel = useRef<HTMLDivElement>(null);
  const [copied, setCopied] = useState(true);
  const svg = useMemo(() => qrSvg(url), [url]);

  useClickOutside(controlsRef ? [panel, controlsRef] : [panel], onClose);

  const copy = async () => {
    try {
      await navigator.clipboard.writeText(url);
      setCopied(true);
    } catch {
      toast.error("Couldn't copy link");
    }
  };

  return (
    <div
      ref={panel}
      className={`absolute top-14 right-0 z-[var(--z-panel)] max-h-[calc(100dvh-6rem)] w-[260px] overflow-y-auto ${glassPanel()}`}
      data-testid="share-panel"
    >
      <div className="flex flex-col gap-3 p-4">
        <div className="flex items-center justify-between">
          <span className="text-sm font-medium">Share this view</span>
          <button
            aria-label="Close share panel"
            className="text-default-400 hover:text-foreground"
            type="button"
            onClick={onClose}
          >
            ✕
          </button>
        </div>
        {/* White behind it whatever the theme: a QR code needs the contrast */}
        <div
          aria-label="QR code of the link to this view"
          className="rounded-lg bg-white p-2 [&>svg]:h-auto [&>svg]:w-full"
          data-testid="share-qr"
          role="img"
          // The SVG is generated here from the URL, never from user markup
          // eslint-disable-next-line react/no-danger
          dangerouslySetInnerHTML={{ __html: svg }}
        />
        <div className="flex items-center gap-2">
          <input
            readOnly
            aria-label="Link to this view"
            className="min-w-0 flex-1 truncate rounded-md bg-default-100 px-2 py-1 text-xs text-default-600"
            value={url}
            onFocus={(e) => e.currentTarget.select()}
          />
          <button
            className="shrink-0 rounded-md bg-primary px-2 py-1 text-xs font-medium text-primary-foreground"
            type="button"
            onClick={copy}
          >
            {copied ? "Copied" : "Copy"}
          </button>
        </div>
        <p className="text-[11px] text-default-400">
          Scan to open this view on another device.
        </p>
      </div>
    </div>
  );
}
