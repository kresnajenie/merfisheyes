"use client";

import { useEffect, useState, useSyncExternalStore } from "react";

import { glassPanel } from "@/components/primitives";
import {
  describeLoad,
  loadSnapshot,
  subscribeLoad,
  type LoadSnapshot,
} from "@/lib/net/load-meter";

/** A burst shorter than this never shows: panning an image must not flicker. */
const SHOW_AFTER_MS = 400;

const serverSnapshot = (): LoadSnapshot => loadSnapshot();

/** The current download burst (lib/net/load-meter.ts), re-rendering as it moves. */
export function useLoadMeter(): LoadSnapshot {
  return useSyncExternalStore(subscribeLoad, loadSnapshot, serverSnapshot);
}

// Loading screens say the same thing in their own place; while one is up the
// floating bar stays away
let inlineCount = 0;
const inlineListeners = new Set<() => void>();
const subscribeInline = (l: () => void) => {
  inlineListeners.add(l);

  return () => inlineListeners.delete(l);
};
const hasInline = () => inlineCount > 0;

/**
 * The same megabytes as the bar, as a line of text for a loading screen
 * (which has a progress bar of its own). Empty until something is read.
 */
export function LoadMeterText({ className }: { className?: string }) {
  const load = useLoadMeter();

  useEffect(() => {
    inlineCount++;
    inlineListeners.forEach((l) => l());

    return () => {
      inlineCount--;
      inlineListeners.forEach((l) => l());
    };
  }, []);

  return (
    <p className={className} data-testid="load-meter-text">
      {load.started ? describeLoad(load).text : "\u00a0"}
    </p>
  );
}

/**
 * How much is being downloaded, at the bottom of the page: megabytes received
 * of megabytes announced, and files finished of files asked for. Mounted once
 * at the app root; it shows for any data read that lasts.
 */
export function LoadMeterBar() {
  const load = useLoadMeter();
  const [shown, setShown] = useState(false);
  const inline = useSyncExternalStore(subscribeInline, hasInline, () => false);

  useEffect(() => {
    if (!load.active) {
      setShown(false);

      return;
    }
    const wait = SHOW_AFTER_MS - (performance.now() - load.since);
    const timer = setTimeout(() => setShown(true), Math.max(0, wait));

    return () => clearTimeout(timer);
  }, [load.active, load.since]);

  if (!shown || !load.active || inline) return null;
  const { fraction, text } = describeLoad(load);

  return (
    <div
      data-ui-overlay
      aria-live="off"
      className={`pointer-events-none fixed bottom-4 left-1/2 z-[var(--z-chrome)] w-72 max-w-[calc(100vw-2rem)] -translate-x-1/2 px-4 py-2 ${glassPanel()}`}
      data-testid="load-meter"
      role="status"
    >
      <div
        aria-label="Download progress"
        aria-valuemax={100}
        aria-valuemin={0}
        aria-valuenow={Math.round(fraction * 100)}
        className="h-1 w-full overflow-hidden rounded-full bg-white/15"
        role="progressbar"
      >
        <div
          className="h-full rounded-full bg-white/80 transition-[width] duration-150"
          style={{ width: `${fraction * 100}%` }}
        />
      </div>
      <p className="mt-1 text-center text-[11px] tabular-nums text-white/80">
        {text}
      </p>
    </div>
  );
}
