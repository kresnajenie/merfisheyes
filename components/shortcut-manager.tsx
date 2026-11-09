"use client";

import type { ImageLayer } from "@/lib/spatialdata/ImageLayer";
import type { ShapesLayer } from "@/lib/spatialdata/ShapesLayer";

import { useCallback, useEffect, useState } from "react";
import { createPortal } from "react-dom";

import { glassPanel } from "@/components/primitives";
import {
  usePanelLayers,
  usePanelSingleMoleculeVisualizationStore,
  usePanelVisualizationStore,
} from "@/lib/hooks/usePanelStores";

/**
 * Keyboard shortcuts for what the scene draws, and the list of them (`?`).
 *
 * Mounted once at the app root, so it reads the global stores — the main
 * viewer, which is the left panel in a split. The split's right panel has
 * its own stores (PanelProvider) and keeps its own buttons; keys do not
 * reach it.
 *
 * `H` (hide every overlay) lives in HideUiManager, which owns that state.
 * It is listed here so one place answers "what can I press".
 */

/** A key, what it toggles, and whether this dataset has that to toggle. */
interface Shortcut {
  key: string;
  label: string;
}

const SHORTCUTS: Shortcut[] = [
  { key: "C", label: "Cells" },
  { key: "M", label: "Molecules" },
  { key: "I", label: "Image" },
  { key: "P", label: "Cell polygons" },
  { key: "H", label: "Hide the whole interface" },
  { key: "?", label: "This list" },
];

/** Whether a keypress is the user typing rather than a shortcut. */
function isTyping(target: EventTarget | null): boolean {
  const el = target as HTMLElement | null;

  if (!el) return false;
  const tag = el.tagName;

  return (
    tag === "INPUT" ||
    tag === "TEXTAREA" ||
    tag === "SELECT" ||
    el.isContentEditable
  );
}

export function ShortcutManager() {
  const [listOpen, setListOpen] = useState(false);
  const scLayerVisible = usePanelVisualizationStore((s) => s.scLayerVisible);
  const setScLayerVisible = usePanelVisualizationStore(
    (s) => s.setScLayerVisible,
  );
  const smLayerVisible = usePanelSingleMoleculeVisualizationStore(
    (s) => s.smLayerVisible,
  );
  const setSmLayerVisible = usePanelSingleMoleculeVisualizationStore(
    (s) => s.setSmLayerVisible,
  );
  const layers = usePanelLayers();
  // Molecules are drawn either by a points layer (a store's tiled points) or
  // the legacy molecule layer; both follow the store's smLayerVisible flag.
  const hasMolecules = !!(layers.points ?? layers.molecules);
  const image = layers.image as ImageLayer | undefined;
  const shapes = layers.shapes as ShapesLayer | undefined;

  // What the dataset on screen actually has, so the list can grey out the rest
  const available: Record<string, boolean> = {
    C: true,
    M: hasMolecules,
    I: !!image,
    P: !!shapes,
    H: true,
    "?": true,
  };

  const handler = useCallback(
    (e: KeyboardEvent) => {
      if (e.ctrlKey || e.metaKey || e.altKey) return;
      if (isTyping(e.target)) return;

      if (e.key === "Escape" && listOpen) {
        e.preventDefault();
        setListOpen(false);

        return;
      }
      // `?` is shift+/ on most layouts, so read the character, not the code
      if (e.key === "?") {
        e.preventDefault();
        setListOpen((open) => !open);

        return;
      }

      switch (e.key.toLowerCase()) {
        case "c":
          e.preventDefault();
          setScLayerVisible(!scLayerVisible);
          break;
        case "m":
          if (!hasMolecules) return;
          e.preventDefault();
          setSmLayerVisible(!smLayerVisible);
          break;
        case "i":
          if (!image) return;
          e.preventDefault();
          image.visible = !image.visible;
          break;
        case "p":
          if (!shapes) return;
          e.preventDefault();
          shapes.visible = !shapes.visible;
          break;
        default:
      }
    },
    [
      hasMolecules,
      image,
      listOpen,
      scLayerVisible,
      setScLayerVisible,
      setSmLayerVisible,
      shapes,
      smLayerVisible,
    ],
  );

  useEffect(() => {
    window.addEventListener("keydown", handler);

    return () => window.removeEventListener("keydown", handler);
  }, [handler]);

  if (!listOpen || typeof document === "undefined") return null;

  return createPortal(
    <div
      data-ui-overlay
      className="fixed inset-0 z-[var(--z-modal-top)] flex items-center justify-center bg-black/60 p-6 backdrop-blur-sm"
      role="presentation"
      onClick={(e) => {
        // Only the backdrop closes it, so a click inside stays inside
        if (e.target === e.currentTarget) setListOpen(false);
      }}
    >
      <div
        aria-label="Keyboard shortcuts"
        className={`${glassPanel()} w-full max-w-sm p-5`}
        role="dialog"
      >
        <div className="mb-3 flex items-center justify-between">
          <h2 className="text-sm font-semibold">Keyboard shortcuts</h2>
          <button
            aria-label="Close"
            className="text-white/50 hover:text-white"
            type="button"
            onClick={() => setListOpen(false)}
          >
            ✕
          </button>
        </div>
        <ul className="flex flex-col gap-1.5">
          {SHORTCUTS.map(({ key, label }) => (
            <li
              key={key}
              className={`flex items-center gap-3 text-xs ${
                available[key] ? "text-white/90" : "text-white/35"
              }`}
            >
              <kbd className="min-w-7 rounded border border-white/20 bg-white/10 px-1.5 py-0.5 text-center font-mono">
                {key}
              </kbd>
              <span>{label}</span>
              {!available[key] && (
                <span className="ml-auto text-white/30">not in this view</span>
              )}
            </li>
          ))}
        </ul>
        <p className="mt-3 text-[11px] leading-snug text-white/40">
          Keys act on the main viewer. A split screen&apos;s right panel keeps
          its own controls.
        </p>
      </div>
    </div>,
    document.body,
  );
}
