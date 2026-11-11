"use client";

import { Popover, PopoverContent, PopoverTrigger } from "@heroui/react";
import { useState } from "react";

import {
  ColorPicker,
  ColorPickerEyeDropper,
  ColorPickerHue,
  ColorPickerOutput,
  ColorPickerSelection,
} from "@/components/ui/shadcn-io/color-picker";
import { shuffleColors, spreadColors } from "@/lib/utils/color-assign";

/**
 * What to do with several legend chips picked at once (useLegendSelection):
 * give them all one colour, or distinct colours again — spread, which keeps
 * away from the colours around them, or shuffle, for when a spread lands on
 * two that look alike.
 *
 * Shown by both legends; each passes `apply`, which writes one colour where
 * that legend keeps it (a cluster's palette, a gene's entry).
 */
export function LegendColorActions({
  selected,
  taken,
  apply,
  onDone,
}: {
  /** The chips picked, in the order they are drawn. */
  selected: string[];
  /** Every colour in use in this legend, so a spread can avoid them. */
  taken: Iterable<string>;
  /** Write one item's colour. Called once per item. */
  apply: (key: string, color: string) => void;
  /** Called after the colours are written, to drop the selection. */
  onDone: () => void;
}) {
  const [pickerOpen, setPickerOpen] = useState(false);
  const [draft, setDraft] = useState("#ffffff");

  if (selected.length < 2) return null;

  const applyAll = (colors: Record<string, string>) => {
    for (const [key, color] of Object.entries(colors)) apply(key, color);
  };

  return (
    <div
      className="flex items-center gap-1.5 rounded-full bg-black/50 px-2 py-1 text-[11px] text-white/80 backdrop-blur-md"
      role="toolbar"
    >
      <span className="px-1 tabular-nums">{selected.length} picked</span>
      <Popover
        isOpen={pickerOpen}
        placement="left"
        onOpenChange={setPickerOpen}
      >
        <PopoverTrigger>
          <button
            className="rounded-full px-2 py-0.5 hover:bg-white/15"
            title="Give all of them one colour"
            type="button"
          >
            One colour
          </button>
        </PopoverTrigger>
        <PopoverContent className="glass-panel w-64 p-3">
          <ColorPicker
            className="rounded-md p-2"
            value={draft}
            onChange={(value) => {
              setDraft(value);
              applyAll(Object.fromEntries(selected.map((key) => [key, value])));
            }}
          >
            <ColorPickerSelection />
            <div className="mt-2 flex items-center gap-2">
              <ColorPickerEyeDropper />
              <div className="grid w-full gap-1">
                <ColorPickerHue />
              </div>
            </div>
            <div className="mt-2 text-center">
              <ColorPickerOutput />
            </div>
          </ColorPicker>
        </PopoverContent>
      </Popover>
      <button
        className="rounded-full px-2 py-0.5 hover:bg-white/15"
        title="Give each of them a different colour"
        type="button"
        onClick={() => {
          applyAll(spreadColors(selected, taken));
          onDone();
        }}
      >
        Spread
      </button>
      <button
        className="rounded-full px-2 py-0.5 hover:bg-white/15"
        title="Try another set of different colours"
        type="button"
        onClick={() => applyAll(shuffleColors(selected))}
      >
        Shuffle
      </button>
      <button
        className="rounded-full px-1.5 py-0.5 text-white/50 hover:bg-white/15 hover:text-white"
        title="Clear the selection"
        type="button"
        onClick={onDone}
      >
        ✕
      </button>
    </div>
  );
}
