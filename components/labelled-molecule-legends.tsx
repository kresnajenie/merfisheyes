"use client";

import type { StandardizedDataset } from "@/lib/StandardizedDataset";
import type { LmMenu } from "@/lib/stores/createLabelledMoleculeVisualizationStore";

import { Button } from "@heroui/button";
import { Popover, PopoverContent, PopoverTrigger } from "@heroui/popover";
import { Slider } from "@heroui/react";
import { Eye, EyeOff, X } from "lucide-react";
import { useMemo, useState } from "react";

import { glassPanel } from "@/components/primitives";
import {
  ColorPicker,
  ColorPickerEyeDropper,
  ColorPickerHue,
  ColorPickerOutput,
  ColorPickerSelection,
} from "@/components/ui/shadcn-io/color-picker";
import {
  LM_MENUS,
  resolveValueColor,
} from "@/lib/stores/createLabelledMoleculeVisualizationStore";
import { useLabelledMoleculeVisualizationStore } from "@/lib/stores/labelledMoleculeVisualizationStore";

const MENU_LABEL: Record<LmMenu, string> = {
  gene: "Gene",
  domain: "Domain",
  cell: "Cell",
};

interface Props {
  dataset: StandardizedDataset;
  clusterVersion?: number;
}

/**
 * Active selections, one group per menu, plus how unselected molecules draw.
 * Chips follow the single-cell legend: an eye to hide without deselecting, the
 * label to recolour, an X to remove.
 */
export default function LabelledMoleculeLegends({
  dataset,
  clusterVersion = 0,
}: Props) {
  const s = useLabelledMoleculeVisualizationStore();
  const [openPicker, setOpenPicker] = useState<string | null>(null);

  const palettes = useMemo(() => {
    return LM_MENUS.reduce(
      (acc, menu) => {
        acc[menu] =
          dataset.clusters?.find((c) => c.column === menu)?.palette ?? null;

        return acc;
      },
      {} as Record<LmMenu, Record<string, string> | null>,
    );
  }, [dataset, clusterVersion]);

  const active = LM_MENUS.filter((m) => s.selections[m].size > 0);

  return (
    <div
      data-ui-overlay
      className="absolute right-6 top-24 z-[var(--z-legends)] flex max-w-xs flex-col items-end gap-3"
    >
      {/* Which column drives colour. Always shown, so the scene is never
          ambiguous — and clicking a name switches it. */}
      <div className={`w-full p-3 ${glassPanel()}`}>
        <div className="mb-2 text-[11px] uppercase tracking-wide text-default-500">
          Colouring by
        </div>
        <div className="flex gap-1">
          {LM_MENUS.map((menu) => (
            <Button
              key={menu}
              className="flex-1"
              color={s.colorBy === menu ? "primary" : "default"}
              size="sm"
              variant={s.colorBy === menu ? "flat" : "light"}
              onPress={() => s.setColorBy(menu)}
            >
              {MENU_LABEL[menu]}
            </Button>
          ))}
        </div>
      </div>

      {active.map((menu) => (
        <div key={menu} className={`w-full p-3 ${glassPanel()}`}>
          <div className="mb-2 flex items-center justify-between">
            <span className="text-xs text-default-500">{MENU_LABEL[menu]}</span>
            <button
              className="text-[11px] text-default-500 hover:text-danger"
              type="button"
              onClick={() => s.clearMenu(menu)}
            >
              Clear
            </button>
          </div>

          <div className="flex flex-wrap justify-end gap-1.5">
            {[...s.selections[menu]].map((value) => {
              const color = resolveValueColor(
                menu,
                value,
                s.colorOverrides[menu],
                s.geneColorSlots,
                palettes[menu],
              );
              const hidden = s.hiddenValues[menu].has(value);
              const key = `${menu}::${value}`;

              return (
                <div
                  key={value}
                  className="group flex items-center gap-1.5 rounded-full bg-default-100/80 px-2 py-1 text-xs transition-all hover:scale-105"
                  style={{ opacity: hidden ? 0.45 : 1 }}
                >
                  <button
                    className="flex items-center justify-center rounded-full p-0.5 text-default-500 transition-colors hover:bg-default-200 hover:text-foreground"
                    title={
                      hidden
                        ? "Hidden — click to show"
                        : "Visible — click to hide"
                    }
                    type="button"
                    onClick={() => s.toggleValueVisibility(menu, value)}
                  >
                    {hidden ? (
                      <EyeOff className="h-3 w-3" />
                    ) : (
                      <Eye className="h-3 w-3" />
                    )}
                  </button>

                  <Popover
                    isOpen={openPicker === key}
                    placement="left"
                    onOpenChange={(o) => !o && setOpenPicker(null)}
                  >
                    <PopoverTrigger>
                      <button
                        className="flex cursor-pointer items-center gap-1.5"
                        title={`${value} — click to recolour`}
                        type="button"
                        onClick={() => setOpenPicker(key)}
                      >
                        <span
                          className="h-2.5 w-2.5 shrink-0 rounded-full ring-1 ring-black/20"
                          style={{ backgroundColor: color }}
                        />
                        <span className="max-w-[130px] truncate">
                          {menu === "gene" ? value.split(" (")[0] : value}
                        </span>
                      </button>
                    </PopoverTrigger>
                    <PopoverContent className="p-2">
                      <ColorPicker
                        key={key}
                        className="rounded-md p-2"
                        value={color}
                        onChange={(v) => s.setColorOverride(menu, value, v)}
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
                    className="flex items-center justify-center rounded-full p-0.5 text-default-500 transition-colors hover:bg-default-200 hover:text-danger"
                    title="Remove"
                    type="button"
                    onClick={() => s.toggleValue(menu, value)}
                  >
                    <X className="h-3 w-3" />
                  </button>
                </div>
              );
            })}
          </div>
        </div>
      ))}

      {active.length > 0 && (
        <div className={`w-full p-3 ${glassPanel()}`}>
          <div className="mb-2 text-xs text-default-500">Unselected</div>
          <div className="flex gap-1">
            {(
              [
                ["hidden", "Hidden"],
                ["grey", "Transparent"],
              ] as const
            ).map(([mode, label]) => (
              <Button
                key={mode}
                className="flex-1"
                color={s.unselectedMode === mode ? "primary" : "default"}
                size="sm"
                variant={s.unselectedMode === mode ? "flat" : "light"}
                onPress={() => s.setUnselectedMode(mode)}
              >
                {label}
              </Button>
            ))}
          </div>

          {s.unselectedMode === "grey" && (
            <Slider
              aria-label="Unselected opacity"
              className="mt-3"
              label="Opacity"
              maxValue={1}
              minValue={0}
              size="sm"
              step={0.02}
              value={s.unselectedAlpha}
              onChange={(v) => s.setUnselectedAlpha(Number(v))}
            />
          )}
        </div>
      )}
    </div>
  );
}
