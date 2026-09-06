"use client";

import type { StandardizedDataset } from "@/lib/StandardizedDataset";
import type { LmMenu } from "@/lib/stores/createLabelledMoleculeVisualizationStore";

import { Button } from "@heroui/button";
import { Slider } from "@heroui/react";
import { X } from "lucide-react";
import { useMemo } from "react";

import { glassPanel } from "@/components/primitives";
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
 * Active selections, one group per menu, plus the controls for how the
 * molecules that fail the filter are drawn.
 *
 * Auto-hides when nothing is selected — at which point every molecule passes,
 * so there is no "unselected" appearance to configure either.
 */
export default function LabelledMoleculeLegends({
  dataset,
  clusterVersion = 0,
}: Props) {
  const s = useLabelledMoleculeVisualizationStore();

  const palettes = useMemo(() => {
    const columnFor: Record<LmMenu, string> = {
      gene: "gene",
      domain: "domain",
      cell: "cell",
    };

    return LM_MENUS.reduce(
      (acc, menu) => {
        acc[menu] =
          dataset.clusters?.find((c) => c.column === columnFor[menu])
            ?.palette ?? null;

        return acc;
      },
      {} as Record<LmMenu, Record<string, string> | null>,
    );
  }, [dataset, clusterVersion]);

  const active = LM_MENUS.filter((m) => s.selections[m].size > 0);

  if (active.length === 0) return null;

  return (
    <div
      data-ui-overlay
      className="absolute right-6 top-24 z-[var(--z-legends)] flex max-w-xs flex-col items-end gap-3"
    >
      {active.map((menu) => (
        <div key={menu} className={`w-full p-3 ${glassPanel()}`}>
          <div className="mb-2 flex items-center justify-between">
            <button
              className={`text-xs ${
                s.colorBy === menu
                  ? "font-medium text-foreground"
                  : "text-default-500 hover:text-foreground"
              }`}
              title={
                s.colorBy === menu
                  ? "Colouring the scene"
                  : "Click to colour the scene by this"
              }
              type="button"
              onClick={() => s.setColorBy(menu)}
            >
              {MENU_LABEL[menu]}
              {s.colorBy === menu ? " ·" : ""}
            </button>
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

              return (
                <div
                  key={value}
                  className="group flex cursor-default items-center gap-1.5 rounded-full px-3 py-1 text-xs transition-all hover:scale-105"
                  style={{ backgroundColor: color }}
                  title={value}
                >
                  <span className="max-w-[140px] truncate text-black">
                    {menu === "gene" ? value.split(" (")[0] : value}
                  </span>
                  <button
                    className="rounded-full p-0.5 text-black/60 transition-colors hover:bg-black/15 hover:text-black"
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

      {/* How the molecules that fail the filter are drawn. */}
      <div className={`w-full p-3 ${glassPanel()}`}>
        <div className="mb-2 text-xs text-default-500">Unselected</div>
        <div className="flex gap-1">
          {(
            [
              ["grey", "Transparent"],
              ["hidden", "Hidden"],
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
    </div>
  );
}
