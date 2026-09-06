"use client";

import type { StandardizedDataset } from "@/lib/StandardizedDataset";
import type { LmMenu } from "@/lib/stores/createLabelledMoleculeVisualizationStore";

import { Button } from "@heroui/button";
import { Checkbox } from "@heroui/checkbox";
import { Input } from "@heroui/input";
import { Slider } from "@heroui/react";
import { Tooltip } from "@heroui/tooltip";
import { useEffect, useMemo } from "react";

import { glassButton, glassPanel } from "@/components/primitives";
import { useSliderRangeLocal } from "@/components/slider-range-popover";
import { LM_MENUS } from "@/lib/stores/createLabelledMoleculeVisualizationStore";
import { useLabelledMoleculeVisualizationStore } from "@/lib/stores/labelledMoleculeVisualizationStore";
import {
  buildSelectionLut,
  countVisible,
} from "@/lib/webgl/labelled-molecule-lut";

const MENU_LABEL: Record<LmMenu, string> = {
  gene: "Gene",
  domain: "Domain",
  cell: "Cell",
};

const MENU_TOOLTIP: Record<LmMenu, string> = {
  gene: "Filter by gene",
  domain: "Filter by RNA domain",
  cell: "Filter by cell",
};

// Same rail geometry as visualization-controls.tsx.
const buttonBaseClass = "w-14 h-14 min-w-0 rounded-full font-medium text-xs";

interface Props {
  dataset: StandardizedDataset;
  clusterVersion?: number;
}

export default function LabelledMoleculeControls({
  dataset,
  clusterVersion = 0,
}: Props) {
  const s = useLabelledMoleculeVisualizationStore();
  // Right-click the size slider to widen its range past the default.
  const sizeRange = useSliderRangeLocal(0.1, 2, s.globalScale);

  const columnFor: Record<LmMenu, string> = {
    gene: "gene",
    domain: "domain",
    cell: "cell",
  };

  const columns = useMemo(() => {
    const out = {} as Record<
      LmMenu,
      {
        uniqueValues: string[];
        valueIndices: ArrayLike<number>;
        palette: Record<string, string> | null;
        counts: number[];
      } | null
    >;

    for (const menu of LM_MENUS) {
      const c = dataset.clusters?.find((cl) => cl.column === columnFor[menu]);

      if (!c?.uniqueValues || !c.valueIndices) {
        out[menu] = null;
        continue;
      }
      // Per-value totals. One pass per column, memoised — not recomputed on
      // every checkbox toggle.
      const counts = new Array(c.uniqueValues.length).fill(0);

      for (let i = 0; i < c.valueIndices.length; i++)
        counts[c.valueIndices[i]]++;

      out[menu] = {
        uniqueValues: c.uniqueValues,
        valueIndices: c.valueIndices,
        palette: c.palette,
        counts,
      };
    }

    return out;
  }, [dataset, clusterVersion]);

  const visible = useMemo(() => {
    if (LM_MENUS.some((m) => !columns[m])) return null;

    const cols = LM_MENUS.map((m) => ({
      indices: columns[m]!.valueIndices,
      lut: buildSelectionLut(
        columns[m]!.uniqueValues,
        s.selections[m],
        s.hiddenValues[m],
      ),
    }));

    return countVisible(cols, dataset.getPointCount());
  }, [columns, s.selections, s.hiddenValues, dataset]);

  const open = s.openMenu;
  const openCol = open ? columns[open] : null;

  const filtered = useMemo(() => {
    if (!open || !openCol) return [];
    const term = s.searchTerm[open].toLowerCase();

    return openCol.uniqueValues
      .map((_, i) => i)
      .filter(
        (i) => !term || openCol.uniqueValues[i].toLowerCase().includes(term),
      )
      .sort((a, b) => {
        // Only the gene list is sortable; the short lists read best A→Z.
        if (open !== "gene") {
          return openCol.uniqueValues[a].localeCompare(openCol.uniqueValues[b]);
        }
        const sign = s.geneSortDir === "asc" ? 1 : -1;

        return s.geneSortBy === "count"
          ? sign * (openCol.counts[a] - openCol.counts[b])
          : sign *
              openCol.uniqueValues[a].localeCompare(openCol.uniqueValues[b]);
      });
  }, [open, openCol, s.searchTerm, s.geneSortBy, s.geneSortDir]);

  // G / D / C switch the colouring column. Ignored while typing (the value
  // search box is a text input) and when a modifier is held, matching how
  // HideUiManager guards its own H shortcut.
  const setColorBy = s.setColorBy;

  useEffect(() => {
    const keys: Record<string, LmMenu> = { g: "gene", d: "domain", c: "cell" };
    const handler = (e: KeyboardEvent) => {
      if (e.ctrlKey || e.metaKey || e.altKey) return;
      const menu = keys[e.key.toLowerCase()];

      if (!menu) return;
      const el = e.target as HTMLElement | null;

      if (
        el &&
        (el.tagName === "INPUT" ||
          el.tagName === "TEXTAREA" ||
          el.isContentEditable)
      ) {
        return;
      }
      setColorBy(menu);
    };

    window.addEventListener("keydown", handler);

    return () => window.removeEventListener("keydown", handler);
  }, [setColorBy]);

  const openPanel = (menu: LmMenu) => {
    s.setOpenMenu(open === menu ? null : menu);
  };

  return (
    <>
      <div
        data-ui-overlay
        className="absolute top-28 left-4 z-[var(--z-rail)] flex flex-col gap-2"
      >
        {LM_MENUS.map((menu) => {
          const isOpen = open === menu;
          const isColoring = s.colorBy === menu;
          const n = s.selections[menu].size;

          return (
            <Tooltip
              key={menu}
              content={
                `${MENU_TOOLTIP[menu]}${isColoring ? " · colouring the scene" : ""}` +
                (n ? ` · ${n} selected` : "") +
                ` · ${MENU_LABEL[menu][0]} to colour by this`
              }
              placement="right"
            >
              <Button
                className={`${buttonBaseClass} ${isOpen ? "" : glassButton()}`}
                color={isOpen ? "primary" : "default"}
                variant={isOpen ? "shadow" : "light"}
                onPress={() => openPanel(menu)}
              >
                {MENU_LABEL[menu]}
              </Button>
            </Tooltip>
          );
        })}

        <Tooltip
          content="Molecule size (right-click to edit range)"
          placement="right"
        >
          <div
            className={`w-14 h-32 rounded-full border-2 border-default-200 p-2 flex flex-col items-center justify-center ${glassButton()}`}
            onContextMenu={sizeRange.onContextMenu}
          >
            <Slider
              aria-label="Molecule size"
              className="h-full"
              maxValue={sizeRange.max}
              minValue={sizeRange.min}
              orientation="vertical"
              size="sm"
              step={0.05}
              value={s.globalScale}
              onChange={(v) => s.setGlobalScale(v as number)}
            />
            {sizeRange.popover}
          </div>
        </Tooltip>

        {/* Value panel — same placement and shell as VisualizationPanel. */}
        {open && openCol && (
          <div
            className={`absolute top-0 left-16 z-[var(--z-panel)] w-[300px] ${glassPanel()}`}
          >
            <div className="p-4 space-y-3">
              <div className="flex items-center justify-between">
                <span className="text-sm font-medium">{MENU_LABEL[open]}</span>
                <Button
                  color={s.colorBy === open ? "primary" : "default"}
                  size="sm"
                  variant={s.colorBy === open ? "flat" : "light"}
                  onPress={() => s.setColorBy(open)}
                >
                  {s.colorBy === open ? "Colouring" : "Colour by"}
                </Button>
              </div>

              {/* Genes are numerous enough to need ordering; domain and cell
                  are short lists and stay alphabetical. */}
              {open === "gene" && (
                <div className="flex gap-1">
                  <Button
                    className="flex-1 text-xs"
                    color={s.geneSortBy === "alpha" ? "primary" : "default"}
                    size="sm"
                    variant={s.geneSortBy === "alpha" ? "flat" : "light"}
                    onPress={() =>
                      s.setGeneSort(
                        "alpha",
                        s.geneSortBy === "alpha" && s.geneSortDir === "asc"
                          ? "desc"
                          : "asc",
                      )
                    }
                  >
                    A→Z{" "}
                    {s.geneSortBy === "alpha"
                      ? s.geneSortDir === "asc"
                        ? "↑"
                        : "↓"
                      : ""}
                  </Button>
                  <Button
                    className="flex-1 text-xs"
                    color={s.geneSortBy === "count" ? "primary" : "default"}
                    size="sm"
                    variant={s.geneSortBy === "count" ? "flat" : "light"}
                    onPress={() =>
                      s.setGeneSort(
                        "count",
                        s.geneSortBy === "count" && s.geneSortDir === "desc"
                          ? "asc"
                          : "desc",
                      )
                    }
                  >
                    # Count{" "}
                    {s.geneSortBy === "count"
                      ? s.geneSortDir === "asc"
                        ? "↑"
                        : "↓"
                      : ""}
                  </Button>
                </div>
              )}

              <Input
                classNames={{ input: "text-sm" }}
                placeholder={`Search ${MENU_LABEL[open].toLowerCase()}`}
                value={s.searchTerm[open]}
                onValueChange={(v) => s.setSearchTerm(open, v)}
              />

              <div className="flex gap-1">
                <Button
                  className="flex-1"
                  color="danger"
                  variant="ghost"
                  onPress={() => {
                    s.setSearchTerm(open, "");
                    s.clearMenu(open);
                  }}
                >
                  Clear
                </Button>
              </div>

              <div className="text-xs text-default-500 text-center">
                {s.selections[open].size === 0
                  ? "No filter — showing all"
                  : `${s.selections[open].size} of ${openCol.uniqueValues.length} selected`}
              </div>

              <div className="max-h-[400px] overflow-y-auto flex flex-col gap-0">
                {filtered.map((i) => {
                  const value = openCol.uniqueValues[i];

                  return (
                    <Checkbox
                      key={value}
                      classNames={{
                        base: "w-full max-w-full",
                        label: "w-full min-w-0",
                      }}
                      isSelected={s.selections[open].has(value)}
                      size="sm"
                      onValueChange={() => s.toggleValue(open, value)}
                    >
                      <div className="flex w-full min-w-0 items-center gap-2">
                        <span className="min-w-0 flex-1 truncate" title={value}>
                          {open === "gene" ? value.split(" (")[0] : value}
                        </span>
                        {/* Fixed-width column so counts line up down the list
                            rather than trailing each name. */}
                        <span className="w-16 shrink-0 text-right text-[10px] tabular-nums text-default-500">
                          {openCol.counts[i].toLocaleString()}
                        </span>
                      </div>
                    </Checkbox>
                  );
                })}
              </div>
            </div>
          </div>
        )}
      </div>

      {/* How many molecules survive the intersection. */}
      <div
        data-ui-overlay
        className={`absolute bottom-6 left-4 z-[var(--z-legends)] rounded-full px-4 py-2 text-xs ${glassButton()}`}
        data-testid="lm-visible-count"
      >
        <span className="font-medium">
          {visible === null ? "…" : visible.toLocaleString()}
        </span>
        <span className="text-default-500">
          {" "}
          / {dataset.getPointCount().toLocaleString()} molecules
        </span>
      </div>
    </>
  );
}
