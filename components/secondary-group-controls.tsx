"use client";

import { useEffect, useMemo, useRef, useState } from "react";

import type { StandardizedDataset } from "@/lib/StandardizedDataset";
import { getColorFromPalette } from "@/lib/utils/color-palette";
import { getEffectiveColumnType } from "@/lib/utils/column-type-utils";
import { loadClusterColumn } from "@/lib/utils/load-cluster-column";
import { getOrderedSecondaryValues } from "@/lib/utils/secondary-order";

interface Props {
  dataset: StandardizedDataset;
  primaryColumn: string | null;
  columnTypeOverrides: Record<string, "categorical" | "numerical">;
  clusterVersion: number;
  incrementClusterVersion: () => void;
  secondaryColumn: string | null;
  selectedSecondaryValues: Set<string>;
  secondaryPaletteOverrides: Record<string, string>;
  secondaryValueOrder: string[];
  setSecondaryColumn: (column: string | null) => void;
  toggleSecondaryValue: (value: string) => void;
  setSecondaryPaletteOverride: (value: string, color: string) => void;
  setSecondaryValueOrder: (values: string[]) => void;
}

export function SecondaryGroupControls({
  dataset,
  primaryColumn,
  columnTypeOverrides,
  clusterVersion,
  incrementClusterVersion,
  secondaryColumn,
  selectedSecondaryValues,
  secondaryPaletteOverrides,
  secondaryValueOrder,
  setSecondaryColumn,
  toggleSecondaryValue,
  setSecondaryPaletteOverride,
  setSecondaryValueOrder,
}: Props) {
  const [loading, setLoading] = useState(false);

  // List of categorical columns excluding the primary.
  const candidates = useMemo(() => {
    if (!dataset) return [] as string[];
    const all =
      dataset.allClusterColumnNames?.length > 0
        ? dataset.allClusterColumnNames
        : (dataset.clusters ?? []).map((c) => c.column);
    return all.filter((name) => {
      if (name === primaryColumn) return false;
      return (
        getEffectiveColumnType(name, dataset, columnTypeOverrides) ===
        "categorical"
      );
    });
  }, [dataset, primaryColumn, columnTypeOverrides, clusterVersion]);

  // Lazy-load the secondary column on selection.
  useEffect(() => {
    if (!secondaryColumn) return;
    const cluster = dataset.clusters?.find((c) => c.column === secondaryColumn);
    if (cluster) return;
    let cancelled = false;
    setLoading(true);
    loadClusterColumn(dataset, secondaryColumn)
      .then((added) => {
        if (cancelled) return;
        if (added) incrementClusterVersion();
      })
      .catch((e) => {
        // eslint-disable-next-line no-console
        console.warn(`Failed to load secondary column ${secondaryColumn}`, e);
      })
      .finally(() => {
        if (!cancelled) setLoading(false);
      });
    return () => {
      cancelled = true;
    };
  }, [dataset, secondaryColumn, incrementClusterVersion]);

  // Unique values for the loaded secondary column, respecting the user's
  // custom ordering when set (falls back to natural order otherwise).
  const values = useMemo(() => {
    if (!secondaryColumn) return [] as string[];
    const cluster = dataset.clusters?.find((c) => c.column === secondaryColumn);
    if (!cluster) return [];
    const natural =
      cluster.uniqueValues && cluster.uniqueValues.length > 0
        ? cluster.uniqueValues
        : Array.from(new Set(cluster.values.map(String)));
    return getOrderedSecondaryValues(natural, secondaryValueOrder);
  }, [dataset, secondaryColumn, clusterVersion, secondaryValueOrder]);

  const moveValue = (idx: number, direction: -1 | 1) => {
    const target = idx + direction;
    if (target < 0 || target >= values.length) return;
    const next = [...values];
    [next[idx], next[target]] = [next[target], next[idx]];
    setSecondaryValueOrder(next);
  };

  return (
    <div
      className="flex items-center gap-2 px-2 py-1 border-b border-white/10"
      onMouseDown={(e) => e.stopPropagation()}
    >
      <label className="text-[10px] text-default-500 flex-shrink-0">
        Group by
      </label>
      <select
        className="text-xs bg-default-100/70 border border-default-300/40 rounded px-1 py-0.5 outline-none focus:border-primary"
        value={secondaryColumn ?? ""}
        onChange={(e) => setSecondaryColumn(e.target.value || null)}
      >
        <option value="">— none —</option>
        {candidates.map((name) => (
          <option key={name} value={name}>
            {name}
          </option>
        ))}
      </select>

      {loading && (
        <span className="text-[10px] text-default-500">loading…</span>
      )}

      {!loading && secondaryColumn && values.length > 0 && (
        <div className="flex items-center gap-1 flex-wrap">
          {values.map((v, i) => (
            <SecondaryValueChip
              key={v}
              canMoveLeft={i > 0}
              canMoveRight={i < values.length - 1}
              color={
                secondaryPaletteOverrides[v] ??
                dataset.clusters?.find((c) => c.column === secondaryColumn)
                  ?.palette?.[v] ??
                getColorFromPalette(i)
              }
              label={v}
              selected={selectedSecondaryValues.has(v)}
              onColorChange={(color) => setSecondaryPaletteOverride(v, color)}
              onMoveLeft={() => moveValue(i, -1)}
              onMoveRight={() => moveValue(i, 1)}
              onToggle={() => toggleSecondaryValue(v)}
            />
          ))}
        </div>
      )}
    </div>
  );
}

function SecondaryValueChip({
  label,
  color,
  selected,
  canMoveLeft,
  canMoveRight,
  onToggle,
  onColorChange,
  onMoveLeft,
  onMoveRight,
}: {
  label: string;
  color: string;
  selected: boolean;
  canMoveLeft: boolean;
  canMoveRight: boolean;
  onToggle: () => void;
  onColorChange: (color: string) => void;
  onMoveLeft: () => void;
  onMoveRight: () => void;
}) {
  const inputRef = useRef<HTMLInputElement>(null);
  return (
    <span
      className={`flex items-center gap-1 text-[10px] rounded-full px-1.5 py-0.5 border ${
        selected
          ? "border-white/30 bg-white/10"
          : "border-white/10 bg-transparent text-default-500"
      }`}
    >
      <button
        aria-label={`Pick colour for ${label}`}
        className="w-3 h-3 rounded-full border border-white/30 flex-shrink-0"
        style={{ backgroundColor: color }}
        type="button"
        onClick={(e) => {
          e.stopPropagation();
          inputRef.current?.click();
        }}
      />
      <input
        ref={inputRef}
        className="hidden"
        type="color"
        value={toHex(color)}
        onChange={(e) => onColorChange(e.target.value)}
      />
      <button
        aria-label={`Move ${label} earlier`}
        className="px-0.5 text-default-400 enabled:hover:text-default-700 disabled:opacity-30 disabled:cursor-not-allowed leading-none"
        disabled={!canMoveLeft}
        type="button"
        onClick={(e) => {
          e.stopPropagation();
          onMoveLeft();
        }}
      >
        ‹
      </button>
      <button
        aria-label={`Move ${label} later`}
        className="px-0.5 text-default-400 enabled:hover:text-default-700 disabled:opacity-30 disabled:cursor-not-allowed leading-none"
        disabled={!canMoveRight}
        type="button"
        onClick={(e) => {
          e.stopPropagation();
          onMoveRight();
        }}
      >
        ›
      </button>
      <button
        className="select-none"
        type="button"
        onClick={(e) => {
          e.stopPropagation();
          onToggle();
        }}
      >
        {label}
      </button>
    </span>
  );
}

// Native color input requires #rrggbb. Best-effort convert rgb()/named colours.
function toHex(color: string): string {
  if (color.startsWith("#")) return color.length === 7 ? color : "#000000";
  const m = color.match(/rgba?\(\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)/);
  if (!m) return "#888888";
  const [, r, g, b] = m;
  const toHex = (x: string) => Number(x).toString(16).padStart(2, "0");
  return `#${toHex(r)}${toHex(g)}${toHex(b)}`;
}
