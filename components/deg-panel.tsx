"use client";

import type { StandardizedDataset } from "@/lib/StandardizedDataset";

import { Input } from "@heroui/input";
import { Autocomplete, AutocompleteItem } from "@heroui/autocomplete";
import { useEffect, useMemo, useRef, useState } from "react";

import {
  usePanelDatasetStore,
  usePanelVisualizationStore,
} from "@/lib/hooks/usePanelStores";
import { glassButton } from "@/components/primitives";
import {
  rankDegsForCelltype,
  isDeStatsInFlight,
  type RankedDeg,
} from "@/lib/utils/de-stats";

interface DegPanelProps {
  onClose: () => void;
  controlsRef?: React.RefObject<HTMLDivElement>;
}

type SortKey = "log2FC" | "meanIn" | "pctIn";

export function DegPanel({ onClose, controlsRef }: DegPanelProps) {
  const panelRef = useRef<HTMLDivElement>(null);
  const { getCurrentDataset } = usePanelDatasetStore();
  const {
    selectedColumn,
    selectedCelltypes,
    selectedGene,
    setSelectedGene,
    deStatsVersion,
  } = usePanelVisualizationStore();

  const rawDataset = getCurrentDataset();
  const dataset =
    rawDataset && "clusters" in rawDataset
      ? (rawDataset as StandardizedDataset)
      : null;

  // Prefer the column the user currently has selected in the celltype panel,
  // fall back to whatever column was eagerly computed at parse time.
  const activeColumn =
    selectedColumn && dataset?.deStatsByColumn.has(selectedColumn)
      ? selectedColumn
      : dataset?.deStats?.column ?? null;
  const deStats = activeColumn
    ? dataset?.deStatsByColumn.get(activeColumn) ?? null
    : null;
  const computing = !!(
    dataset &&
    selectedColumn &&
    !dataset.deStatsByColumn.has(selectedColumn) &&
    isDeStatsInFlight(dataset.id, selectedColumn)
  );
  // Touch deStatsVersion so this component re-renders when a compute finishes.
  void deStatsVersion;

  // Default target celltype: first currently selected, else first in deStats.
  const initialTarget = useMemo(() => {
    if (!deStats || deStats.celltypes.length === 0) return null;
    const cts = new Set(deStats.celltypes);
    for (const ct of selectedCelltypes) {
      if (cts.has(ct)) return ct;
    }
    return deStats.celltypes[0];
  }, [deStats, selectedCelltypes]);

  const [target, setTarget] = useState<string | null>(initialTarget);
  const [searchTerm, setSearchTerm] = useState("");
  const [sortKey, setSortKey] = useState<SortKey>("log2FC");
  const [sortDesc, setSortDesc] = useState(true);

  // Reset target whenever the active column changes (celltypes are a different
  // set per column).
  useEffect(() => {
    setTarget(initialTarget);
  }, [activeColumn, initialTarget]);

  // Click outside to close
  useEffect(() => {
    const handleClickOutside = (event: MouseEvent) => {
      const t = event.target as Node;
      if (panelRef.current && panelRef.current.contains(t)) return;
      if (controlsRef?.current && controlsRef.current.contains(t)) return;
      if (
        t instanceof HTMLElement &&
        t.closest(
          '[role="listbox"], [data-slot="content"], [role="dialog"], [data-slot="backdrop"]',
        )
      )
        return;
      onClose();
    };
    document.addEventListener("mousedown", handleClickOutside);
    return () => document.removeEventListener("mousedown", handleClickOutside);
  }, [onClose, controlsRef]);

  const ranked: RankedDeg[] = useMemo(() => {
    if (!deStats || !target) return [];
    return rankDegsForCelltype(deStats, target);
  }, [deStats, target]);

  const sorted = useMemo(() => {
    if (sortKey === "log2FC" && sortDesc) return ranked; // already sorted
    const copy = ranked.slice();
    copy.sort((a, b) => {
      const av = a[sortKey];
      const bv = b[sortKey];
      return sortDesc ? bv - av : av - bv;
    });
    return copy;
  }, [ranked, sortKey, sortDesc]);

  const filtered = useMemo(() => {
    if (!searchTerm.trim()) return sorted;
    const q = searchTerm.toLowerCase();
    return sorted.filter((r) => r.gene.toLowerCase().includes(q));
  }, [sorted, searchTerm]);

  const targetCount = useMemo(() => {
    if (!deStats || !target) return 0;
    const i = deStats.celltypes.indexOf(target);
    return i === -1 ? 0 : deStats.cellCounts[i];
  }, [deStats, target]);
  const totalCount = useMemo(
    () => deStats?.cellCounts.reduce((a, b) => a + b, 0) ?? 0,
    [deStats],
  );

  const handleSort = (key: SortKey) => {
    if (key === sortKey) setSortDesc((d) => !d);
    else {
      setSortKey(key);
      setSortDesc(true);
    }
  };

  return (
    <div
      ref={panelRef}
      className={`absolute top-0 left-16 z-50 w-[420px] border-2 border-white/20 rounded-3xl shadow-lg ${glassButton()}`}
    >
      <div className="p-4 space-y-3">
        {computing && !deStats ? (
          <div className="text-sm text-default-400 py-6 text-center">
            Computing DEG for{" "}
            <span className="text-default-200 font-medium">
              {selectedColumn}
            </span>
            …
          </div>
        ) : !deStats ? (
          <div className="text-sm text-default-400 py-4 text-center">
            DEG stats are not available for this dataset. Drop an H5AD file to
            precompute per-celltype expression.
          </div>
        ) : (
          <>
            <div className="text-xs text-default-400 flex items-center gap-2">
              Column:{" "}
              <span className="text-default-200 font-medium">
                {deStats.column}
              </span>
              {computing && (
                <span className="text-default-500 italic">
                  · recomputing for "{selectedColumn}"…
                </span>
              )}
            </div>

            <Autocomplete
              className="max-w-full"
              color="primary"
              label="Target celltype"
              placeholder="Pick a celltype"
              selectedKey={target ?? undefined}
              onSelectionChange={(key) =>
                setTarget((key as string) || null)
              }
            >
              {deStats.celltypes.map((ct, i) => (
                <AutocompleteItem
                  key={ct}
                  endContent={
                    <span className="text-xs text-default-400">
                      n={deStats.cellCounts[i]}
                    </span>
                  }
                >
                  {ct}
                </AutocompleteItem>
              ))}
            </Autocomplete>

            <div className="text-xs text-default-400 -mt-1">
              {targetCount.toLocaleString()} cells in target ·{" "}
              {(totalCount - targetCount).toLocaleString()} in rest
            </div>

            <Input
              classNames={{ input: "text-sm" }}
              placeholder="Search gene"
              size="sm"
              value={searchTerm}
              onValueChange={setSearchTerm}
            />

            {/* Header */}
            <div className="grid grid-cols-[1.4fr_1fr_0.8fr_0.8fr] gap-2 text-[11px] uppercase tracking-wide text-default-400 px-2">
              <button
                className="text-left hover:text-default-200"
                onClick={() => handleSort("log2FC")}
              >
                Gene
              </button>
              <SortHeader
                active={sortKey === "log2FC"}
                desc={sortDesc}
                label="log2FC"
                onClick={() => handleSort("log2FC")}
              />
              <SortHeader
                active={sortKey === "meanIn"}
                desc={sortDesc}
                label="mean_in"
                onClick={() => handleSort("meanIn")}
              />
              <SortHeader
                active={sortKey === "pctIn"}
                desc={sortDesc}
                label="pct_in"
                onClick={() => handleSort("pctIn")}
              />
            </div>

            {/* Rows */}
            <div className="max-h-[420px] overflow-y-auto flex flex-col">
              {filtered.length === 0 ? (
                <div className="text-xs text-default-400 py-4 text-center">
                  No genes match.
                </div>
              ) : (
                filtered.slice(0, 500).map((r) => {
                  const isSelected = r.gene === selectedGene;
                  return (
                    <button
                      key={r.gene}
                      className={`grid grid-cols-[1.4fr_1fr_0.8fr_0.8fr] gap-2 px-2 py-1.5 rounded text-left text-xs hover:bg-white/10 transition-colors ${
                        isSelected ? "bg-primary/30" : ""
                      }`}
                      onClick={() => setSelectedGene(r.gene)}
                      title={`mean_out=${r.meanOut.toFixed(4)} · pct_out=${r.pctOut.toFixed(3)}`}
                    >
                      <span className="font-mono truncate">{r.gene}</span>
                      <span
                        className={`font-mono ${
                          r.log2FC > 0
                            ? "text-rose-300"
                            : r.log2FC < 0
                              ? "text-sky-300"
                              : "text-default-300"
                        }`}
                      >
                        {formatLog2FC(r.log2FC)}
                      </span>
                      <span className="font-mono text-default-300">
                        {r.meanIn.toFixed(3)}
                      </span>
                      <span className="font-mono text-default-300">
                        {(r.pctIn * 100).toFixed(1)}%
                      </span>
                    </button>
                  );
                })
              )}
              {filtered.length > 500 && (
                <div className="text-[10px] text-default-500 text-center py-1">
                  Showing top 500 of {filtered.length}
                </div>
              )}
            </div>
          </>
        )}
      </div>
    </div>
  );
}

function SortHeader({
  active,
  desc,
  label,
  onClick,
}: {
  active: boolean;
  desc: boolean;
  label: string;
  onClick: () => void;
}) {
  return (
    <button
      className={`text-left ${active ? "text-default-200" : "hover:text-default-200"}`}
      onClick={onClick}
    >
      {label}
      {active ? (desc ? " ↓" : " ↑") : ""}
    </button>
  );
}

function formatLog2FC(v: number): string {
  if (!isFinite(v)) return v > 0 ? "+∞" : "-∞";
  const sign = v > 0 ? "+" : "";
  return sign + v.toFixed(2);
}
