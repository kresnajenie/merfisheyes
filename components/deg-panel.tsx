"use client";

import type { StandardizedDataset } from "@/lib/StandardizedDataset";

import { Input } from "@heroui/input";
import { Button } from "@heroui/button";
import { Autocomplete, AutocompleteItem } from "@heroui/autocomplete";
import { useMemo, useState } from "react";

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

// Sentinel Autocomplete key for the "vs Rest" reference option. Anything
// starting with "__" can't collide with a real celltype label produced by
// the obs JSON normalizer (which strings + fillna("") values).
const REST_KEY = "__rest__";

export function DegPanel({ onClose: _onClose, controlsRef: _controlsRef }: DegPanelProps) {
  const { getCurrentDataset } = usePanelDatasetStore();
  const {
    selectedColumn,
    selectedGene,
    setSelectedGene,
    deStatsVersion,
    degTarget,
    setDegTarget,
    degReference,
    setDegReference,
    degTargetAuto,
    setDegTargetAuto,
    degReferenceAuto,
    setDegReferenceAuto,
    degSearchTerm,
    setDegSearchTerm,
    degSortKey,
    setDegSortKey,
    degSortDesc,
    setDegSortDesc,
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

  const ranked: RankedDeg[] = useMemo(() => {
    if (!deStats || !degTarget) return [];
    return rankDegsForCelltype(deStats, degTarget, {
      reference: degReference,
    });
  }, [deStats, degTarget, degReference]);

  const sorted = useMemo(() => {
    if (degSortKey === "log2FC" && degSortDesc) return ranked; // already sorted
    const copy = ranked.slice();
    copy.sort((a, b) => {
      const av = a[degSortKey];
      const bv = b[degSortKey];
      return degSortDesc ? bv - av : av - bv;
    });
    return copy;
  }, [ranked, degSortKey, degSortDesc]);

  const filtered = useMemo(() => {
    if (!degSearchTerm.trim()) return sorted;
    const q = degSearchTerm.toLowerCase();
    return sorted.filter((r) => r.gene.toLowerCase().includes(q));
  }, [sorted, degSearchTerm]);

  const targetCount = useMemo(() => {
    if (!deStats || !degTarget) return 0;
    const i = deStats.celltypes.indexOf(degTarget);
    return i === -1 ? 0 : deStats.cellCounts[i];
  }, [deStats, degTarget]);
  const totalCount = useMemo(
    () => deStats?.cellCounts.reduce((a, b) => a + b, 0) ?? 0,
    [deStats],
  );
  const referenceCount = useMemo(() => {
    if (!deStats) return 0;
    if (!degReference) return totalCount - targetCount;
    const i = deStats.celltypes.indexOf(degReference);
    return i === -1 ? 0 : deStats.cellCounts[i];
  }, [deStats, degReference, targetCount, totalCount]);

  const handleSort = (key: SortKey) => {
    if (key === degSortKey) setDegSortDesc(!degSortDesc);
    else {
      setDegSortKey(key);
      setDegSortDesc(true);
    }
  };

  const [copyState, setCopyState] = useState<"idle" | "ok" | "err">("idle");

  const buildCsv = (): string => {
    const header = "gene,log2FC,mean_in,mean_out,pct_in,pct_out,n_in,n_out";
    const lines = [header];
    const nIn = targetCount;
    const nOut = referenceCount;
    for (const r of filtered) {
      lines.push(
        [
          csvCell(r.gene),
          csvNum(r.log2FC),
          csvNum(r.meanIn),
          csvNum(r.meanOut),
          csvNum(r.pctIn),
          csvNum(r.pctOut),
          nIn,
          nOut,
        ].join(","),
      );
    }
    return lines.join("\n");
  };

  const exportFilename = (): string => {
    const stamp = new Date().toISOString().slice(0, 19).replace(/[:T]/g, "-");
    const target = sanitizeForFilename(degTarget ?? "target");
    const ref = sanitizeForFilename(degReference ?? "rest");
    return `deg_${target}_vs_${ref}_${stamp}.csv`;
  };

  const onCopy = async () => {
    if (!filtered.length) return;
    try {
      await navigator.clipboard.writeText(buildCsv());
      setCopyState("ok");
    } catch {
      setCopyState("err");
    }
    setTimeout(() => setCopyState("idle"), 1500);
  };

  const onDownload = () => {
    if (!filtered.length) return;
    const blob = new Blob([buildCsv()], { type: "text/csv;charset=utf-8" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = exportFilename();
    document.body.appendChild(a);
    a.click();
    a.remove();
    URL.revokeObjectURL(url);
  };

  return (
    <div
      className={`absolute top-0 left-16 z-50 w-[420px] border-2 border-white/20 rounded-3xl shadow-lg ${glassButton()}`}
    >
      <div className="p-4 space-y-3">
        {computing && !deStats ? (
          <div className="text-sm text-default-400 py-6 text-center">
            Loading DEG for{" "}
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

            <div className="flex items-end gap-2">
              <Autocomplete
                className="flex-1"
                color="primary"
                isDisabled={degTargetAuto}
                label="Target celltype"
                placeholder="Pick a celltype"
                selectedKey={degTarget ?? undefined}
                onSelectionChange={(key) =>
                  setDegTarget((key as string) || null)
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
              <Button
                className="mb-1 shrink-0"
                color={degTargetAuto ? "primary" : "default"}
                size="sm"
                title={
                  degTargetAuto
                    ? "Target follows the most recently selected celltype"
                    : "Target is set manually"
                }
                variant={degTargetAuto ? "solid" : "flat"}
                onPress={() => setDegTargetAuto(!degTargetAuto)}
              >
                {degTargetAuto ? "Auto" : "Manual"}
              </Button>
            </div>

            <div className="flex items-end gap-2">
              <Autocomplete
                className="flex-1"
                isDisabled={degReferenceAuto}
                label="Reference"
                placeholder="Rest"
                selectedKey={degReference ?? REST_KEY}
                onSelectionChange={(key) => {
                  const k = key as string | null;
                  setDegReference(!k || k === REST_KEY ? null : k);
                }}
              >
                <AutocompleteItem
                  key={REST_KEY}
                  endContent={
                    <span className="text-xs text-default-400">
                      n={(totalCount - targetCount).toLocaleString()}
                    </span>
                  }
                >
                  Rest
                </AutocompleteItem>
                <>
                  {deStats.celltypes
                    .map((ct, i) => ({ ct, i }))
                    .filter(({ ct }) => ct !== degTarget)
                    .map(({ ct, i }) => (
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
                </>
              </Autocomplete>
              <Button
                className="mb-1 shrink-0"
                color={degReferenceAuto ? "primary" : "default"}
                size="sm"
                title={
                  degReferenceAuto
                    ? "Reference follows the 2nd most recently selected celltype (vs Rest when there isn't one)"
                    : "Reference is set manually"
                }
                variant={degReferenceAuto ? "solid" : "flat"}
                onPress={() => setDegReferenceAuto(!degReferenceAuto)}
              >
                {degReferenceAuto ? "2nd selected" : "Manual"}
              </Button>
            </div>

            <div className="text-xs text-default-400 -mt-1">
              {targetCount.toLocaleString()} in target ·{" "}
              {referenceCount.toLocaleString()} in{" "}
              {degReference ?? "rest"}
            </div>

            <div className="flex items-center gap-2">
              <Input
                classNames={{ input: "text-sm" }}
                className="flex-1"
                placeholder="Search gene"
                size="sm"
                value={degSearchTerm}
                onValueChange={setDegSearchTerm}
              />
              <Button
                className="shrink-0"
                color={copyState === "ok" ? "success" : copyState === "err" ? "danger" : "default"}
                isDisabled={!filtered.length}
                size="sm"
                title="Copy filtered rows as CSV to clipboard"
                variant="flat"
                onPress={onCopy}
              >
                {copyState === "ok" ? "Copied" : copyState === "err" ? "Error" : "Copy"}
              </Button>
              <Button
                className="shrink-0"
                isDisabled={!filtered.length}
                size="sm"
                title="Download filtered rows as CSV"
                variant="flat"
                onPress={onDownload}
              >
                CSV
              </Button>
            </div>

            {/* Header */}
            <div className="grid grid-cols-[1.4fr_1fr_0.8fr_0.8fr] gap-2 text-[11px] uppercase tracking-wide text-default-400 px-2">
              <button
                className="text-left hover:text-default-200"
                onClick={() => handleSort("log2FC")}
              >
                Gene
              </button>
              <SortHeader
                active={degSortKey === "log2FC"}
                desc={degSortDesc}
                label="log2FC"
                onClick={() => handleSort("log2FC")}
              />
              <SortHeader
                active={degSortKey === "meanIn"}
                desc={degSortDesc}
                label="mean_in"
                onClick={() => handleSort("meanIn")}
              />
              <SortHeader
                active={degSortKey === "pctIn"}
                desc={degSortDesc}
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

function csvCell(s: string): string {
  if (/[",\n]/.test(s)) return `"${s.replace(/"/g, '""')}"`;
  return s;
}

function csvNum(v: number): string {
  if (Number.isNaN(v)) return "NaN";
  if (v === Infinity) return "Inf";
  if (v === -Infinity) return "-Inf";
  return String(v);
}

function sanitizeForFilename(s: string): string {
  return s.replace(/[^A-Za-z0-9._-]+/g, "_").slice(0, 64);
}
