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
// The overlay is rendered by three-scene from the GLOBAL single-molecule stores
// (it is not panel-scoped), so read/write the same globals here to stay in sync.
import { useSingleMoleculeStore } from "@/lib/stores/singleMoleculeStore";
import { useSingleMoleculeVisualizationStore } from "@/lib/stores/singleMoleculeVisualizationStore";
import { glassPanel } from "@/components/primitives";
import {
  rankDegsForCelltype,
  isDeStatsInFlight,
  type RankedDeg,
} from "@/lib/utils/de-stats";

interface DegPanelProps {
  onClose: () => void;
  controlsRef?: React.RefObject<HTMLDivElement>;
}

type SortKey = "fc" | "mean1" | "pct";

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
    degDataMode,
    setDegDataMode,
  } = usePanelVisualizationStore();

  // Single-molecule overlay (present only when this SC dataset has one linked).
  // When present, gene clicks can be routed to the molecule overlay instead of
  // the single-cell expression colouring.
  const smDataset = useSingleMoleculeStore((s) => s.getCurrentDataset());
  const smSelectedGenes = useSingleMoleculeVisualizationStore(
    (s) => s.selectedGenes,
  );
  const smAddGene = useSingleMoleculeVisualizationStore((s) => s.addGene);
  const smRemoveGene = useSingleMoleculeVisualizationStore((s) => s.removeGene);

  const [showExprInfo, setShowExprInfo] = useState(false);
  const hasSmOverlay = !!smDataset;
  const [clickTargetRaw, setClickTarget] = useState<"cell" | "molecule">("cell");
  // Fall back to "cell" whenever there's no overlay so a stale "molecule"
  // choice can't silently swallow clicks.
  const clickTarget = hasSmOverlay ? clickTargetRaw : "cell";
  const smGeneSet = useMemo(
    () => (smDataset ? new Set(smDataset.uniqueGenes) : null),
    [smDataset],
  );

  const onGeneClick = (gene: string) => {
    if (clickTarget === "molecule") {
      if (smSelectedGenes.has(gene)) {
        smRemoveGene(gene);
      } else if (!smGeneSet || smGeneSet.has(gene)) {
        smAddGene(gene);
      }

      return;
    }
    setSelectedGene(gene);
  };

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

  // Detect raw counts vs log-normalized from the stored means. Log-normalized
  // matrices (log1p of library-size-normalized counts) rarely exceed ~15;
  // raw-count means of well-expressed genes run much higher. The user can
  // override this in the panel.
  const detectedMode: "counts" | "log" = useMemo(() => {
    if (!deStats) return "counts";
    let max = 0;
    const m = deStats.means;
    for (let i = 0; i < m.length; i++) if (m[i] > max) max = m[i];
    return max > 30 ? "counts" : "log";
  }, [deStats]);
  const mode: "counts" | "log" =
    degDataMode === "auto" ? detectedMode : degDataMode;

  // Per-gene display values: for log data, un-log the means so mean1/mean2 and
  // the fold change are in linear (interpretable) units.
  const display = useMemo(() => {
    const lin = (v: number) => (mode === "log" ? Math.expm1(v) : v);
    return ranked.map((r) => {
      const mean1 = lin(r.meanIn);
      const mean2 = lin(r.meanOut);
      const fc = mean2 > 0 ? mean1 / mean2 : mean1 > 0 ? Infinity : 1;
      return { ...r, mean1, mean2, fc };
    });
  }, [ranked, mode]);

  const sorted = useMemo(() => {
    const key =
      degSortKey === "mean1" ? "mean1" : degSortKey === "pct" ? "pctIn" : "fc";
    const copy = display.slice();
    copy.sort((a, b) => {
      const av = a[key as keyof typeof a] as number;
      const bv = b[key as keyof typeof b] as number;
      return degSortDesc ? bv - av : av - bv;
    });
    return copy;
  }, [display, degSortKey, degSortDesc]);

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
    const header = `# expression: ${mode === "log" ? "log-normalized (means un-logged for FC)" : "raw counts"}\ngene,mean1,mean2,fold_change,pct_cells_grp1,pct_cells_grp2,n_grp1,n_grp2`;
    const lines = [header];
    const nIn = targetCount;
    const nOut = referenceCount;
    for (const r of filtered) {
      lines.push(
        [
          csvCell(r.gene),
          csvNum(r.mean1),
          csvNum(r.mean2),
          csvNum(r.fc),
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
      data-testid="deg-panel"
      data-deg-column={deStats?.column ?? ""}
      data-deg-target={degTarget ?? ""}
      data-deg-reference={degReference ?? ""}
      className={`absolute top-0 left-16 z-[var(--z-panel)] w-[420px] ${glassPanel()}`}
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
                classNames={{
                  base: "font-semibold",
                }}
                color="primary"
                data-testid="deg-target"
                inputProps={{
                  classNames: {
                    input: "text-foreground font-semibold text-sm",
                    label: "text-primary-400",
                  },
                }}
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
                color="secondary"
                data-testid="deg-reference"
                inputProps={{
                  classNames: {
                    input: "text-foreground font-semibold text-sm",
                    label: "text-secondary-400",
                  },
                }}
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

            {/* Data-mode indicator + override */}
            <div className="flex items-center gap-2 px-2 text-[11px] text-default-400">
              <span
                title={`We detected this dataset as ${detectedMode === "log" ? "log-normalized" : "raw counts"} from its expression values (log-normalized data stays small, ~<15; raw counts run much higher). Fold change = mean of group 1 ÷ mean of group 2${mode === "log" ? " (means shown are un-logged so the ratio is a real fold change)" : ""}. Change it here if it's wrong.`}
              >
                Expression:
              </span>
              <select
                className="rounded bg-content2 px-1.5 py-0.5 text-[11px]"
                value={degDataMode}
                onChange={(e) =>
                  setDegDataMode(e.target.value as "auto" | "counts" | "log")
                }
              >
                <option value="auto">
                  auto ({detectedMode === "log" ? "log-norm" : "counts"})
                </option>
                <option value="counts">raw counts</option>
                <option value="log">log-normalized</option>
              </select>

              <button
                aria-label="How expression mode is detected"
                className={`flex h-4 w-4 items-center justify-center rounded-full border text-[9px] font-semibold transition-colors ${
                  showExprInfo
                    ? "border-primary text-primary"
                    : "border-default-400 text-default-400 hover:border-default-200 hover:text-default-200"
                }`}
                onClick={() => setShowExprInfo((v) => !v)}
              >
                i
              </button>

              {hasSmOverlay && (
                <div className="ml-auto flex items-center gap-1">
                  <span title="Choose what clicking a gene does: colour the cells by its single-cell expression, or add/remove it as a molecule layer in the overlay.">
                    Click adds to:
                  </span>
                  <div className="flex rounded bg-content2 p-0.5">
                    <button
                      className={`rounded px-2 py-0.5 transition-colors ${
                        clickTarget === "cell"
                          ? "bg-primary text-white"
                          : "text-default-400 hover:text-default-200"
                      }`}
                      onClick={() => setClickTarget("cell")}
                    >
                      Cell
                    </button>
                    <button
                      className={`rounded px-2 py-0.5 transition-colors ${
                        clickTarget === "molecule"
                          ? "bg-primary text-white"
                          : "text-default-400 hover:text-default-200"
                      }`}
                      onClick={() => setClickTarget("molecule")}
                    >
                      Molecule
                    </button>
                  </div>
                </div>
              )}
            </div>

            {showExprInfo && (
              <div className="mx-2 rounded-lg bg-content2/60 p-2.5 text-[11px] leading-relaxed text-default-400">
                <span className="text-default-200">
                  How the mode is detected.
                </span>{" "}
                We scan the mean expression of every gene in the current
                comparison and take the largest, currently{" "}
                <span className="font-mono text-default-200">
                  {formatMean(
                    deStats
                      ? deStats.means.reduce((a, b) => (b > a ? b : a), 0)
                      : 0,
                  )}
                </span>
                . Log-normalized data (log1p of size-normalized counts) almost
                never exceeds ~15, while raw counts of well-expressed genes run
                much higher — so a max above{" "}
                <span className="font-mono text-default-200">30</span> is read
                as raw counts, otherwise log-normalized. The mode only affects
                mean1 / mean2 and the fold change: for log data the means are
                un-logged (
                <span className="font-mono">expm1</span>) first so the ratio is
                a true linear fold change. Override it with the dropdown if the
                guess is wrong.
              </div>
            )}

            {/* Header */}
            <div className="grid grid-cols-[1.3fr_0.9fr_0.9fr_0.9fr] gap-2 text-[11px] tracking-wide text-default-300 px-2">
              <button
                className="text-left hover:text-default-200"
                onClick={() => handleSort("mean1")}
              >
                Gene
              </button>
              <SortHeader
                active={degSortKey === "mean1"}
                desc={degSortDesc}
                label="mean1"
                title={`Mean expression in ${degTarget ?? "group 1"}${mode === "log" ? " (un-logged)" : ""}`}
                onClick={() => handleSort("mean1")}
              />
              <span
                className="text-left text-default-400"
                title={`Mean expression in ${degReference ?? "the rest"}${mode === "log" ? " (un-logged)" : ""}`}
              >
                mean2
              </span>
              <SortHeader
                active={degSortKey === "fc"}
                desc={degSortDesc}
                label="FC"
                title="Fold change = mean1 ÷ mean2. >1 = higher in group 1."
                onClick={() => handleSort("fc")}
              />
            </div>

            {/* Rows */}
            <div className="max-h-[420px] overflow-y-auto flex flex-col">
              {filtered.length === 0 ? (
                <div className="text-xs text-default-400 py-4 text-center">
                  No genes match.
                </div>
              ) : (
                filtered.slice(0, 500).map((r, rowIndex) => {
                  const smViz =
                    clickTarget === "molecule"
                      ? smSelectedGenes.get(r.gene)
                      : undefined;
                  const isSelected =
                    clickTarget === "molecule"
                      ? !!smViz
                      : r.gene === selectedGene;
                  // In molecule mode, genes absent from the overlay dataset
                  // can't be shown — dim + disable them.
                  const notInSm =
                    clickTarget === "molecule" &&
                    !!smGeneSet &&
                    !smGeneSet.has(r.gene);
                  return (
                    <button
                      key={r.gene}
                      data-testid="deg-row"
                      data-gene={r.gene}
                      data-row-index={rowIndex}
                      disabled={notInSm}
                      className={`grid grid-cols-[1.3fr_0.9fr_0.9fr_0.9fr] gap-2 px-2 py-1.5 rounded text-left text-xs transition-colors ${
                        notInSm
                          ? "opacity-40 cursor-not-allowed"
                          : "hover:bg-white/10"
                      } ${isSelected ? "bg-primary/30" : ""}`}
                      onClick={() => onGeneClick(r.gene)}
                      title={
                        notInSm
                          ? "Not present in the molecule overlay dataset"
                          : `${r.pctIn > 0 || r.pctOut > 0 ? `% of cells expressing — group 1: ${(r.pctIn * 100).toFixed(1)}%, group 2: ${(r.pctOut * 100).toFixed(1)}%` : ""}`
                      }
                    >
                      <span className="font-mono truncate flex items-center gap-1.5">
                        {smViz && (
                          <span
                            className="inline-block h-2 w-2 shrink-0 rounded-full"
                            style={{ backgroundColor: smViz.color }}
                          />
                        )}
                        {r.gene}
                      </span>
                      <span className="font-mono text-foreground">
                        {formatMean(r.mean1)}
                      </span>
                      <span className="font-mono text-foreground/80">
                        {formatMean(r.mean2)}
                      </span>
                      <span
                        className={`font-mono ${
                          r.fc > 1
                            ? "text-rose-300"
                            : r.fc < 1
                              ? "text-sky-300"
                              : "text-default-300"
                        }`}
                      >
                        {formatFC(r.fc)}
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
  title,
  onClick,
}: {
  active: boolean;
  desc: boolean;
  label: string;
  title?: string;
  onClick: () => void;
}) {
  return (
    <button
      className={`text-left ${active ? "text-default-200" : "hover:text-default-200"}`}
      title={title}
      onClick={onClick}
    >
      {label}
      {active ? (desc ? " ↓" : " ↑") : ""}
    </button>
  );
}

function formatFC(v: number): string {
  if (!isFinite(v)) return "∞×";
  if (v === 0) return "0×";
  // Two sig-figs of scale: 12× / 3.2× / 0.45×.
  const digits = v >= 10 ? 0 : v >= 1 ? 1 : 2;

  return `${v.toFixed(digits)}×`;
}

function formatMean(v: number): string {
  if (!isFinite(v)) return "∞";
  if (v === 0) return "0";
  if (v >= 100) return v.toFixed(0);
  if (v >= 1) return v.toFixed(2);

  return v.toFixed(3);
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
