"use client";

import { Button } from "@heroui/button";
import { Input } from "@heroui/input";
import { Checkbox } from "@heroui/checkbox";
import { Switch } from "@heroui/switch";
import { Tooltip } from "@heroui/tooltip";
import { useMemo, useState } from "react";
import { toast } from "react-toastify";

import {
  usePanelSingleMoleculeStore,
  usePanelSingleMoleculeVisualizationStore,
} from "@/lib/hooks/usePanelStores";

function formatCount(n: number): string {
  if (n >= 1_000_000) return `${(n / 1_000_000).toFixed(1)}m`;
  if (n >= 1_000) return `${(n / 1_000).toFixed(1)}k`;

  return n.toString();
}

// Inline single-molecule gene picker: search, sort, Assigned/Unassigned toggles,
// scrollable checkbox list. Renders content only — caller provides any outer
// chrome (card, positioning).
export function SingleMoleculeGenePicker() {
  const [searchTerm, setSearchTerm] = useState("");
  const [sortBy, setSortBy] = useState<"alpha" | "count">("alpha");
  const [sortDir, setSortDir] = useState<"asc" | "desc">("asc");

  const dataset = usePanelSingleMoleculeStore((state) => {
    const id = state.currentDatasetId;

    return id ? state.datasets.get(id) : null;
  });

  const {
    selectedGenes,
    selectedGenesLegend,
    addGene,
    removeGene,
    clearGenes,
    showAssigned,
    setShowAssigned,
    showUnassigned,
    setShowUnassigned,
  } = usePanelSingleMoleculeVisualizationStore();

  const filteredGenes = useMemo(() => {
    if (!dataset) return [];

    const filtered = dataset.uniqueGenes.filter((gene) =>
      gene.toLowerCase().includes(searchTerm.toLowerCase()),
    );
    const counts = dataset.moleculeCounts;

    if (sortBy === "count" && counts) {
      const getTotal = (gene: string) => {
        const c = counts[gene];

        return c ? c.assigned + (c.unassigned ?? 0) : 0;
      };

      filtered.sort((a, b) =>
        sortDir === "asc"
          ? getTotal(a) - getTotal(b)
          : getTotal(b) - getTotal(a),
      );
    } else {
      filtered.sort((a, b) =>
        sortDir === "asc" ? a.localeCompare(b) : b.localeCompare(a),
      );
    }

    return filtered;
  }, [dataset, searchTerm, sortBy, sortDir]);

  if (!dataset) {
    return (
      <p className="text-sm text-default-400 text-center py-4">
        No single molecule dataset loaded
      </p>
    );
  }

  return (
    <div className="space-y-3">
      {/* Title */}
      <div className="flex items-center justify-between">
        <h3 className="text-lg font-semibold">Select Genes</h3>
        <div className="flex items-center gap-1">
          <span className="text-sm text-default-500">
            {selectedGenesLegend.size} selected
          </span>
          <Tooltip content="Copy all gene names" delay={300} placement="top">
            <button
              className="p-1 rounded hover:bg-default-200 transition-colors"
              onClick={() => {
                navigator.clipboard.writeText(dataset.uniqueGenes.join(","));
                toast.success("Gene names copied to clipboard");
              }}
            >
              <svg
                className="w-3.5 h-3.5 text-default-500"
                fill="none"
                stroke="currentColor"
                strokeWidth={1.5}
                viewBox="0 0 24 24"
              >
                <path
                  d="M15.75 17.25v3.375c0 .621-.504 1.125-1.125 1.125h-9.75a1.125 1.125 0 01-1.125-1.125V7.875c0-.621.504-1.125 1.125-1.125H6.75a9.06 9.06 0 011.5.124m7.5 10.376h3.375c.621 0 1.125-.504 1.125-1.125V11.25c0-4.46-3.243-8.161-7.5-8.876a9.06 9.06 0 00-1.5-.124H9.375c-.621 0-1.125.504-1.125 1.125v3.5m7.5 10.375H9.375a1.125 1.125 0 01-1.125-1.125v-9.25m12 6.625v-1.875a3.375 3.375 0 00-3.375-3.375h-1.5a1.125 1.125 0 01-1.125-1.125v-1.5a3.375 3.375 0 00-3.375-3.375H9.75"
                  strokeLinecap="round"
                  strokeLinejoin="round"
                />
              </svg>
            </button>
          </Tooltip>
        </div>
      </div>

      {/* Global Assigned/Unassigned Toggles */}
      {dataset.hasUnassigned && (
        <div className="flex flex-col gap-1 py-1">
          <div className="flex items-center justify-between">
            <span className="text-xs text-default-400">Show Assigned</span>
            <Switch
              isSelected={showAssigned}
              size="sm"
              onValueChange={setShowAssigned}
            />
          </div>
          <div className="flex items-center justify-between">
            <span className="text-xs text-default-400">Show Unassigned</span>
            <Switch
              isSelected={showUnassigned}
              size="sm"
              onValueChange={setShowUnassigned}
            />
          </div>
        </div>
      )}

      {/* Search */}
      <Input
        classNames={{ input: "text-sm" }}
        placeholder="Search genes..."
        size="sm"
        value={searchTerm}
        onValueChange={setSearchTerm}
      />

      {/* Clear All */}
      <Button
        className="w-full"
        color="danger"
        size="sm"
        variant="ghost"
        onPress={() => {
          clearGenes();
          setSearchTerm("");
        }}
      >
        Clear All
      </Button>

      {/* Sort buttons */}
      <div className="flex gap-1">
        <Button
          className="flex-1 text-xs"
          color={sortBy === "alpha" ? "primary" : "default"}
          size="sm"
          variant={sortBy === "alpha" ? "flat" : "light"}
          onPress={() => {
            if (sortBy === "alpha") {
              setSortDir(sortDir === "asc" ? "desc" : "asc");
            } else {
              setSortBy("alpha");
              setSortDir("asc");
            }
          }}
        >
          A→Z {sortBy === "alpha" ? (sortDir === "asc" ? "↑" : "↓") : ""}
        </Button>
        {dataset.moleculeCounts && (
          <Button
            className="flex-1 text-xs"
            color={sortBy === "count" ? "primary" : "default"}
            size="sm"
            variant={sortBy === "count" ? "flat" : "light"}
            onPress={() => {
              if (sortBy === "count") {
                setSortDir(sortDir === "asc" ? "desc" : "asc");
              } else {
                setSortBy("count");
                setSortDir("desc");
              }
            }}
          >
            # Count {sortBy === "count" ? (sortDir === "asc" ? "↑" : "↓") : ""}
          </Button>
        )}
      </div>

      {/* Gene List */}
      <div className="max-h-[400px] overflow-y-auto flex flex-col gap-0">
        {filteredGenes.length > 0 ? (
          filteredGenes.map((gene) => {
            const counts = dataset.moleculeCounts?.[gene];

            return (
              <Tooltip
                key={gene}
                content={
                  counts
                    ? `${counts.assigned.toLocaleString()} assigned${counts.unassigned != null ? ` | ${counts.unassigned.toLocaleString()} unassigned` : ""}`
                    : gene
                }
                delay={300}
                placement="right"
              >
                <div>
                  <Checkbox
                    className="w-full"
                    isSelected={selectedGenes.has(gene)}
                    size="sm"
                    onValueChange={() => {
                      if (selectedGenes.has(gene)) {
                        removeGene(gene);
                      } else {
                        addGene(gene);
                      }
                    }}
                  >
                    <span className="text-sm">
                      {gene}
                      {counts && (
                        <span className="text-default-400 ml-1">
                          {formatCount(counts.assigned)}
                          {counts.unassigned != null &&
                            ` | ${formatCount(counts.unassigned)}`}
                        </span>
                      )}
                    </span>
                  </Checkbox>
                </div>
              </Tooltip>
            );
          })
        ) : (
          <p className="text-sm text-default-400 text-center py-4">
            {searchTerm ? "No genes found" : "No genes available"}
          </p>
        )}
      </div>

      {/* Gene count info */}
      <div className="text-xs text-default-400 text-center pt-2 border-t border-default-200">
        Showing {filteredGenes.length} of {dataset.uniqueGenes.length || 0}{" "}
        genes
      </div>
    </div>
  );
}
