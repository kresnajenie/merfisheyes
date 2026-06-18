"use client";

import { Button } from "@heroui/button";
import { Input } from "@heroui/input";
import { Checkbox } from "@heroui/checkbox";
import { Switch } from "@heroui/switch";
import { Slider } from "@heroui/react";
import { Tooltip } from "@heroui/tooltip";
import { useState, useMemo, useRef, useEffect } from "react";
import { toast } from "react-toastify";

import { CameraPanel } from "./camera-panel";

import {
  usePanelSingleMoleculeStore,
  usePanelSingleMoleculeVisualizationStore,
  usePanelId,
} from "@/lib/hooks/usePanelStores";
import { useSplitScreenStore } from "@/lib/stores/splitScreenStore";
import { glassButton } from "@/components/primitives";
import { VISUALIZATION_CONFIG } from "@/lib/config/visualization.config";

function formatCount(n: number): string {
  if (n >= 1_000_000) return `${(n / 1_000_000).toFixed(1)}m`;
  if (n >= 1_000) return `${(n / 1_000).toFixed(1)}k`;

  return n.toString();
}

export function SingleMoleculeControls() {
  const { isSplitMode, enableSplit } = useSplitScreenStore();
  const panelId = usePanelId();
  const [isPanelOpen, setIsPanelOpen] = useState(false);
  const [searchTerm, setSearchTerm] = useState("");
  const [sortBy, setSortBy] = useState<"alpha" | "count">("alpha");
  const [sortDir, setSortDir] = useState<"asc" | "desc">("asc");

  // Get dataset
  const dataset = usePanelSingleMoleculeStore((state) => {
    const id = state.currentDatasetId;

    return id ? state.datasets.get(id) : null;
  });

  // Get visualization state
  const {
    selectedGenes,
    selectedGenesLegend,
    addGene,
    removeGene,
    clearGenes,
    globalScale,
    setGlobalScale,
    viewMode,
    setViewMode,
    showAssigned,
    setShowAssigned,
    showUnassigned,
    setShowUnassigned,
    sceneRotation,
    setSceneRotation,
    flipX,
    setFlipX,
    flipY,
    setFlipY,
  } = usePanelSingleMoleculeVisualizationStore();

  const [isCameraOpen, setIsCameraOpen] = useState(false);
  const controlsRef = useRef<HTMLDivElement>(null);

  // Track shift key for 45° snap on rotation slider
  useEffect(() => {
    const handleKeyDown = (e: KeyboardEvent) => {
      if (e.key === "Shift") (window as any).__shiftHeld = true;
    };
    const handleKeyUp = (e: KeyboardEvent) => {
      if (e.key === "Shift") (window as any).__shiftHeld = false;
    };
    window.addEventListener("keydown", handleKeyDown);
    window.addEventListener("keyup", handleKeyUp);
    return () => {
      window.removeEventListener("keydown", handleKeyDown);
      window.removeEventListener("keyup", handleKeyUp);
    };
  }, []);

  // Filter and sort genes
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
        sortDir === "asc"
          ? a.localeCompare(b)
          : b.localeCompare(a),
      );
    }

    return filtered;
  }, [dataset, searchTerm, sortBy, sortDir]);

  return (
    <div ref={controlsRef} className="absolute top-28 left-4 z-50 flex flex-col gap-2">
      {/* Gene Selection Button */}
      <Button
        className={`w-14 h-14 min-w-0 rounded-full font-medium text-xs ${
          isPanelOpen ? "" : glassButton()
        }`}
        color={isPanelOpen ? "primary" : "default"}
        variant={isPanelOpen ? "shadow" : "light"}
        onPress={() => setIsPanelOpen(!isPanelOpen)}
      >
        Genes
      </Button>

      {/* Dot Size Slider */}
      <Tooltip content="Change dotsize" placement="right">
        <div
          className={`w-14 h-32 rounded-full border-2 border-default-200 p-2 flex flex-col items-center justify-center ${glassButton()}`}
        >
          <Slider
            aria-label="Dot size"
            className="h-full"
            maxValue={VISUALIZATION_CONFIG.SINGLE_MOLECULE_GLOBAL_SCALE_MAX}
            minValue={VISUALIZATION_CONFIG.SINGLE_MOLECULE_GLOBAL_SCALE_MIN}
            orientation="vertical"
            size="sm"
            step={VISUALIZATION_CONFIG.SINGLE_MOLECULE_GLOBAL_SCALE_STEP}
            value={globalScale}
            onChange={(value) => setGlobalScale(value as number)}
          />
        </div>
      </Tooltip>

      {/* Camera Button */}
      <Tooltip content="Camera controls" placement="right">
        <Button
          className={`w-14 h-14 min-w-0 rounded-full font-medium text-xs ${isCameraOpen ? "" : glassButton()}`}
          color={isCameraOpen ? "primary" : "default"}
          variant={isCameraOpen ? "shadow" : "light"}
          onPress={() => {
            setIsCameraOpen(!isCameraOpen);
            if (!isCameraOpen) setIsPanelOpen(false);
          }}
        >
          <svg className="w-5 h-5" fill="none" stroke="currentColor" strokeWidth={1.5} viewBox="0 0 24 24">
            <path d="M6.827 6.175A2.31 2.31 0 015.186 7.23c-.38.054-.757.112-1.134.175C2.999 7.58 2.25 8.507 2.25 9.574V18a2.25 2.25 0 002.25 2.25h15A2.25 2.25 0 0021.75 18V9.574c0-1.067-.75-1.994-1.802-2.169a47.865 47.865 0 00-1.134-.175 2.31 2.31 0 01-1.64-1.055l-.822-1.316a2.192 2.192 0 00-1.736-1.039 48.774 48.774 0 00-5.232 0 2.192 2.192 0 00-1.736 1.039l-.821 1.316z" strokeLinecap="round" strokeLinejoin="round" />
            <path d="M16.5 12.75a4.5 4.5 0 11-9 0 4.5 4.5 0 019 0z" strokeLinecap="round" strokeLinejoin="round" />
          </svg>
        </Button>
      </Tooltip>

      {/* Split Screen Button */}
      {!isSplitMode && !panelId && (
        <Tooltip content="Split screen" placement="right">
          <Button
            className={`w-14 h-14 min-w-0 rounded-full font-medium text-xs ${glassButton()}`}
            color="default"
            variant="light"
            onPress={enableSplit}
          >
            <svg
              className="w-5 h-5"
              fill="none"
              stroke="currentColor"
              strokeWidth={1.5}
              viewBox="0 0 24 24"
            >
              <path
                d="M9 4.5v15m6-15v15M4.5 19.5h15a1.5 1.5 0 001.5-1.5V6a1.5 1.5 0 00-1.5-1.5h-15A1.5 1.5 0 003 6v12a1.5 1.5 0 001.5 1.5z"
                strokeLinecap="round"
                strokeLinejoin="round"
              />
            </svg>
          </Button>
        </Tooltip>
      )}

      {/* Gene Selection Panel */}
      {isPanelOpen && (
        <div
          className={`absolute top-0 left-16 z-50 w-[320px] border-2 border-white/20 rounded-3xl shadow-lg ${glassButton()}`}
        >
          <div className="p-4 space-y-3">
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
                      if (dataset) {
                        navigator.clipboard.writeText(dataset.uniqueGenes.join(","));
                        toast.success("Gene names copied to clipboard");
                      }
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
            {dataset?.hasUnassigned && (
              <div className="flex flex-col gap-1 py-1">
                <div className="flex items-center justify-between">
                  <span className="text-xs text-default-400">
                    Show Assigned
                  </span>
                  <Switch
                    isSelected={showAssigned}
                    size="sm"
                    onValueChange={setShowAssigned}
                  />
                </div>
                <div className="flex items-center justify-between">
                  <span className="text-xs text-default-400">
                    Show Unassigned
                  </span>
                  <Switch
                    isSelected={showUnassigned}
                    size="sm"
                    onValueChange={setShowUnassigned}
                  />
                </div>
              </div>
            )}

            {/* Search Input */}
            <Input
              classNames={{
                input: "text-sm",
              }}
              placeholder="Search genes..."
              size="sm"
              value={searchTerm}
              onValueChange={setSearchTerm}
            />

            {/* Clear Button */}
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

            {/* Sort Buttons */}
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
              {dataset?.moleculeCounts && (
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
                  const counts = dataset?.moleculeCounts?.[gene];

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
              Showing {filteredGenes.length} of{" "}
              {dataset?.uniqueGenes.length || 0} genes
            </div>
          </div>
        </div>
      )}

      {/* Camera Panel */}
      {isCameraOpen && (
        <CameraPanel
          controlsRef={controlsRef}
          onClose={() => setIsCameraOpen(false)}
          sceneRotation={sceneRotation}
          setSceneRotation={setSceneRotation}
          flipX={flipX}
          setFlipX={setFlipX}
          flipY={flipY}
          setFlipY={setFlipY}
          viewMode={viewMode}
          setViewMode={setViewMode}
        />
      )}
    </div>
  );
}
