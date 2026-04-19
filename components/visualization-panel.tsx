"use client";

import type { VisualizationMode } from "@/lib/stores/createVisualizationStore";
import type { StandardizedDataset } from "@/lib/StandardizedDataset";

import { Input } from "@heroui/input";
import { Button } from "@heroui/button";
import { Checkbox } from "@heroui/checkbox";
import { RadioGroup, Radio } from "@heroui/radio";
import { useState, useMemo, useRef, useEffect, useCallback } from "react";

import { Autocomplete, AutocompleteItem } from "@heroui/autocomplete";
import { Slider, Textarea } from "@heroui/react";
import { Modal, ModalContent, ModalHeader, ModalBody, ModalFooter } from "@heroui/modal";
import { Tooltip } from "@heroui/tooltip";
import { toast } from "react-toastify";

import {
  usePanelDatasetStore,
  usePanelVisualizationStore,
} from "@/lib/hooks/usePanelStores";
import { glassButton } from "@/components/primitives";
import {
  getEffectiveColumnType,
  canTreatAsNumerical,
} from "@/lib/utils/column-type-utils";

interface VisualizationPanelProps {
  mode: VisualizationMode; // This is panelMode, not the visualization mode array
  onClose: () => void;
  controlsRef?: React.RefObject<HTMLDivElement>;
}

export function VisualizationPanel({
  mode,
  onClose,
  controlsRef,
}: VisualizationPanelProps) {
  const [searchTerm, setSearchTerm] = useState("");
  const [currentPage, setCurrentPage] = useState(1);
  const itemsPerPage = 1000; // Show 1000 items per page
  const panelRef = useRef<HTMLDivElement>(null);
  const [isSequenceModalOpen, setIsSequenceModalOpen] = useState(false);
  const [sequenceText, setSequenceText] = useState("");
  const { getCurrentDataset } = usePanelDatasetStore();
  const {
    selectedColumn,
    setSelectedColumn,
    selectedCelltypes,
    setCelltypes,
    selectedGene,
    toggleCelltype,
    celltypePlayback: isPlaying,
    celltypePlaybackInterval: playInterval,
    celltypePlaybackSequence,
    setCelltypePlayback,
    setCelltypePlaybackInterval,
    setCelltypePlaybackSequence,
    setSelectedGene,
    setMode,
    celltypeSearchTerm,
    geneSearchTerm,
    setCelltypeSearchTerm,
    setGeneSearchTerm,
    clusterVersion,
    incrementClusterVersion,
    columnTypeOverrides,
    toggleColumnType,
    pinnedTooltipColumns,
    togglePinnedTooltipColumn,
  } = usePanelVisualizationStore();

  const currentSearchTerm =
    mode === "celltype" ? celltypeSearchTerm : geneSearchTerm;
  const updateSearchTerm =
    mode === "celltype" ? setCelltypeSearchTerm : setGeneSearchTerm;

  // Handle click outside to close panel
  useEffect(() => {
    const handleClickOutside = (event: MouseEvent) => {
      // Don't close if a modal is open
      if (isSequenceModalOpen) return;

      const target = event.target as Node;

      // Don't close if clicking inside the panel
      if (panelRef.current && panelRef.current.contains(target)) {
        return;
      }

      // Don't close if clicking on the controls
      if (controlsRef?.current && controlsRef.current.contains(target)) {
        return;
      }

      // Don't close if clicking inside a popover/listbox/modal (rendered via portal)
      if (
        target instanceof HTMLElement &&
        target.closest('[role="listbox"], [data-slot="content"], [role="dialog"], .nextui-modal-backdrop, [data-slot="backdrop"]')
      ) {
        return;
      }

      // Close if clicking outside both panel and controls
      onClose();
    };

    // Add event listener
    document.addEventListener("mousedown", handleClickOutside);

    // Cleanup
    return () => {
      document.removeEventListener("mousedown", handleClickOutside);
    };
  }, [onClose, controlsRef]);

  // Type guard: this component only works with StandardizedDataset
  const rawDataset = getCurrentDataset();
  const dataset =
    rawDataset && "clusters" in rawDataset
      ? (rawDataset as StandardizedDataset)
      : null;

  // Get all cluster columns for dropdown (including unloaded ones)
  const clusterColumns = useMemo(() => {
    if (!dataset) return [];

    const loadedColumns = new Set(
      (dataset.clusters || []).map((c) => c.column),
    );

    // Use allClusterColumnNames if available (S3/chunked datasets)
    if (dataset.allClusterColumnNames?.length > 0) {
      return dataset.allClusterColumnNames.map((name) => ({
        key: name,
        label: name,
        type: getEffectiveColumnType(name, dataset, columnTypeOverrides),
        loaded: loadedColumns.has(name),
      }));
    }

    // Fallback: use loaded clusters only (local file uploads)
    return (dataset.clusters || []).map((cluster) => ({
      key: cluster.column,
      label: cluster.column,
      type: getEffectiveColumnType(
        cluster.column,
        dataset,
        columnTypeOverrides,
      ),
      loaded: true,
    }));
  }, [dataset, clusterVersion, columnTypeOverrides]);

  // Check if the selected column is numerical
  const isNumericalColumn = useMemo(() => {
    if (!dataset || !selectedColumn) return false;

    return (
      getEffectiveColumnType(selectedColumn, dataset, columnTypeOverrides) ===
      "numerical"
    );
  }, [dataset, selectedColumn, clusterVersion, columnTypeOverrides]);

  // Whether the selected column can be toggled to numerical
  // (only show toggle when the values are actually number-like)
  const showTypeToggle = useMemo(() => {
    if (!dataset || !selectedColumn) return false;

    // Always allow switching back from numerical to categorical
    if (isNumericalColumn) return true;

    // Only offer "Treat as Numerical" when values parse as numbers
    return canTreatAsNumerical(selectedColumn, dataset);
  }, [dataset, selectedColumn, isNumericalColumn, clusterVersion]);

  // Get items based on mode and dataset
  const items = useMemo(() => {
    if (!dataset) return [];

    switch (mode) {
      case "celltype": {
        // Get unique cell types from the selected cluster column
        if (!dataset.clusters || !selectedColumn) return [];

        // Find the cluster data for the selected column
        const selectedCluster = dataset.clusters.find(
          (c) => c.column === selectedColumn,
        );

        if (!selectedCluster) return [];

        // Skip for numerical columns (list is hidden — only column dropdown shown)
        if (
          getEffectiveColumnType(
            selectedColumn,
            dataset,
            columnTypeOverrides,
          ) === "numerical"
        )
          return [];

        // Get unique values then always sort alphabetically
        const palette = selectedCluster.palette || {};
        const uniqueVals = (selectedCluster.uniqueValues
          ?? [...new Set(selectedCluster.values.map(String))]
        ).slice().sort((a, b) => a.localeCompare(b, undefined, { numeric: true }));

        return uniqueVals.map((val) => ({
          id: val,
          label: val,
          color: palette[val] || "#808080",
        }));
      }

      case "gene": {
        // Get all genes in deterministic, human-friendly order
        const sortedGenes = [...dataset.genes].sort((a, b) =>
          a.localeCompare(b, undefined, {
            numeric: true,
            sensitivity: "base",
          }),
        );

        return sortedGenes.map((gene) => ({
          id: gene,
          label: gene,
          color: "#FFFFFF",
        }));
      }

      default:
        return [];
    }
  }, [dataset, mode, selectedColumn, clusterVersion]);

  const filteredItems = items.filter((item) =>
    item.label.toLowerCase().includes(currentSearchTerm.toLowerCase()),
  );

  // Reset to page 1 when search term changes
  useEffect(() => {
    setCurrentPage(1);
  }, [currentSearchTerm, mode]);

  // Pagination logic
  const totalPages = Math.ceil(filteredItems.length / itemsPerPage);
  const startIndex = (currentPage - 1) * itemsPerPage;
  const endIndex = startIndex + itemsPerPage;
  const paginatedItems = filteredItems.slice(startIndex, endIndex);

  const stopPlayback = useCallback(() => {
    setCelltypePlayback(false);
  }, [setCelltypePlayback]);

  const getTitle = () => {
    switch (mode) {
      case "celltype":
        return "celltype";
      case "gene":
        return "gene";
      default:
        return mode;
    }
  };

  return (
    <>
    <div
      ref={panelRef}
      className={`absolute top-0 left-16 z-50 w-[300px] border-2 border-white/20 rounded-3xl shadow-lg ${glassButton()}`}
    >
      {/* Content */}
      <div className="p-4 space-y-3">
        {/* Select celltypes - only show for celltype mode */}
        {mode === "celltype" && (
          <>
          <Autocomplete
            className="max-w-xs"
            color="primary"
            label="Cluster Column"
            placeholder="Select a column"
            selectedKey={selectedColumn}
            onSelectionChange={async (key) => {
              const columnKey = (key as string) || null;

              if (!columnKey || !dataset) {
                setSelectedColumn(null);
                return;
              }

              // Check if this column is already loaded
              const loadedCluster = dataset.clusters?.find(
                (c) => c.column === columnKey,
              );

              if (loadedCluster) {
                // Already loaded — use effective type (respects overrides)
                const isNumerical =
                  getEffectiveColumnType(
                    columnKey,
                    dataset,
                    columnTypeOverrides,
                  ) === "numerical";

                setSelectedColumn(columnKey, isNumerical);
                setMode(["celltype"]);
              } else if (dataset.adapter) {
                // Not loaded yet — fetch on demand
                const isNumerical =
                  getEffectiveColumnType(
                    columnKey,
                    dataset,
                    columnTypeOverrides,
                  ) === "numerical";

                setSelectedColumn(columnKey, isNumerical);
                setMode(["celltype"]);

                const toastId = toast.loading(
                  `Loading cluster column "${columnKey}"...`,
                );

                try {
                  let newClusters: Array<{
                    column: string;
                    type: string;
                    values: any[];
                    palette: Record<string, string> | null;
                    uniqueValues?: string[];
                  }> | null = null;

                  if (dataset.adapter.mode === "local") {
                    // Local datasets: load on main thread
                    newClusters = await dataset.adapter.loadClusters([
                      columnKey,
                    ]);
                  } else {
                    // S3 / custom S3: load in worker to avoid freezing UI
                    const { getStandardizedDatasetWorker } = await import(
                      "@/lib/workers/standardizedDatasetWorkerManager"
                    );
                    const worker = await getStandardizedDatasetWorker();

                    newClusters = await worker.loadClusterFromS3(
                      dataset.id,
                      [columnKey],
                      dataset.metadata?.customS3BaseUrl,
                    );
                  }

                  if (newClusters && newClusters.length > 0) {
                    dataset.addClusters(newClusters);
                    incrementClusterVersion();
                  }
                  toast.dismiss(toastId);
                } catch (error) {
                  console.warn(
                    `Failed to load cluster column ${columnKey}:`,
                    error,
                  );
                  toast.update(toastId, {
                    render: `Failed to load "${columnKey}"`,
                    type: "error",
                    isLoading: false,
                    autoClose: 3000,
                  });
                }
              }
            }}
          >
            {clusterColumns.map((column) => (
              <AutocompleteItem
                key={column.key}
                startContent={
                  <Tooltip
                    content={
                      column.type === "numerical" ? "Numerical" : "Categorical"
                    }
                    delay={300}
                    placement="left"
                  >
                    {column.type === "numerical" ? (
                      <svg
                        className="w-4 h-4 text-default-500 flex-shrink-0"
                        fill="none"
                        stroke="currentColor"
                        strokeWidth={2}
                        viewBox="0 0 24 24"
                      >
                        <path
                          d="M3 3v18h18"
                          strokeLinecap="round"
                          strokeLinejoin="round"
                        />
                        <path
                          d="M7 16l4-8 4 5 5-9"
                          strokeLinecap="round"
                          strokeLinejoin="round"
                        />
                      </svg>
                    ) : (
                      <svg
                        className="w-4 h-4 text-default-500 flex-shrink-0"
                        fill="none"
                        stroke="currentColor"
                        strokeWidth={2}
                        viewBox="0 0 24 24"
                      >
                        <circle cx="8" cy="8" r="3" />
                        <circle cx="16" cy="8" r="3" />
                        <circle cx="8" cy="16" r="3" />
                        <circle cx="16" cy="16" r="3" />
                      </svg>
                    )}
                  </Tooltip>
                }
                endContent={
                  <div className="flex items-center gap-1.5">
                    {column.loaded && (
                      <span className="w-2 h-2 rounded-full bg-green-500 flex-shrink-0" />
                    )}
                    <Tooltip content={pinnedTooltipColumns.has(column.key) ? "Unpin from tooltip" : "Pin to tooltip"} delay={300} placement="right">
                      <button
                        className="p-0.5 rounded hover:bg-white/10 transition-colors"
                        onClick={(e) => {
                          e.stopPropagation();
                          e.preventDefault();
                          togglePinnedTooltipColumn(column.key);
                        }}
                        onMouseDown={(e) => e.stopPropagation()}
                      >
                        <svg
                          className={`w-3.5 h-3.5 flex-shrink-0 transition-colors ${pinnedTooltipColumns.has(column.key) ? "text-primary" : "text-default-400"}`}
                          fill={pinnedTooltipColumns.has(column.key) ? "currentColor" : "none"}
                          stroke="currentColor"
                          strokeWidth={2}
                          viewBox="0 0 24 24"
                        >
                          <path d="M12 2l-2 7H4l6 4.5L8 21l4-3 4 3-2-7.5L20 9h-6z" strokeLinecap="round" strokeLinejoin="round" />
                        </svg>
                      </button>
                    </Tooltip>
                  </div>
                }
              >
                {column.label}
              </AutocompleteItem>
            ))}
          </Autocomplete>

          {/* Celltype playback controls — only for categorical columns */}
          {selectedColumn && !isNumericalColumn && items.length > 0 && (
            <div className="flex flex-col gap-2">
              <div className="flex items-center gap-2">
                <Tooltip content={isPlaying ? "Stop" : "Play through celltypes"} placement="bottom">
                  <Button
                    className="min-w-0 w-10 h-8"
                    color={isPlaying ? "danger" : "primary"}
                    size="sm"
                    variant={isPlaying ? "solid" : "flat"}
                    onPress={() => {
                      if (isPlaying) {
                        setCelltypePlayback(false);
                      } else {
                        setCelltypePlayback(true);
                        onClose();
                      }
                    }}
                  >
                    {isPlaying ? (
                      <svg className="w-4 h-4" fill="currentColor" viewBox="0 0 24 24">
                        <rect x="6" y="4" width="4" height="16" />
                        <rect x="14" y="4" width="4" height="16" />
                      </svg>
                    ) : (
                      <svg className="w-4 h-4" fill="currentColor" viewBox="0 0 24 24">
                        <polygon points="5,3 19,12 5,21" />
                      </svg>
                    )}
                  </Button>
                </Tooltip>
                <Tooltip content="Edit playback sequence" placement="bottom">
                  <Button
                    className="min-w-0 w-10 h-8"
                    size="sm"
                    variant="flat"
                    onPress={() => {
                      setSequenceText(celltypePlaybackSequence.join(", "));
                      setIsSequenceModalOpen(true);
                    }}
                  >
                    <svg className="w-4 h-4" fill="none" stroke="currentColor" strokeWidth={2} viewBox="0 0 24 24">
                      <path d="M4 6h16M4 12h16M4 18h7" strokeLinecap="round" />
                    </svg>
                  </Button>
                </Tooltip>
                <div className="flex-1 flex items-center gap-2">
                  <Slider
                    aria-label="Playback speed"
                    className="flex-1"
                    maxValue={5}
                    minValue={0.2}
                    size="sm"
                    step={0.1}
                    value={playInterval}
                    onChange={(v) => setCelltypePlaybackInterval(v as number)}
                  />
                  <Input
                    aria-label="Interval seconds"
                    className="w-14"
                    classNames={{ input: "text-xs text-center" }}
                    size="sm"
                    type="number"
                    value={playInterval.toFixed(1)}
                    onValueChange={(v) => {
                      const num = parseFloat(v);
                      if (!isNaN(num) && num >= 0.1) setCelltypePlaybackInterval(num);
                    }}
                  />
                  <span className="text-xs text-default-400">s</span>
                </div>
              </div>
              {celltypePlaybackSequence.length > 0 && (
                <div className="text-xs text-default-400">
                  Sequence: {celltypePlaybackSequence.length} celltypes
                  <Button
                    className="min-w-0 h-5 ml-1"
                    color="danger"
                    size="sm"
                    variant="light"
                    onPress={() => setCelltypePlaybackSequence([])}
                  >
                    Clear
                  </Button>
                </div>
              )}
            </div>
          )}

          {/* Column type toggle for the selected column */}
          {showTypeToggle && selectedColumn && (
            <Button
              className="w-full"
              size="sm"
              variant="flat"
              startContent={
                isNumericalColumn ? (
                  <svg
                    className="w-4 h-4 flex-shrink-0"
                    fill="none"
                    stroke="currentColor"
                    strokeWidth={2}
                    viewBox="0 0 24 24"
                  >
                    <path
                      d="M3 3v18h18"
                      strokeLinecap="round"
                      strokeLinejoin="round"
                    />
                    <path
                      d="M7 16l4-8 4 5 5-9"
                      strokeLinecap="round"
                      strokeLinejoin="round"
                    />
                  </svg>
                ) : (
                  <svg
                    className="w-4 h-4 flex-shrink-0"
                    fill="none"
                    stroke="currentColor"
                    strokeWidth={2}
                    viewBox="0 0 24 24"
                  >
                    <circle cx="8" cy="8" r="3" />
                    <circle cx="16" cy="8" r="3" />
                    <circle cx="8" cy="16" r="3" />
                    <circle cx="16" cy="16" r="3" />
                  </svg>
                )
              }
              onPress={() => {
                const currentType = isNumericalColumn
                  ? "numerical"
                  : "categorical";
                const newType =
                  currentType === "categorical" ? "numerical" : "categorical";

                toggleColumnType(selectedColumn, currentType as "categorical" | "numerical");
                toast.info(`"${selectedColumn}" set to ${newType}`, {
                  autoClose: 2000,
                });
              }}
            >
              Treat as {isNumericalColumn ? "Categorical" : "Numerical"}
            </Button>
          )}
          </>
        )}

        {/* Only show search, clear, and list if not numerical column in celltype mode */}
        {!(mode === "celltype" && isNumericalColumn) && (
          <>
            {/* Search Input */}
            <Input
              classNames={{
                input: "text-sm",
              }}
              placeholder={`Search ${getTitle()}`}
              value={currentSearchTerm}
              onValueChange={updateSearchTerm}
            />

            {/* Clear Button + Copy Genes */}
            <div className="flex gap-1">
              <Button
                className="flex-1"
                color="danger"
                variant="ghost"
                onPress={() => {
                  updateSearchTerm("");
                  if (mode === "celltype") {
                    // Clear all selected celltypes
                    selectedCelltypes.forEach((ct) => toggleCelltype(ct));
                  } else if (mode === "gene") {
                    // Clear selected gene
                    setSelectedGene(null);
                  }
                }}
              >
                Clear
              </Button>
              {mode === "gene" && dataset && (
                <Tooltip content="Copy all gene names" delay={300} placement="top">
                  <Button
                    className="min-w-0 px-3"
                    variant="ghost"
                    onPress={() => {
                      navigator.clipboard.writeText(dataset.genes.join(","));
                      toast.success("Gene names copied to clipboard");
                    }}
                  >
                    <svg
                      className="w-4 h-4 text-default-500"
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
                  </Button>
                </Tooltip>
              )}
            </div>

            {/* Pagination Info */}
            {filteredItems.length > itemsPerPage && (
              <div className="text-xs text-default-500 text-center">
                Showing {startIndex + 1}-
                {Math.min(endIndex, filteredItems.length)} of{" "}
                {filteredItems.length}
              </div>
            )}

            {/* List */}
            <div className="max-h-[400px] overflow-y-auto flex flex-col gap-0">
              {mode === "celltype" ? (
                // Checkbox list for celltype mode
                paginatedItems.map((item) => (
                  <Checkbox
                    key={item.id}
                    className="w-full"
                    isSelected={selectedCelltypes.has(item.id)}
                    size="sm"
                    onValueChange={() => {
                      if (isPlaying) stopPlayback();
                      toggleCelltype(item.id);
                    }}
                  >
                    <span style={{ color: item.color }}>{item.label}</span>
                  </Checkbox>
                ))
              ) : mode === "gene" ? (
                // Radio group for gene mode
                <RadioGroup
                  value={selectedGene || ""}
                  onValueChange={(value) => {
                    setSelectedGene(value || null);
                    // Mode is now automatically updated by setSelectedGene
                  }}
                >
                  {paginatedItems.map((item) => (
                    <Radio key={item.id} size="sm" value={item.id}>
                      <span style={{ color: item.color }}>{item.label}</span>
                    </Radio>
                  ))}
                </RadioGroup>
              ) : (
                // For other modes, just display items
                paginatedItems.map((item) => (
                  <div key={item.id} className="p-2">
                    <span style={{ color: item.color }}>{item.label}</span>
                  </div>
                ))
              )}
            </div>

            {/* Pagination Controls */}
            {totalPages > 1 && (
              <div className="flex justify-between items-center gap-2">
                <Button
                  isDisabled={currentPage === 1}
                  size="sm"
                  variant="flat"
                  onPress={() => setCurrentPage((p) => Math.max(1, p - 1))}
                >
                  Previous
                </Button>
                <span className="text-xs text-default-500">
                  Page {currentPage} of {totalPages}
                </span>
                <Button
                  isDisabled={currentPage === totalPages}
                  size="sm"
                  variant="flat"
                  onPress={() =>
                    setCurrentPage((p) => Math.min(totalPages, p + 1))
                  }
                >
                  Next
                </Button>
              </div>
            )}
          </>
        )}
      </div>
      {/* Slider Section
      <div className="p-4 border-t border-default-200">
        <div className="space-y-2">
          <input
            type="range"
            min="0"
            max="3.39"
            step="0.01"
            defaultValue="0.92"
            className="w-full"
          />
          <div className="flex justify-between text-xs text-default-500">
            <span>0.00</span>
            <span className="text-default-900 font-semibold">0.92</span>
            <span>3.39</span>
          </div>
        </div>
      </div> */}
    </div>

    {/* Playback sequence modal */}
    <Modal
      isOpen={isSequenceModalOpen}
      onClose={() => setIsSequenceModalOpen(false)}
      size="lg"
    >
      <ModalContent>
        <ModalHeader>Playback Sequence</ModalHeader>
        <ModalBody>
          <p className="text-sm text-default-500 mb-2">
            Enter celltype names separated by commas. Playback will iterate through them in order.
            Leave empty to play all celltypes.
          </p>
          <Textarea
            label="Celltype sequence"
            placeholder="e.g. Neuron, Astrocyte, Microglia"
            value={sequenceText}
            onValueChange={setSequenceText}
            minRows={3}
            maxRows={10}
          />
          {sequenceText.trim() && (
            <div className="text-xs text-default-400 mt-1">
              {sequenceText.split(",").map((s) => s.trim()).filter(Boolean).length} celltypes in sequence
            </div>
          )}
        </ModalBody>
        <ModalFooter>
          <Button
            color="danger"
            variant="light"
            onPress={() => {
              setSequenceText("");
              setCelltypePlaybackSequence([]);
              setIsSequenceModalOpen(false);
            }}
          >
            Clear
          </Button>
          <Button
            color="primary"
            onPress={() => {
              const sequence = sequenceText
                .split(",")
                .map((s) => s.trim())
                .filter(Boolean);

              // Validate against available items
              const validSet = new Set(items.map((i) => i.id));
              const invalid = sequence.filter((s) => !validSet.has(s));
              if (invalid.length > 0) {
                toast.warning(
                  `Skipping unknown celltypes: ${invalid.join(", ")}`,
                  { autoClose: 4000 },
                );
              }

              const valid = sequence.filter((s) => validSet.has(s));
              setCelltypePlaybackSequence(valid);
              setIsSequenceModalOpen(false);
            }}
          >
            Save
          </Button>
        </ModalFooter>
      </ModalContent>
    </Modal>
    </>
  );
}
