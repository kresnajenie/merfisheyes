"use client";

import type { VisualizationMode } from "@/lib/stores/createVisualizationStore";
import type { StandardizedDataset } from "@/lib/StandardizedDataset";

import { Input } from "@heroui/input";
import { Button } from "@heroui/button";
import { Checkbox } from "@heroui/checkbox";
import { RadioGroup, Radio } from "@heroui/radio";
import { useState, useMemo, useRef, useEffect } from "react";
import { Autocomplete, AutocompleteItem } from "@heroui/autocomplete";
import { Tooltip } from "@heroui/tooltip";
import { toast } from "react-toastify";

import {
  usePanelDatasetStore,
  usePanelVisualizationStore,
} from "@/lib/hooks/usePanelStores";
import { glassButton } from "@/components/primitives";

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
  const { getCurrentDataset } = usePanelDatasetStore();
  const {
    selectedColumn,
    setSelectedColumn,
    selectedCelltypes,
    selectedGene,
    toggleCelltype,
    setSelectedGene,
    setMode,
    celltypeSearchTerm,
    geneSearchTerm,
    setCelltypeSearchTerm,
    setGeneSearchTerm,
    clusterVersion,
    incrementClusterVersion,
  } = usePanelVisualizationStore();

  const currentSearchTerm =
    mode === "celltype" ? celltypeSearchTerm : geneSearchTerm;
  const updateSearchTerm =
    mode === "celltype" ? setCelltypeSearchTerm : setGeneSearchTerm;

  // Handle click outside to close panel
  useEffect(() => {
    const handleClickOutside = (event: MouseEvent) => {
      const target = event.target as Node;

      // Don't close if clicking inside the panel
      if (panelRef.current && panelRef.current.contains(target)) {
        return;
      }

      // Don't close if clicking on the controls
      if (controlsRef?.current && controlsRef.current.contains(target)) {
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
        type: (dataset.allClusterColumnTypes?.[name] || "categorical") as
          | "categorical"
          | "numerical",
        loaded: loadedColumns.has(name),
      }));
    }

    // Fallback: use loaded clusters only (local file uploads)
    return (dataset.clusters || []).map((cluster) => ({
      key: cluster.column,
      label: cluster.column,
      type: cluster.type as "categorical" | "numerical",
      loaded: true,
    }));
  }, [dataset, clusterVersion]);

  // Check if the selected column is numerical
  const isNumericalColumn = useMemo(() => {
    if (!dataset?.clusters || !selectedColumn) return false;

    const selectedCluster = dataset.clusters.find(
      (c) => c.column === selectedColumn,
    );

    return selectedCluster?.type === "numerical";
  }, [dataset, selectedColumn, clusterVersion]);

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

        // Skip expensive unique-value computation for numerical columns
        // (the list is hidden anyway — only the column dropdown is shown)
        if (selectedCluster.type === "numerical") return [];

        const uniqueTypes = new Set<string>();
        const itemsMap = new Map<
          string,
          { id: string; label: string; color: string }
        >();

        selectedCluster.values.forEach((value) => {
          const typeStr = String(value);

          if (!uniqueTypes.has(typeStr)) {
            uniqueTypes.add(typeStr);
            const palette = selectedCluster.palette || {};

            itemsMap.set(typeStr, {
              id: typeStr,
              label: typeStr,
              color: palette[typeStr] || "#808080",
            });
          }
        });

        return Array.from(itemsMap.values()).sort((a, b) =>
          a.label.localeCompare(b.label, undefined, { numeric: true }),
        );
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
    <div
      ref={panelRef}
      className={`absolute top-0 left-16 z-50 w-[300px] border-2 border-white/20 rounded-3xl shadow-lg ${glassButton()}`}
    >
      {/* Content */}
      <div className="p-4 space-y-3">
        {/* Select celltypes - only show for celltype mode */}
        {mode === "celltype" && (
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
                // Already loaded — use it directly
                const isNumerical = loadedCluster.type === "numerical";

                setSelectedColumn(columnKey, isNumerical);
                setMode(["celltype"]);
              } else if (dataset.adapter) {
                // Not loaded yet — fetch on demand
                const isNumerical =
                  dataset.allClusterColumnTypes?.[columnKey] === "numerical";

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
                  column.loaded ? (
                    <span className="w-2 h-2 rounded-full bg-green-500 flex-shrink-0" />
                  ) : null
                }
              >
                {column.label}
              </AutocompleteItem>
            ))}
          </Autocomplete>
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

            {/* Clear Button */}
            <Button
              className="w-full"
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
                      toggleCelltype(item.id);
                      // Mode is now automatically updated by toggleCelltype
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
  );
}
