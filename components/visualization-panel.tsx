"use client";

import type { VisualizationMode } from "@/lib/stores/visualizationStore";
import type { StandardizedDataset } from "@/lib/StandardizedDataset";

import { Input } from "@heroui/input";
import { Button } from "@heroui/button";
import { Checkbox } from "@heroui/checkbox";
import { RadioGroup, Radio } from "@heroui/radio";
import { useState, useMemo } from "react";
import { Autocomplete, AutocompleteItem } from "@heroui/autocomplete";

import { useDatasetStore } from "@/lib/stores/datasetStore";
import { useVisualizationStore } from "@/lib/stores/visualizationStore";

interface VisualizationPanelProps {
  mode: VisualizationMode; // This is panelMode, not the visualization mode array
  onClose: () => void;
}

export function VisualizationPanel({ mode }: VisualizationPanelProps) {
  const [searchTerm, setSearchTerm] = useState("");
  const { getCurrentDataset } = useDatasetStore();
  const {
    selectedColumn,
    setSelectedColumn,
    selectedCelltypes,
    selectedGene,
    toggleCelltype,
    setSelectedGene,
    setMode,
  } = useVisualizationStore();

  // Type guard: this component only works with StandardizedDataset
  const rawDataset = getCurrentDataset();
  const dataset =
    rawDataset && "clusters" in rawDataset
      ? (rawDataset as StandardizedDataset)
      : null;

  // Get cluster columns for dropdown
  const clusterColumns = useMemo(() => {
    if (!dataset?.clusters) return [];

    return dataset.clusters.map((cluster) => ({
      key: cluster.column,
      label: `${cluster.column} (${cluster.type})`,
    }));
  }, [dataset]);

  // Check if the selected column is numerical
  const isNumericalColumn = useMemo(() => {
    if (!dataset?.clusters || !selectedColumn) return false;

    const selectedCluster = dataset.clusters.find(
      (c) => c.column === selectedColumn,
    );

    return selectedCluster?.type === "numerical";
  }, [dataset, selectedColumn]);

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

        return Array.from(itemsMap.values());
      }

      case "gene": {
        // Get all genes
        return dataset.genes.map((gene) => ({
          id: gene,
          label: gene,
          color: "#FFFFFF",
        }));
      }

      default:
        return [];
    }
  }, [dataset, mode, selectedColumn]);

  const filteredItems = items.filter((item) =>
    item.label.toLowerCase().includes(searchTerm.toLowerCase()),
  );

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
    <div className="fixed top-20 left-20 z-50 w-[300px] bg-content1 border-2 border-default-200 rounded-lg shadow-lg">
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
            onSelectionChange={(key) => {
              const columnKey = (key as string) || null;
              const selectedCluster = dataset?.clusters?.find(
                (c) => c.column === columnKey,
              );
              const isNumerical = selectedCluster?.type === "numerical";

              setSelectedColumn(columnKey, isNumerical);
              // When a cluster column is selected, switch visualization mode to celltype
              if (key) {
                setMode(["celltype"]);
              }
            }}
          >
            {clusterColumns.map((column) => (
              <AutocompleteItem key={column.key}>
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
              value={searchTerm}
              onValueChange={setSearchTerm}
            />

            {/* Clear Button */}
            <Button
              className="w-full"
              color="danger"
              variant="ghost"
              onPress={() => {
                setSearchTerm("");
                if (mode === "celltype") {
                  // Clear all selected celltypes
                  useVisualizationStore.setState({
                    selectedCelltypes: new Set<string>(),
                  });
                } else if (mode === "gene") {
                  // Clear selected gene
                  setSelectedGene(null);
                }
              }}
            >
              Clear
            </Button>

            {/* List */}
            <div className="max-h-[400px] overflow-y-auto flex flex-col gap-0">
              {mode === "celltype" ? (
                // Checkbox list for celltype mode
                filteredItems.map((item) => (
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
                  {filteredItems.map((item) => (
                    <Radio key={item.id} size="sm" value={item.id}>
                      <span style={{ color: item.color }}>{item.label}</span>
                    </Radio>
                  ))}
                </RadioGroup>
              ) : (
                // For other modes, just display items
                filteredItems.map((item) => (
                  <div key={item.id} className="p-2">
                    <span style={{ color: item.color }}>{item.label}</span>
                  </div>
                ))
              )}
            </div>
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
