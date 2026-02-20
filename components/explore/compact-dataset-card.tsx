"use client";

import { useState } from "react";
import { Card, CardBody } from "@heroui/card";
import { Chip } from "@heroui/chip";

import type { CatalogDatasetItem, CatalogDatasetEntry } from "./types";

interface CompactDatasetCardProps {
  dataset: CatalogDatasetItem;
  onSelect: (dataset: CatalogDatasetItem, entry: CatalogDatasetEntry) => void;
}

export function CompactDatasetCard({ dataset, onSelect }: CompactDatasetCardProps) {
  const [expanded, setExpanded] = useState(false);
  const entries = dataset.entries ?? [];
  const hasMultiple = entries.length > 1;
  const uniqueTypes = [...new Set(entries.map((e) => e.datasetType))];

  const formatCount = (n: number | null) => {
    if (n == null) return null;
    if (n >= 1_000_000) return `${(n / 1_000_000).toFixed(1)}M`;
    if (n >= 1_000) return `${(n / 1_000).toFixed(0)}K`;
    return String(n);
  };

  const handlePress = () => {
    if (entries.length === 1) {
      onSelect(dataset, entries[0]);
      return;
    }
    if (hasMultiple) {
      setExpanded((p) => !p);
    }
  };

  return (
    <div>
      <Card
        isPressable={entries.length > 0}
        className="border-none bg-default-100/50 hover:bg-default-200/50 transition-colors"
        shadow="none"
        onPress={handlePress}
      >
        <CardBody className="p-3 flex flex-row gap-3 items-center">
          {/* Thumbnail — only when image exists */}
          {dataset.thumbnailUrl && (
            <div className="w-10 h-10 rounded-lg bg-default-200/50 flex items-center justify-center shrink-0 overflow-hidden">
              <img
                alt=""
                className="w-full h-full object-cover"
                src={dataset.thumbnailUrl}
              />
            </div>
          )}

          {/* Info */}
          <div className="flex-1 min-w-0">
            <p className="text-sm font-medium truncate">{dataset.title}</p>
            <div className="flex items-center gap-2 text-xs text-default-500">
              {uniqueTypes.map((type) => (
                <Chip
                  key={type}
                  color={type === "single_cell" ? "primary" : "secondary"}
                  size="sm"
                  variant="flat"
                >
                  {type === "single_cell" ? "SC" : "SM"}
                </Chip>
              ))}
              {hasMultiple && (
                <span className="text-default-400">{entries.length} entries</span>
              )}
              {dataset.numCells != null && (
                <span>{formatCount(dataset.numCells)} cells</span>
              )}
            </div>
          </div>
        </CardBody>
      </Card>

      {/* Expanded entries */}
      {expanded && hasMultiple && (
        <div className="ml-6 mt-1 flex flex-col gap-1">
          {entries.map((entry) => (
            <button
              key={entry.id}
              className="flex items-center gap-2 px-3 py-1.5 rounded-lg text-left text-xs hover:bg-default-100 dark:hover:bg-default-200/50 transition-colors"
              onClick={() => onSelect(dataset, entry)}
            >
              <Chip
                color={entry.datasetType === "single_cell" ? "primary" : "secondary"}
                size="sm"
                variant="flat"
              >
                {entry.datasetType === "single_cell" ? "SC" : "SM"}
              </Chip>
              <span className="truncate">{entry.label}</span>
            </button>
          ))}
        </div>
      )}
    </div>
  );
}
