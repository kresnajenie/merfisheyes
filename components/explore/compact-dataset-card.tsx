"use client";

import { Card, CardBody } from "@heroui/card";
import { Chip } from "@heroui/chip";

import type { CatalogDatasetItem } from "./types";

interface CompactDatasetCardProps {
  dataset: CatalogDatasetItem;
  onSelect: (dataset: CatalogDatasetItem) => void;
}

export function CompactDatasetCard({ dataset, onSelect }: CompactDatasetCardProps) {
  const formatCount = (n: number | null) => {
    if (n == null) return null;
    if (n >= 1_000_000) return `${(n / 1_000_000).toFixed(1)}M`;
    if (n >= 1_000) return `${(n / 1_000).toFixed(0)}K`;
    return String(n);
  };

  return (
    <Card
      isPressable
      className="border-none bg-default-100/50 hover:bg-default-200/50 transition-colors"
      shadow="none"
      onPress={() => onSelect(dataset)}
    >
      <CardBody className="p-3 flex flex-row gap-3 items-center">
        {/* Thumbnail */}
        <div className="w-10 h-10 rounded-lg bg-default-200/50 flex items-center justify-center shrink-0 overflow-hidden">
          {dataset.thumbnailUrl ? (
            <img
              alt=""
              className="w-full h-full object-cover"
              src={dataset.thumbnailUrl}
            />
          ) : (
            <span className="text-lg">
              {dataset.datasetType === "single_cell" ? "🔬" : "🧬"}
            </span>
          )}
        </div>

        {/* Info */}
        <div className="flex-1 min-w-0">
          <p className="text-sm font-medium truncate">{dataset.title}</p>
          <div className="flex items-center gap-2 text-xs text-default-500">
            <Chip
              color={dataset.datasetType === "single_cell" ? "primary" : "secondary"}
              size="sm"
              variant="flat"
            >
              {dataset.datasetType === "single_cell" ? "SC" : "SM"}
            </Chip>
            {dataset.numCells != null && (
              <span>{formatCount(dataset.numCells)} cells</span>
            )}
          </div>
        </div>
      </CardBody>
    </Card>
  );
}
