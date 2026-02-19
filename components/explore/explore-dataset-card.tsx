"use client";

import { Card, CardBody } from "@heroui/card";
import { Chip } from "@heroui/chip";
import { Image } from "@heroui/image";
import { useRouter } from "next/navigation";

import type { CatalogDatasetItem } from "./types";

interface ExploreDatasetCardProps {
  dataset: CatalogDatasetItem;
}

export function ExploreDatasetCard({ dataset }: ExploreDatasetCardProps) {
  const router = useRouter();

  const handlePress = () => {
    if (dataset.externalLink) {
      window.open(dataset.externalLink, "_blank", "noopener,noreferrer");
      return;
    }
    if (dataset.s3BaseUrl) {
      const base = dataset.datasetType === "single_molecule" ? "/sm-viewer" : "/viewer";
      router.push(`${base}/from-s3?url=${encodeURIComponent(dataset.s3BaseUrl)}`);
      return;
    }
    if (dataset.datasetId) {
      const base = dataset.datasetType === "single_molecule" ? "/sm-viewer" : "/viewer";
      router.push(`${base}/${dataset.datasetId}`);
    }
  };

  const formatCount = (n: number | null) => {
    if (n == null) return null;
    if (n >= 1_000_000) return `${(n / 1_000_000).toFixed(1)}M`;
    if (n >= 1_000) return `${(n / 1_000).toFixed(0)}K`;
    return String(n);
  };

  return (
    <Card
      isBlurred
      isPressable
      className="border-none bg-background/60 dark:bg-default-100/50 hover:bg-default-200/50 hover:scale-[1.02] transition-all duration-200 cursor-pointer"
      shadow="sm"
      onPress={handlePress}
    >
      <CardBody className="p-0 overflow-hidden">
        {/* Thumbnail */}
        <div className="relative h-40 w-full bg-default-200/30">
          {dataset.thumbnailUrl ? (
            <Image
              alt={dataset.title}
              className="w-full h-40 object-cover"
              radius="none"
              src={dataset.thumbnailUrl}
            />
          ) : (
            <div className="w-full h-full flex items-center justify-center text-default-400">
              <span className="text-3xl">
                {dataset.datasetType === "single_cell" ? "🔬" : "🧬"}
              </span>
            </div>
          )}
          {/* Type badge overlay */}
          <div className="absolute top-2 left-2">
            <Chip
              color={dataset.datasetType === "single_cell" ? "primary" : "secondary"}
              size="sm"
              variant="solid"
            >
              {dataset.datasetType === "single_cell" ? "Single Cell" : "Single Molecule"}
            </Chip>
          </div>
        </div>

        {/* Content */}
        <div className="p-4 flex flex-col gap-2">
          <h3 className="font-semibold text-foreground text-sm line-clamp-2">
            {dataset.title}
          </h3>

          {/* Metadata row */}
          <div className="flex flex-wrap gap-1">
            {dataset.species && (
              <Chip size="sm" variant="flat">{dataset.species}</Chip>
            )}
            {dataset.tissue && (
              <Chip size="sm" variant="flat">{dataset.tissue}</Chip>
            )}
          </div>

          {/* Platform & stats */}
          <div className="flex flex-col gap-0.5 text-xs text-default-500">
            {dataset.platform && <span>Platform: {dataset.platform}</span>}
            <div className="flex gap-3">
              {dataset.numCells != null && (
                <span>
                  {formatCount(dataset.numCells)}{" "}
                  {dataset.datasetType === "single_molecule" ? "molecules" : "cells"}
                </span>
              )}
              {dataset.numGenes != null && (
                <span>{formatCount(dataset.numGenes)} genes</span>
              )}
            </div>
          </div>

          {/* Tags */}
          {dataset.tags.length > 0 && (
            <div className="flex flex-wrap gap-1 mt-1">
              {dataset.tags.slice(0, 3).map((tag) => (
                <Chip key={tag} className="text-[10px]" size="sm" variant="dot">
                  {tag}
                </Chip>
              ))}
            </div>
          )}
        </div>
      </CardBody>
    </Card>
  );
}
