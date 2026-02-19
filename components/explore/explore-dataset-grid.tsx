"use client";

import { Pagination } from "@heroui/pagination";
import { Spinner } from "@heroui/spinner";

import { ExploreDatasetCard } from "./explore-dataset-card";
import type { CatalogDatasetItem } from "./types";

interface ExploreDatasetGridProps {
  datasets: CatalogDatasetItem[];
  loading: boolean;
  total: number;
  page: number;
  limit: number;
  onPageChange: (page: number) => void;
}

export function ExploreDatasetGrid({
  datasets,
  loading,
  total,
  page,
  limit,
  onPageChange,
}: ExploreDatasetGridProps) {
  const totalPages = Math.max(1, Math.ceil(total / limit));

  if (loading) {
    return (
      <div className="flex justify-center py-12">
        <Spinner size="lg" />
      </div>
    );
  }

  if (datasets.length === 0) {
    return (
      <p className="text-default-500 text-center py-12">
        No datasets found. Try adjusting your search or filters.
      </p>
    );
  }

  return (
    <div>
      <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 gap-6">
        {datasets.map((dataset) => (
          <ExploreDatasetCard key={dataset.id} dataset={dataset} />
        ))}
      </div>

      {totalPages > 1 && (
        <div className="flex justify-center mt-8">
          <Pagination
            showControls
            page={page}
            total={totalPages}
            onChange={onPageChange}
          />
        </div>
      )}
    </div>
  );
}
