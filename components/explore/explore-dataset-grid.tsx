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

  const rangeStart = (page - 1) * limit + 1;
  const rangeEnd = Math.min(page * limit, total);

  return (
    <div>
      <div className="flex items-center justify-between mb-4">
        <p className="text-sm text-default-500">
          Showing {rangeStart}–{rangeEnd} of {total} dataset{total !== 1 ? "s" : ""}
        </p>
        {totalPages > 1 && (
          <Pagination
            showControls
            page={page}
            size="sm"
            total={totalPages}
            onChange={onPageChange}
          />
        )}
      </div>

      <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 gap-6">
        {datasets.map((dataset) => (
          <ExploreDatasetCard key={dataset.id} dataset={dataset} />
        ))}
      </div>

      {totalPages > 1 && (
        <div className="flex items-center justify-between mt-8">
          <p className="text-sm text-default-500">
            Page {page} of {totalPages}
          </p>
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
