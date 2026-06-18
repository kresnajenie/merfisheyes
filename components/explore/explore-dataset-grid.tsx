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
  geneHighlight?: string;
}

export function ExploreDatasetGrid({
  datasets,
  loading,
  total,
  page,
  limit,
  onPageChange,
  geneHighlight,
}: ExploreDatasetGridProps) {
  const totalPages = Math.max(1, Math.ceil(total / limit));

  if (!loading && datasets.length === 0) {
    return (
      <p className="text-default-500 text-center py-12">
        No datasets found. Try adjusting your search or filters.
      </p>
    );
  }

  const rangeStart = (page - 1) * limit + 1;
  const rangeEnd = Math.min(page * limit, total);

  return (
    <div className="relative">
      {/* Loading overlay — keeps content visible */}
      {loading && (
        <div className="absolute inset-0 z-10 flex items-start justify-center pt-12 bg-background/50 rounded-xl">
          <Spinner size="lg" />
        </div>
      )}

      <div className={loading ? "opacity-50 pointer-events-none" : ""}>
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

        <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 xl:grid-cols-4 gap-4">
          {datasets.map((dataset) => (
            <ExploreDatasetCard key={dataset.id} dataset={dataset} geneHighlight={geneHighlight} />
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
    </div>
  );
}
