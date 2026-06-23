"use client";

import type { CatalogDatasetItem } from "./types";

import { ExploreDatasetCard } from "./explore-dataset-card";

interface BilDatasetsProps {
  datasets: CatalogDatasetItem[];
  onViewAll?: () => void;
  onCardClick?: (dataset: CatalogDatasetItem) => void;
}

export function BilDatasets({
  datasets,
  onViewAll,
  onCardClick,
}: BilDatasetsProps) {
  if (datasets.length === 0) return null;

  return (
    <section>
      <h2 className="text-lg font-semibold mb-4">From Brain Image Library</h2>
      <div className="grid grid-cols-2 gap-3">
        {datasets.slice(0, 5).map((dataset) => (
          <ExploreDatasetCard
            key={dataset.id}
            dataset={dataset}
            onCardClick={onCardClick}
          />
        ))}
      </div>
      {datasets.length > 5 && onViewAll && (
        <button
          className="mt-3 text-sm text-primary hover:underline"
          onClick={onViewAll}
        >
          View all {datasets.length} BIL datasets →
        </button>
      )}
    </section>
  );
}
