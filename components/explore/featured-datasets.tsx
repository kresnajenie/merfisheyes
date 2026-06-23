"use client";

import type { CatalogDatasetItem } from "./types";

import { ExploreDatasetCard } from "./explore-dataset-card";

interface FeaturedDatasetsProps {
  datasets: CatalogDatasetItem[];
  onViewAll?: () => void;
  onCardClick?: (dataset: CatalogDatasetItem) => void;
}

export function FeaturedDatasets({
  datasets,
  onViewAll,
  onCardClick,
}: FeaturedDatasetsProps) {
  if (datasets.length === 0) return null;

  return (
    <section>
      <h2 className="text-lg font-semibold mb-4">Featured Datasets</h2>
      <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 xl:grid-cols-4 gap-4">
        {datasets.slice(0, 4).map((dataset) => (
          <ExploreDatasetCard
            key={dataset.id}
            dataset={dataset}
            onCardClick={onCardClick}
          />
        ))}
      </div>
      {datasets.length > 4 && onViewAll && (
        <button
          className="mt-3 text-sm text-primary hover:underline"
          onClick={onViewAll}
        >
          View all {datasets.length} featured →
        </button>
      )}
    </section>
  );
}
