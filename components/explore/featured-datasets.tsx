"use client";

import { ExploreDatasetCard } from "./explore-dataset-card";
import type { CatalogDatasetItem } from "./types";

interface FeaturedDatasetsProps {
  datasets: CatalogDatasetItem[];
  onViewAll?: () => void;
}

export function FeaturedDatasets({ datasets, onViewAll }: FeaturedDatasetsProps) {
  if (datasets.length === 0) return null;

  return (
    <section>
      <h2 className="text-lg font-semibold mb-4">Featured Datasets</h2>
      <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 xl:grid-cols-4 gap-4">
        {datasets.slice(0, 4).map((dataset) => (
          <ExploreDatasetCard key={dataset.id} dataset={dataset} />
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
