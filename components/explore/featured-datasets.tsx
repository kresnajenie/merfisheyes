"use client";

import { ScrollShadow } from "@heroui/scroll-shadow";

import { ExploreDatasetCard } from "./explore-dataset-card";
import type { CatalogDatasetItem } from "./types";

interface FeaturedDatasetsProps {
  datasets: CatalogDatasetItem[];
}

export function FeaturedDatasets({ datasets }: FeaturedDatasetsProps) {
  if (datasets.length === 0) return null;

  return (
    <section className="mb-8">
      <h2 className="text-lg font-semibold mb-4">Featured Datasets</h2>
      <ScrollShadow orientation="horizontal" className="pb-4">
        <div className="flex gap-4">
          {datasets.map((dataset) => (
            <div key={dataset.id} className="w-72 shrink-0">
              <ExploreDatasetCard dataset={dataset} usePopover />
            </div>
          ))}
        </div>
      </ScrollShadow>
    </section>
  );
}
