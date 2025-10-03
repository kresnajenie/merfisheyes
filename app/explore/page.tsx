"use client";

import { title } from "@/components/primitives";
import { DatasetCard } from "@/components/dataset-card";
import { exampleDatasets } from "@/lib/data/example-datasets";
import LightRays from "@/components/react-bits/LightRays";

export default function ExplorePage() {
  return (
    <>
      <div className="fixed inset-0 w-full h-full z-0">
        <LightRays
          raysOrigin="top-left"
          raysColor="#FF1CF7"
          rayLength={10}
          raysSpeed={1.0}
          lightSpread={1.2}
          pulsating={false}
          mouseInfluence={0.1}
        />
      </div>
      <div className="w-full">
        <div className="mb-8">
          <h1 className={title({ size: "lg" })}>
            <span className={title({ color: "violet", size: "lg" })}>
              Explore
            </span>{" "}
            Example Datasets
          </h1>
          <p className="text-default-500 mt-4">
            Browse and visualize our curated collection of spatial
            transcriptomics datasets
          </p>
        </div>

        <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 gap-6">
          {exampleDatasets.map((dataset) => (
            <DatasetCard key={dataset.id} dataset={dataset} />
          ))}
        </div>
      </div>
    </>
  );
}
