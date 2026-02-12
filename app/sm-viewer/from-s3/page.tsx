"use client";

import type { SingleMoleculeDataset } from "@/lib/SingleMoleculeDataset";

import { Suspense, useEffect, useState } from "react";
import { useRouter, useSearchParams } from "next/navigation";
import { Button, Spinner } from "@heroui/react";

import { SingleMoleculeThreeScene } from "@/components/single-molecule-three-scene";
import { SingleMoleculeControls } from "@/components/single-molecule-controls";
import { SingleMoleculeLegends } from "@/components/single-molecule-legends";
import { useSingleMoleculeStore } from "@/lib/stores/singleMoleculeStore";
import { useSingleMoleculeVisualizationStore } from "@/lib/stores/singleMoleculeVisualizationStore";
import LightRays from "@/components/react-bits/LightRays";
import { subtitle, title } from "@/components/primitives";

function SingleMoleculeViewerFromS3Content() {
  const searchParams = useSearchParams();
  const router = useRouter();
  const { addDataset } = useSingleMoleculeStore();
  const { addGene } = useSingleMoleculeVisualizationStore();
  const [dataset, setDataset] = useState<SingleMoleculeDataset | null>(null);
  const [isLoading, setIsLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  const s3Url = searchParams.get("url");

  useEffect(() => {
    if (!s3Url) {
      setError("No S3 URL provided");
      setIsLoading(false);

      return;
    }

    loadDatasetFromCustomS3(decodeURIComponent(s3Url));
  }, [s3Url]);

  const loadDatasetFromCustomS3 = async (baseUrl: string) => {
    try {
      setIsLoading(true);
      setError(null);

      console.log("Loading single molecule dataset from custom S3:", baseUrl);

      // Import SingleMoleculeDataset
      const { SingleMoleculeDataset } = await import(
        "@/lib/SingleMoleculeDataset"
      );

      // Load dataset using fromCustomS3 method with lazy loading
      const smDataset = await SingleMoleculeDataset.fromCustomS3(
        baseUrl,
        (progress, message) => {
          console.log(`${progress}%: ${message}`);
        },
      );

      console.log(
        "SingleMoleculeDataset loaded from custom S3:",
        smDataset.getSummary(),
      );

      // Store dataset in both local state and global store
      setDataset(smDataset);
      addDataset(smDataset);
      console.log("Dataset added to singleMoleculeStore");

      // Auto-select first 3 genes
      const genesToSelect = smDataset.uniqueGenes.slice(0, 3);

      console.log("Auto-selecting genes:", genesToSelect);

      genesToSelect.forEach((gene) => {
        const geneProps = smDataset.geneColors[gene];

        if (!geneProps) {
          console.error(`Missing geneProps for gene: ${gene}`);

          return;
        }

        console.log(
          `Adding gene to visualization: ${gene} with color ${geneProps.color}`,
        );
        addGene(gene, geneProps.color, geneProps.size);
      });

      setIsLoading(false);
    } catch (err) {
      console.error(
        "Error loading single molecule dataset from custom S3:",
        err,
      );
      setError(
        err instanceof Error
          ? err.message
          : "Failed to load dataset from custom S3",
      );
      setIsLoading(false);
    }
  };

  // Loading state
  if (isLoading) {
    return (
      <>
        <div className="fixed inset-0 w-full h-full z-0">
          <LightRays
            lightSpread={1.0}
            mouseInfluence={0.1}
            pulsating={false}
            rayLength={10}
            raysColor="#667eea"
            raysOrigin="top-left"
            raysSpeed={0.8}
          />
        </div>
        <div className="fixed inset-0 w-full h-full z-0">
          <LightRays
            lightSpread={1.0}
            mouseInfluence={0.1}
            pulsating={false}
            rayLength={10}
            raysColor="#764ba2"
            raysOrigin="top-right"
            raysSpeed={0.8}
          />
        </div>
        <div className="relative z-10 flex items-center justify-center h-full">
          <div className="flex flex-col items-center gap-4">
            <Spinner color="secondary" size="lg" />
            <p className={subtitle()}>Loading dataset from S3...</p>
          </div>
        </div>
      </>
    );
  }

  // Error state
  if (error) {
    return (
      <>
        <div className="fixed inset-0 w-full h-full z-0">
          <LightRays
            lightSpread={1.0}
            mouseInfluence={0.1}
            pulsating={false}
            rayLength={10}
            raysColor="#FF72E1"
            raysOrigin="top-center"
            raysSpeed={0.8}
          />
        </div>
        <div className="relative z-10 flex items-center justify-center h-full p-8">
          <div className="flex flex-col items-center gap-6 max-w-2xl w-full">
            <div className="text-center">
              <h2 className={title({ size: "md", color: "pink" })}>
                Failed to load dataset from S3
              </h2>
              <p className={subtitle({ class: "mt-4" })}>{error}</p>
            </div>
            <div className="flex gap-4">
              <Button
                color="secondary"
                onPress={() => router.push("/sm-viewer")}
              >
                Go to Home
              </Button>
              <Button
                color="default"
                variant="bordered"
                onPress={() =>
                  s3Url && loadDatasetFromCustomS3(decodeURIComponent(s3Url))
                }
              >
                Retry
              </Button>
            </div>
          </div>
        </div>
      </>
    );
  }

  // Dataset loaded
  if (!dataset) {
    return null;
  }

  return (
    <>
      <SingleMoleculeControls />
      <SingleMoleculeLegends />
      <SingleMoleculeThreeScene />
    </>
  );
}

export default function SingleMoleculeViewerFromS3Page() {
  return (
    <Suspense
      fallback={
        <div className="flex items-center justify-center h-full">
          <Spinner size="lg" />
        </div>
      }
    >
      <SingleMoleculeViewerFromS3Content />
    </Suspense>
  );
}
