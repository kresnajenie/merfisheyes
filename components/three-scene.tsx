"use client";

import type { StandardizedDataset } from "@/lib/StandardizedDataset";
import type { PointData } from "@/lib/webgl/types";

import { useEffect, useRef, useState } from "react";
import * as THREE from "three";
import { Spinner } from "@heroui/react";

import { initializeScene } from "@/lib/webgl/scene-manager";
import {
  createPointCloud,
  updatePointCloudAttributes,
} from "@/lib/webgl/point-cloud";
import {
  updateGeneVisualization,
  updateCelltypeVisualization,
  updateNumericalCelltypeVisualization,
} from "@/lib/webgl/visualization-utils";
import { useVisualizationStore } from "@/lib/stores/visualizationStore";

interface ThreeSceneProps {
  dataset?: StandardizedDataset | null;
}

export function ThreeScene({ dataset }: ThreeSceneProps) {
  const containerRef = useRef<HTMLDivElement>(null);
  const pointCloudRef = useRef<THREE.Points | null>(null);
  const sceneRef = useRef<THREE.Scene | null>(null);
  const [isLoadingGene, setIsLoadingGene] = useState(false);

  // Get visualization settings from store
  const {
    mode,
    selectedGene,
    selectedColumn,
    selectedCelltypes,
    colorPalette,
    alphaScale,
    sizeScale,
  } = useVisualizationStore();

  // Effect 1: Scene initialization - runs when dataset changes
  useEffect(() => {
    if (!containerRef.current) return;

    // If dataset is provided, use its spatial coordinates
    if (dataset) {
      console.log("dataset");
      console.log("Creating point cloud from dataset:", dataset.name);
      console.log("Point count:", dataset.getPointCount());
      console.log("Spatial dimensions:", dataset.spatial.dimensions);

      // Calculate bounding box and center of spatial data
      const bounds = {
        minX: Infinity,
        maxX: -Infinity,
        minY: Infinity,
        maxY: -Infinity,
        minZ: Infinity,
        maxZ: -Infinity,
      };

      dataset.spatial.coordinates.forEach((coord) => {
        bounds.minX = Math.min(bounds.minX, coord[0]);
        bounds.maxX = Math.max(bounds.maxX, coord[0]);
        bounds.minY = Math.min(bounds.minY, coord[1]);
        bounds.maxY = Math.max(bounds.maxY, coord[1]);
        if (dataset.spatial.dimensions === 3 && coord[2] !== undefined) {
          bounds.minZ = Math.min(bounds.minZ, coord[2]);
          bounds.maxZ = Math.max(bounds.maxZ, coord[2]);
        }
      });

      // Calculate center point (scaled)
      const center = new THREE.Vector3(
        ((bounds.minX + bounds.maxX) / 2) * 500,
        ((bounds.minY + bounds.maxY) / 2) * 500,
        dataset.spatial.dimensions === 3
          ? ((bounds.minZ + bounds.maxZ) / 2) * 500
          : 0,
      );

      // Calculate size of data
      const size = Math.max(
        (bounds.maxX - bounds.minX) * 500,
        (bounds.maxY - bounds.minY) * 500,
        dataset.spatial.dimensions === 3
          ? (bounds.maxZ - bounds.minZ) * 500
          : 0,
      );

      // Position camera at appropriate distance
      const distance = size * 1.5;
      const cameraPos = new THREE.Vector3(
        center.x,
        center.y,
        center.z + distance,
      );

      // Initialize Three.js scene with options
      const { scene, animate, dispose } = initializeScene(
        containerRef.current,
        {
          is2D: dataset.spatial.dimensions === 2,
          cameraPosition: cameraPos,
          lookAtPosition: center,
        },
      );

      // Convert dataset spatial coordinates to PointData format
      const pointData: PointData[] = dataset.spatial.coordinates.map(
        (coord) => ({
          x: coord[0] * 500, // Scale coordinates
          y: coord[1] * 500,
          z: coord[2] !== undefined ? coord[2] * 500 : 0,
          r: Math.random(), // Random colors for now
          g: Math.random(),
          b: Math.random(),
          size: 1.0,
          alpha: 1.0,
        }),
      );

      console.log("Point data created:", pointData.length, "points");

      // Create point cloud mesh with custom shaders
      const pointCloud = createPointCloud(pointData, 5);

      pointCloudRef.current = pointCloud; // Store reference
      sceneRef.current = scene; // Store scene reference
      scene.add(pointCloud);

      // Start animation
      animate();

      // Cleanup on unmount
      return () => {
        scene.remove(pointCloud);
        pointCloud.geometry.dispose();
        (pointCloud.material as any).dispose();
        pointCloudRef.current = null;
        sceneRef.current = null;
        dispose();
      };
    } else {
      console.log("No dataset provided to ThreeScene");

      // Initialize scene without dataset
      const { animate, dispose } = initializeScene(containerRef.current);

      // Start animation even without data
      animate();

      return () => {
        dispose();
      };
    }
  }, [dataset]);

  // Effect 2: Update gene visualization
  useEffect(() => {
    if (!pointCloudRef.current || !dataset || mode !== "gene") return;

    console.log("Updating gene visualization:", selectedGene);

    const updateGene = async () => {
      if (!pointCloudRef.current || !dataset) return;

      // If no gene selected, fall back to celltype visualization
      if (!selectedGene) {
        setIsLoadingGene(false);
        const result = updateCelltypeVisualization(
          dataset,
          selectedColumn,
          selectedCelltypes,
          colorPalette,
          alphaScale,
          sizeScale,
        );

        if (result && pointCloudRef.current) {
          updatePointCloudAttributes(
            pointCloudRef.current,
            result.colors,
            result.sizes,
            result.alphas,
          );
        }

        return;
      }

      try {
        setIsLoadingGene(true);
        const result = await updateGeneVisualization(
          dataset,
          selectedGene,
          alphaScale,
          sizeScale,
        );

        if (result && pointCloudRef.current) {
          updatePointCloudAttributes(
            pointCloudRef.current,
            result.colors,
            result.sizes,
            result.alphas,
          );
        }
      } finally {
        setIsLoadingGene(false);
      }
    };

    updateGene();
  }, [
    dataset,
    selectedGene,
    selectedColumn,
    selectedCelltypes,
    colorPalette,
    alphaScale,
    sizeScale,
  ]);

  // Effect 3: Update celltype visualization
  useEffect(() => {
    if (!pointCloudRef.current || !dataset || mode !== "celltype") return;

    console.log("Updating celltype visualization:", {
      selectedColumn,
      selectedCelltypes: Array.from(selectedCelltypes),
    });

    // Check if the selected column is numerical
    const selectedCluster = dataset.clusters?.find(
      (c) => c.column === selectedColumn,
    );
    const isNumerical = selectedCluster?.type === "numerical";

    // Use appropriate visualization function based on column type
    const result = isNumerical
      ? updateNumericalCelltypeVisualization(
          dataset,
          selectedColumn,
          alphaScale,
          sizeScale,
        )
      : updateCelltypeVisualization(
          dataset,
          selectedColumn,
          selectedCelltypes,
          colorPalette,
          alphaScale,
          sizeScale,
        );

    if (result && pointCloudRef.current) {
      updatePointCloudAttributes(
        pointCloudRef.current,
        result.colors,
        result.sizes,
        result.alphas,
      );
    }
  }, [
    dataset,
    selectedColumn,
    selectedCelltypes,
    colorPalette,
    alphaScale,
    sizeScale,
  ]);

  return (
    <>
      <div
        ref={containerRef}
        className="fixed inset-0 w-full h-full"
        style={{ margin: 0, padding: 0 }}
      />

      {/* Loading overlay for gene expression fetching */}
      {isLoadingGene && (
        <div className="fixed inset-0 flex items-center justify-center bg-black/30 backdrop-blur-sm z-50 pointer-events-none">
          <div className="bg-default-100/90 rounded-lg p-6 shadow-lg flex flex-col items-center gap-3">
            <Spinner color="primary" size="lg" />
            <p className="text-sm font-medium">Loading gene expression...</p>
            {selectedGene && (
              <p className="text-xs text-default-500">{selectedGene}</p>
            )}
          </div>
        </div>
      )}
    </>
  );
}
