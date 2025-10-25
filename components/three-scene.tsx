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
  updateCombinedVisualization,
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

  // Raycaster and interaction state
  const raycasterRef = useRef<THREE.Raycaster>(new THREE.Raycaster());
  const mouseRef = useRef<THREE.Vector2>(new THREE.Vector2());
  const cameraRef = useRef<THREE.Camera | null>(null);
  const rendererRef = useRef<THREE.WebGLRenderer | null>(null);
  const tooltipRef = useRef<HTMLDivElement | null>(null);
  const hoveredPointRef = useRef<number | null>(null);
  const lastCameraPositionRef = useRef<THREE.Vector3>(new THREE.Vector3());

  // Store current visualization data for tooltips
  const geneExpressionRef = useRef<number[] | null>(null);
  const colorPaletteRef = useRef<Record<string, string>>({});
  const clusterValuesRef = useRef<(string | number)[]>([]);
  const isNumericalClusterRef = useRef<boolean>(false);

  // Get visualization settings from store
  const {
    mode,
    selectedGene,
    selectedColumn,
    selectedCelltypes,
    colorPalette,
    alphaScale,
    sizeScale,
    toggleCelltype,
  } = useVisualizationStore();

  // Store current mode and selection in refs to avoid closure issues
  const modeRef = useRef(mode);
  const selectedGeneRef = useRef(selectedGene);
  const selectedColumnRef = useRef(selectedColumn);

  // Update refs when store values change
  useEffect(() => {
    modeRef.current = mode;
    selectedGeneRef.current = selectedGene;
    selectedColumnRef.current = selectedColumn;
  }, [mode, selectedGene, selectedColumn]);

  // Helper function: Create tooltip element
  const createTooltip = (): HTMLDivElement => {
    const tooltip = document.createElement("div");
    tooltip.className =
      "absolute bg-black/80 text-white px-2.5 py-1.5 rounded text-sm font-sans pointer-events-none hidden z-[1000] shadow-lg min-w-[80px]";
    document.body.appendChild(tooltip);
    return tooltip;
  };

  // Helper function: Get color from point cloud at specific index
  const getPointColor = (index: number): string => {
    if (!pointCloudRef.current) return "#808080";
    const colorAttr = pointCloudRef.current.geometry.attributes
      .color as THREE.BufferAttribute;
    if (!colorAttr) return "#808080";

    const r = Math.round(colorAttr.getX(index) * 255);
    const g = Math.round(colorAttr.getY(index) * 255);
    const b = Math.round(colorAttr.getZ(index) * 255);

    return `#${r.toString(16).padStart(2, "0")}${g.toString(16).padStart(2, "0")}${b.toString(16).padStart(2, "0")}`;
  };

  // Helper function: Show tooltip
  const showTooltip = (position: THREE.Vector3, index: number) => {
    if (!tooltipRef.current || !cameraRef.current || !rendererRef.current)
      return;

    // Convert 3D position to screen coordinates
    const vector = position.clone();
    vector.project(cameraRef.current);

    const x = (vector.x * 0.5 + 0.5) * rendererRef.current.domElement.clientWidth;
    const y =
      (-vector.y * 0.5 + 0.5) * rendererRef.current.domElement.clientHeight;

    // Get point color (gene gradient color in gene mode)
    const pointColor = getPointColor(index);

    // Determine what to show based on current mode and cluster type
    const isNumerical = isNumericalClusterRef.current;
    const clusterValue = clusterValuesRef.current[index];
    const geneValue = geneExpressionRef.current?.[index];

    // Get current values from refs
    const currentMode = modeRef.current;
    const currentGene = selectedGeneRef.current;
    const currentColumn = selectedColumnRef.current;

    let tooltipContent = "";

    if (currentMode.includes("gene") && currentGene) {
      // Gene mode
      if (isNumerical) {
        // Numerical cluster + gene: show both values without color circle
        tooltipContent = `
          <div class="flex flex-col gap-1">
            <div>${currentColumn}: ${clusterValue}</div>
            <div>${currentGene}: ${geneValue?.toFixed(2) ?? "N/A"}</div>
          </div>
        `;
      } else {
        // Categorical cluster + gene: show 2 rows with colored circles
        // Row 1: cluster color + cluster name
        // Row 2: gene gradient color + gene value
        const clusterColor =
          colorPaletteRef.current[String(clusterValue)] || "#808080";

        tooltipContent = `
          <div class="flex flex-col gap-1">
            <div class="flex items-center">
              <div style="width: 12px; height: 12px; border-radius: 50%; background-color: ${clusterColor}; margin-right: 6px;"></div>
              <span>${clusterValue}</span>
            </div>
            <div class="flex items-center">
              <div style="width: 12px; height: 12px; border-radius: 50%; background-color: ${pointColor}; margin-right: 6px;"></div>
              <span>${geneValue?.toFixed(2) ?? "N/A"}</span>
            </div>
          </div>
        `;
      }
    } else {
      // Celltype mode or gene mode without gene selected
      if (isNumerical) {
        // Numerical cluster: just show the value
        tooltipContent = `<div>${currentColumn}: ${clusterValue}</div>`;
      } else {
        // Categorical cluster: show with color circle from palette (not rendered color)
        const clusterColor =
          colorPaletteRef.current[String(clusterValue)] || "#808080";
        tooltipContent = `
          <div class="flex items-center">
            <div style="width: 12px; height: 12px; border-radius: 50%; background-color: ${clusterColor}; margin-right: 6px;"></div>
            <span>${clusterValue}</span>
          </div>
        `;
      }
    }

    tooltipRef.current.innerHTML = tooltipContent;
    tooltipRef.current.style.left = `${x + 10}px`;
    tooltipRef.current.style.top = `${y + 10}px`;
    tooltipRef.current.classList.remove("hidden");
  };

  // Helper function: Hide tooltip
  const hideTooltip = () => {
    if (tooltipRef.current) {
      tooltipRef.current.classList.add("hidden");
    }
  };

  // Helper function: Check intersections with adaptive threshold
  const checkIntersections = () => {
    if (
      !pointCloudRef.current ||
      !cameraRef.current ||
      !rendererRef.current ||
      !dataset
    )
      return;

    // Get current camera position to check if we've moved
    const currentCameraPosition = cameraRef.current.position.clone();
    const cameraHasMoved = !currentCameraPosition.equals(
      lastCameraPositionRef.current,
    );
    lastCameraPositionRef.current.copy(currentCameraPosition);

    // Calculate camera distance to determine raycaster parameters
    const cameraDistance = cameraRef.current.position.length();

    // Set adaptive thresholds for raycasting with multiple tiers
    // IMPORTANT: Smaller threshold = more precision needed, Larger threshold = easier selection
    // When zoomed IN (small distance) = need SMALLER threshold for accuracy
    // When zoomed OUT (large distance) = need LARGER threshold for easier selection
    // Note: threshold is in world space units, so it needs to be VERY small
    let threshold;
    if (cameraDistance < 150) {
      // Very close zoom: precise selection
      threshold = 0.1;
    } else if (cameraDistance < 250) {
      // Close zoom: moderately precise
      threshold = 0.2;
    } else if (cameraDistance < 400) {
      // Medium zoom: balanced (your current zoom at ~315)
      threshold = 0.3;
    } else if (cameraDistance < 600) {
      // Far zoom: easier selection
      threshold = 0.5;
    } else if (cameraDistance < 900) {
      // Very far zoom: very easy selection
      threshold = 1.0;
    } else {
      // Extremely far zoom: maximum ease
      threshold = 2.0;
    }

    raycasterRef.current.params.Points!.threshold = threshold;

    // Update the raycaster with the current mouse position and camera
    raycasterRef.current.setFromCamera(mouseRef.current, cameraRef.current);

    // Check for intersections with the points mesh
    const intersects = raycasterRef.current.intersectObject(
      pointCloudRef.current,
    );

    // If we found an intersection
    if (intersects.length > 0) {
      // Sort intersections by distance if there are multiple
      if (intersects.length > 1) {
        intersects.sort((a, b) => a.distance - b.distance);
      }

      // Get the index of the closest point that was intersected
      const index = intersects[0].index!;

      // If this is a different point than the one we were previously hovering over
      if (hoveredPointRef.current !== index) {
        hoveredPointRef.current = index;
        const position = intersects[0].point;
        showTooltip(position, index);
      }
    } else {
      // If we're not hovering over any point, hide the tooltip
      if (hoveredPointRef.current !== null) {
        hoveredPointRef.current = null;
        hideTooltip();
      }
    }
  };

  // Effect 0: Tooltip cleanup on component unmount
  useEffect(() => {
    return () => {
      if (tooltipRef.current) {
        document.body.removeChild(tooltipRef.current);
        tooltipRef.current = null;
      }
    };
  }, []);

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
      const { scene, camera, renderer, animate, dispose } = initializeScene(
        containerRef.current,
        {
          is2D: dataset.spatial.dimensions === 2,
          cameraPosition: cameraPos,
          lookAtPosition: center,
        },
      );

      // Store camera and renderer refs for raycasting
      cameraRef.current = camera;
      rendererRef.current = renderer;

      // Create tooltip
      if (!tooltipRef.current) {
        tooltipRef.current = createTooltip();
      }

      // Throttled intersection checking (50ms)
      let lastCheckTime = 0;
      const throttleDelay = 50;

      // Mouse move handler
      const handleMouseMove = (event: MouseEvent) => {
        const rect = renderer.domElement.getBoundingClientRect();
        mouseRef.current.x = ((event.clientX - rect.left) / rect.width) * 2 - 1;
        mouseRef.current.y =
          -((event.clientY - rect.top) / rect.height) * 2 + 1;

        // Throttle intersection checking
        const now = Date.now();
        if (now - lastCheckTime > throttleDelay) {
          lastCheckTime = now;
          checkIntersections();
        }
      };

      // Double-click handler
      const handleDoubleClick = () => {
        if (hoveredPointRef.current !== null && dataset.clusters) {
          const index = hoveredPointRef.current;
          const clusterValue = clusterValuesRef.current[index];
          const clusterValueStr = String(clusterValue);

          console.log("Double-clicked cluster:", {
            column: selectedColumnRef.current,
            value: clusterValue,
            index: index,
          });

          // Only toggle celltype if it's not a numerical cluster
          if (!isNumericalClusterRef.current) {
            // Toggle the cluster in selectedCelltypes
            toggleCelltype(clusterValueStr);
            // Mode is now automatically updated by toggleCelltype
          }
        }
      };

      // Add event listeners
      renderer.domElement.addEventListener("mousemove", handleMouseMove);
      renderer.domElement.addEventListener("dblclick", handleDoubleClick);

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
        // Remove event listeners
        renderer.domElement.removeEventListener("mousemove", handleMouseMove);
        renderer.domElement.removeEventListener("dblclick", handleDoubleClick);

        // Hide tooltip
        hideTooltip();

        scene.remove(pointCloud);
        pointCloud.geometry.dispose();
        (pointCloud.material as any).dispose();
        pointCloudRef.current = null;
        sceneRef.current = null;
        cameraRef.current = null;
        rendererRef.current = null;
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

  // Effect 2: Update visualization based on mode array
  useEffect(() => {
    if (!pointCloudRef.current || !dataset) return;

    console.log("Updating visualization with mode:", mode);

    const updateVisualization = async () => {
      if (!pointCloudRef.current || !dataset) return;

      // Store cluster data for tooltip
      const selectedCluster = dataset.clusters?.find(
        (c) => c.column === selectedColumn,
      );
      if (selectedCluster) {
        clusterValuesRef.current = selectedCluster.values;
        isNumericalClusterRef.current = selectedCluster.type === "numerical";
        colorPaletteRef.current = selectedCluster.palette || colorPalette;
      }

      // Determine which visualization to use based on mode array
      const hasGeneMode = mode.includes("gene");
      const hasCelltypeMode = mode.includes("celltype");

      // Check if the selected column is numerical
      const isNumerical = selectedCluster?.type === "numerical";

      let result = null;

      if (hasGeneMode && hasCelltypeMode && selectedGene && selectedCelltypes.size > 0) {
        // Combined mode: gene expression on selected celltypes
        console.log("Using combined gene + celltype visualization");

        try {
          setIsLoadingGene(true);

          // Fetch gene expression data for tooltip
          const expression = await dataset.getGeneExpression(selectedGene);
          geneExpressionRef.current = expression;

          result = await updateCombinedVisualization(
            dataset,
            selectedGene,
            selectedColumn,
            selectedCelltypes,
            alphaScale,
            sizeScale,
          );
        } finally {
          setIsLoadingGene(false);
        }
      } else if (hasGeneMode && selectedGene) {
        // Gene mode only
        console.log("Using gene-only visualization");

        try {
          setIsLoadingGene(true);

          // Fetch gene expression data for tooltip
          const expression = await dataset.getGeneExpression(selectedGene);
          geneExpressionRef.current = expression;

          result = await updateGeneVisualization(
            dataset,
            selectedGene,
            alphaScale,
            sizeScale,
          );
        } finally {
          setIsLoadingGene(false);
        }
      } else if (hasCelltypeMode) {
        // Celltype mode only
        console.log("Using celltype-only visualization");
        setIsLoadingGene(false);

        // Use appropriate visualization function based on column type
        result = isNumerical
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
      }

      // Update point cloud if we have a result
      if (result && pointCloudRef.current) {
        updatePointCloudAttributes(
          pointCloudRef.current,
          result.colors,
          result.sizes,
          result.alphas,
        );
      }
    };

    updateVisualization();
  }, [
    dataset,
    selectedGene,
    selectedColumn,
    selectedCelltypes,
    colorPalette,
    alphaScale,
    sizeScale,
    mode,
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
