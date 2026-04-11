"use client";

import type { StandardizedDataset } from "@/lib/StandardizedDataset";
import { getClusterValue } from "@/lib/StandardizedDataset";

import { useEffect, useRef, useState } from "react";
import * as THREE from "three";
import { toast } from "react-toastify";

import { initializeScene } from "@/lib/webgl/scene-manager";
import {
  createPointCloudFromBuffers,
  updatePointCloudAttributes,
  updateDotSize,
} from "@/lib/webgl/point-cloud";
import {
  updateGeneVisualization,
  updateCelltypeVisualization,
  updateNumericalCelltypeVisualization,
  updateCombinedVisualization,
} from "@/lib/webgl/visualization-utils";
import {
  usePanelVisualizationStore,
  usePanelId,
} from "@/lib/hooks/usePanelStores";
import { useSplitScreenStore } from "@/lib/stores/splitScreenStore";
import {
  getDatasetLinkConfig,
  fetchMappingConfig,
} from "@/lib/config/dataset-links";
import { VisualizationLegends } from "@/components/visualization-legends";
import { getEffectiveColumnType } from "@/lib/utils/column-type-utils";
import { SpatialScaleBar } from "@/components/spatial-scale-bar";
import { VISUALIZATION_CONFIG } from "@/lib/config/visualization.config";

interface ThreeSceneProps {
  dataset?: StandardizedDataset | null;
}

export function ThreeScene({ dataset }: ThreeSceneProps) {
  const containerRef = useRef<HTMLDivElement>(null);
  const pointCloudRef = useRef<THREE.Points | null>(null);
  const [pointCloudVersion, setPointCloudVersion] = useState(0);
  const sceneRef = useRef<THREE.Scene | null>(null);
  const geneToastIdRef = useRef<string | number | null>(null);

  // Raycaster and interaction state
  const raycasterRef = useRef<THREE.Raycaster>(new THREE.Raycaster());
  const mouseRef = useRef<THREE.Vector2>(new THREE.Vector2());
  const cameraRef = useRef<THREE.Camera | null>(null);
  const rendererRef = useRef<THREE.WebGLRenderer | null>(null);
  const controlsRef = useRef<any>(null);
  const tooltipRef = useRef<HTMLDivElement | null>(null);
  const hoveredPointRef = useRef<number | null>(null);
  const lastCameraPositionRef = useRef<THREE.Vector3>(new THREE.Vector3());

  // Store current visualization data for tooltips
  const geneExpressionRef = useRef<number[] | null>(null);
  const colorPaletteRef = useRef<Record<string, string>>({});
  const clusterRef = useRef<any>(null);
  const isNumericalClusterRef = useRef<boolean>(false);
  const baseDotSizeRef = useRef<number>(5);

  // Get visualization settings from store
  const {
    mode,
    selectedGene,
    selectedColumn,
    selectedCelltypes,
    colorPalette,
    alphaScale,
    sizeScale,
    geneScaleMin,
    geneScaleMax,
    setGeneScaleMin,
    setGeneScaleMax,
    numericalScaleMin,
    numericalScaleMax,
    setNumericalScaleMin,
    setNumericalScaleMax,
    toggleCelltype,
    clusterVersion,
    incrementClusterVersion,
    columnTypeOverrides,
    viewMode,
  } = usePanelVisualizationStore();

  // Split screen support
  const panelId = usePanelId();
  const { enableSplit, setRightPanelS3 } = useSplitScreenStore();
  const linkConfigRef = useRef(dataset ? getDatasetLinkConfig(dataset) : null);

  // Update link config when dataset changes:
  // 1. Check hardcoded registry first (instant)
  // 2. Try fetching mapping.json from custom S3 (async, background)
  // 3. If mapping found, lazy-load the link column so right-click works immediately
  useEffect(() => {
    if (!dataset) {
      linkConfigRef.current = null;
      return;
    }

    // Check hardcoded registry first
    const registryConfig = getDatasetLinkConfig(dataset);
    linkConfigRef.current = registryConfig;

    // Try fetching mapping.json (only for custom S3 datasets)
    if (!registryConfig && dataset.metadata?.customS3BaseUrl) {
      fetchMappingConfig(dataset).then(async (mappingConfig) => {
        if (!mappingConfig) return;

        linkConfigRef.current = mappingConfig;

        // "__all__" means every cell links to the same SM dataset — no column to preload
        if (mappingConfig.linkColumn === "__all__") return;

        // Check if link column is already loaded
        const alreadyLoaded = dataset.clusters?.some(
          (c) => c.column === mappingConfig.linkColumn,
        );
        if (alreadyLoaded) return;

        // Check if the column exists in the dataset
        if (
          dataset.allClusterColumnNames &&
          !dataset.allClusterColumnNames.includes(mappingConfig.linkColumn)
        ) {
          console.warn(
            `mapping.json linkColumn "${mappingConfig.linkColumn}" not found in dataset columns`,
          );
          return;
        }

        // Lazy-load the link column in the background
        try {
          const { getStandardizedDatasetWorker } = await import(
            "@/lib/workers/standardizedDatasetWorkerManager"
          );
          const worker = await getStandardizedDatasetWorker();

          const newClusters = await worker.loadClusterFromS3(
            dataset.id,
            [mappingConfig.linkColumn],
            dataset.metadata?.customS3BaseUrl,
          );

          if (newClusters && newClusters.length > 0) {
            dataset.addClusters(newClusters);
            incrementClusterVersion();
          }
        } catch (error) {
          console.warn(
            `Failed to preload link column "${mappingConfig.linkColumn}":`,
            error,
          );
        }
      });
    }
  }, [dataset, incrementClusterVersion]);

  // Store current mode and selection in refs to avoid closure issues
  const modeRef = useRef(mode);
  const selectedGeneRef = useRef(selectedGene);
  const selectedColumnRef = useRef(selectedColumn);
  const previousGeneRef = useRef<string | null>(null);
  const previousColumnRef = useRef<string | null>(null);
  const previousColumnTypeRef = useRef<boolean>(false);

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

    // Get canvas position relative to screen coordinates
    const rect = rendererRef.current.domElement.getBoundingClientRect();

    const x = (vector.x * 0.5 + 0.5) * rect.width + rect.left;
    const y = (-vector.y * 0.5 + 0.5) * rect.height + rect.top;

    // Get point color (gene gradient color in gene mode)
    const pointColor = getPointColor(index);

    // Determine what to show based on current mode and cluster type
    const isNumerical = isNumericalClusterRef.current;
    const clusterValue = clusterRef.current
      ? getClusterValue(clusterRef.current, index)
      : undefined;
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
    if (!pointCloudRef.current || !cameraRef.current || !rendererRef.current)
      return;

    // Get current camera position to check if we've moved
    const currentCameraPosition = cameraRef.current.position.clone();
    const cameraHasMoved = !currentCameraPosition.equals(
      lastCameraPositionRef.current,
    );

    lastCameraPositionRef.current.copy(currentCameraPosition);

    // Screen-space pixel threshold: convert a fixed pixel radius to world units
    // This adapts automatically to any coordinate range and zoom level
    const PIXEL_RADIUS = 5; // hover within 5 pixels of a point
    const camera = cameraRef.current as THREE.PerspectiveCamera;
    const target = controlsRef.current?.target;
    const cameraDistance = target
      ? camera.position.distanceTo(target)
      : camera.position.length();
    const canvasHeight = rendererRef.current.domElement.clientHeight;
    const fovRad = (camera.fov / 2) * (Math.PI / 180);
    const pixelSize = (2 * cameraDistance * Math.tan(fovRad)) / canvasHeight;
    const threshold = pixelSize * PIXEL_RADIUS;

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
      // Calculate bounding box and center of spatial data
      const bounds = {
        minX: Infinity,
        maxX: -Infinity,
        minY: Infinity,
        maxY: -Infinity,
        minZ: Infinity,
        maxZ: -Infinity,
      };

      const coords = dataset.spatial.coordinates;
      const dims = dataset.spatial.dimensions;
      const numPts = dataset.getPointCount();

      for (let p = 0; p < numPts; p++) {
        const x =
          coords instanceof Float32Array
            ? coords[p * dims]
            : (coords as number[][])[p][0];
        const y =
          coords instanceof Float32Array
            ? coords[p * dims + 1]
            : (coords as number[][])[p][1];

        bounds.minX = Math.min(bounds.minX, x);
        bounds.maxX = Math.max(bounds.maxX, x);
        bounds.minY = Math.min(bounds.minY, y);
        bounds.maxY = Math.max(bounds.maxY, y);
        if (dims === 3) {
          const z =
            coords instanceof Float32Array
              ? coords[p * dims + 2]
              : (coords as number[][])[p][2];

          if (z !== undefined) {
            bounds.minZ = Math.min(bounds.minZ, z);
            bounds.maxZ = Math.max(bounds.maxZ, z);
          }
        }
      }

      // Scale factor: 500x for normalized [-1,1] data, 1x for raw coordinates
      const coordScale = dataset.normalized === false ? 1 : 500;

      // Calculate center point (scaled)
      const center = new THREE.Vector3(
        ((bounds.minX + bounds.maxX) / 2) * coordScale,
        ((bounds.minY + bounds.maxY) / 2) * coordScale,
        dataset.spatial.dimensions === 3
          ? ((bounds.minZ + bounds.maxZ) / 2) * coordScale
          : 0,
      );

      // Calculate size of data
      const size = Math.max(
        (bounds.maxX - bounds.minX) * coordScale,
        (bounds.maxY - bounds.minY) * coordScale,
        dataset.spatial.dimensions === 3
          ? (bounds.maxZ - bounds.minZ) * coordScale
          : 0,
      );

      // Position camera at appropriate distance
      const distance = size * 1.5;
      const cameraPos = new THREE.Vector3(
        center.x,
        center.y,
        center.z + distance,
      );

      // Set near/far planes based on data scale to avoid clipping
      const maxExtent = Math.max(size, distance);
      const near = maxExtent * 0.001;
      const far = maxExtent * 10;

      // Initialize Three.js scene with options
      const { scene, camera, renderer, controls, animate, dispose } =
        initializeScene(containerRef.current, {
          is2D: viewMode === "2D",
          cameraPosition: cameraPos,
          lookAtPosition: center,
          near,
          far,
        });

      // Store camera, renderer, and controls refs for raycasting + scale bar
      cameraRef.current = camera;
      rendererRef.current = renderer;
      controlsRef.current = controls;

      // Create tooltip
      if (!tooltipRef.current) {
        tooltipRef.current = createTooltip();
      }

      // Throttled intersection checking (50ms)
      let lastCheckTime = 0;
      const throttleDelay = 10;

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
          const clusterValueStr = clusterRef.current
            ? getClusterValue(clusterRef.current, index)
            : "";

          // Only toggle celltype if it's not a numerical cluster
          if (!isNumericalClusterRef.current) {
            // Toggle the cluster in selectedCelltypes
            toggleCelltype(clusterValueStr);
            // Mode is now automatically updated by toggleCelltype
          }
        }
      };

      // Right-click handler: open split screen with linked SM dataset
      const handleContextMenu = (event: MouseEvent) => {
        // Only works on left panel (not inside a split right panel)
        if (panelId) return;
        // Need a hovered point
        if (hoveredPointRef.current === null) return;
        // Need a link config for this dataset
        const linkConfig = linkConfigRef.current;

        if (!linkConfig) return;

        event.preventDefault();

        const index = hoveredPointRef.current;

        // "__all__" linkColumn means every cell maps to the same SM dataset
        if (linkConfig.linkColumn === "__all__") {
          const smUrl = linkConfig.links["__all__"];
          if (!smUrl) return;
          enableSplit();
          setRightPanelS3(smUrl, "sm");
          toast.info("Opening SM dataset");
          return;
        }

        // Find the link column in dataset clusters (always reads the configured column,
        // regardless of which column is currently selected for visualization)
        const linkCluster = dataset!.clusters?.find(
          (c) => c.column === linkConfig.linkColumn,
        );

        if (!linkCluster) {
          console.warn(
            `Link column "${linkConfig.linkColumn}" not found in dataset clusters`,
          );

          return;
        }

        const setValue = getClusterValue(linkCluster, index);
        const smUrl = linkConfig.links[setValue];

        if (!smUrl) {
          toast.warning(`No single molecule data available for "${setValue}"`);

          return;
        }

        enableSplit();
        setRightPanelS3(smUrl, "sm");
        toast.info(`Opening SM dataset for ${setValue}`);
      };

      // Add event listeners
      renderer.domElement.addEventListener("mousemove", handleMouseMove);
      renderer.domElement.addEventListener("dblclick", handleDoubleClick);
      renderer.domElement.addEventListener("contextmenu", handleContextMenu);

      // Build positions Float32Array directly from dataset coordinates
      const spatialCoords = dataset.spatial.coordinates;
      const spatialDims = dataset.spatial.dimensions;
      const ptCount = dataset.getPointCount();
      const positions = new Float32Array(ptCount * 3);

      // Scale factor: 500x for normalized [-1,1] data, 1x for raw coordinates
      const cs = dataset.normalized === false ? 1 : 500;

      if (spatialCoords instanceof Float32Array) {
        // Flat Float32Array path (optimized — no object creation)
        if (spatialDims === 3 && cs === 1) {
          // Best case: 3D raw coords — direct copy, no scaling needed
          positions.set(spatialCoords);
        } else {
          for (let p = 0; p < ptCount; p++) {
            positions[p * 3] = spatialCoords[p * spatialDims] * cs;
            positions[p * 3 + 1] = spatialCoords[p * spatialDims + 1] * cs;
            positions[p * 3 + 2] =
              spatialDims === 3 ? spatialCoords[p * spatialDims + 2] * cs : 0;
          }
        }
      } else {
        // number[][] path
        const coords = spatialCoords as number[][];
        for (let p = 0; p < ptCount; p++) {
          positions[p * 3] = coords[p][0] * cs;
          positions[p * 3 + 1] = coords[p][1] * cs;
          positions[p * 3 + 2] =
            spatialDims === 3 ? (coords[p][2] ?? 0) * cs : 0;
        }
      }

      // Create point cloud mesh with custom shaders
      // The shader computes: gl_PointSize = size * dotSize * proj[1][1] / -mvPosition.z
      // We want ~3px dots at the initial zoom level.
      // Back-calculate: dotSize = targetPx * distance / (baseSize * proj[1][1])
      // where distance ≈ size*1.5, proj[1][1] ≈ 1.3 (75° FOV), baseSize = POINT_BASE_SIZE
      const targetPx = 0.1;
      const proj11 = 1.0 / Math.tan((75 * Math.PI) / 180 / 2); // ~1.3
      const baseDotSize =
        (targetPx * distance) / (VISUALIZATION_CONFIG.POINT_BASE_SIZE * proj11);
      baseDotSizeRef.current = baseDotSize;
      const pointCloud = createPointCloudFromBuffers(
        positions,
        ptCount,
        baseDotSize,
      );

      pointCloudRef.current = pointCloud; // Store reference
      sceneRef.current = scene; // Store scene reference
      scene.add(pointCloud);
      setPointCloudVersion((v) => v + 1); // Trigger visualization update

      // Start animation
      animate();

      // Cleanup on unmount
      return () => {
        // Remove event listeners
        renderer.domElement.removeEventListener("mousemove", handleMouseMove);
        renderer.domElement.removeEventListener("dblclick", handleDoubleClick);
        renderer.domElement.removeEventListener(
          "contextmenu",
          handleContextMenu,
        );

        // Hide tooltip
        hideTooltip();

        scene.remove(pointCloud);
        pointCloud.geometry.dispose();
        (pointCloud.material as any).dispose();
        pointCloudRef.current = null;
        sceneRef.current = null;
        cameraRef.current = null;
        rendererRef.current = null;
        controlsRef.current = null;
        dispose();
      };
    } else {
      // Initialize scene without dataset
      const { animate, dispose } = initializeScene(containerRef.current);

      // Start animation even without data
      animate();

      return () => {
        dispose();
      };
    }
  }, [dataset, viewMode]);

  // Effect 2: Update visualization based on mode array
  useEffect(() => {
    if (!pointCloudRef.current || !dataset) return;

    const updateVisualization = async () => {
      if (!pointCloudRef.current || !dataset) return;

      // Store cluster data for tooltip
      const selectedCluster = dataset.clusters?.find(
        (c) => c.column === selectedColumn,
      );

      if (selectedCluster) {
        clusterRef.current = selectedCluster;
        isNumericalClusterRef.current = selectedColumn
          ? getEffectiveColumnType(
              selectedColumn,
              dataset,
              columnTypeOverrides,
            ) === "numerical"
          : false;
        colorPaletteRef.current = selectedCluster.palette || colorPalette;
      }

      // Determine which visualization to use based on mode array
      const hasGeneMode = mode.includes("gene");
      const hasCelltypeMode = mode.includes("celltype");

      // Check if the selected column is numerical (respects overrides)
      const isNumerical = selectedColumn
        ? getEffectiveColumnType(
            selectedColumn,
            dataset,
            columnTypeOverrides,
          ) === "numerical"
        : false;

      let result = null;

      // Check if gene has changed (for auto-scaling)
      const geneChanged = previousGeneRef.current !== selectedGene;

      if (geneChanged) {
        previousGeneRef.current = selectedGene;
      }

      // Check if column or its effective type has changed (for auto-scaling numerical columns)
      const columnChanged = previousColumnRef.current !== selectedColumn;
      const typeChanged = previousColumnTypeRef.current !== isNumerical;

      if (columnChanged) {
        previousColumnRef.current = selectedColumn;
      }
      if (typeChanged) {
        previousColumnTypeRef.current = isNumerical;
      }

      const shouldAutoScale = columnChanged || typeChanged;

      if (
        hasGeneMode &&
        hasCelltypeMode &&
        selectedGene &&
        selectedCelltypes.size > 0
      ) {
        // Combined mode: gene expression on selected celltypes
        try {
          if (geneChanged) {
            geneToastIdRef.current = toast.loading(
              `Loading expression for "${selectedGene}"...`,
            );
          }

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
            geneScaleMin,
            geneScaleMax,
            // Only auto-set scale when gene changes, not when user manually adjusts
            geneChanged ? setGeneScaleMin : undefined,
            geneChanged ? setGeneScaleMax : undefined,
          );
        } finally {
          if (geneToastIdRef.current != null) {
            toast.dismiss(geneToastIdRef.current);
            geneToastIdRef.current = null;
          }
        }
      } else if (hasGeneMode && selectedGene) {
        // Gene mode only
        try {
          if (geneChanged) {
            geneToastIdRef.current = toast.loading(
              `Loading expression for "${selectedGene}"...`,
            );
          }

          // Fetch gene expression data for tooltip
          const expression = await dataset.getGeneExpression(selectedGene);

          geneExpressionRef.current = expression;

          result = await updateGeneVisualization(
            dataset,
            selectedGene,
            alphaScale,
            sizeScale,
            geneScaleMin,
            geneScaleMax,
            // Only auto-set scale when gene changes, not when user manually adjusts
            geneChanged ? setGeneScaleMin : undefined,
            geneChanged ? setGeneScaleMax : undefined,
          );
        } finally {
          if (geneToastIdRef.current != null) {
            toast.dismiss(geneToastIdRef.current);
            geneToastIdRef.current = null;
          }
        }
      } else if (hasCelltypeMode) {
        // Celltype mode only

        // Use appropriate visualization function based on column type
        result = isNumerical
          ? updateNumericalCelltypeVisualization(
              dataset,
              selectedColumn,
              alphaScale,
              sizeScale,
              numericalScaleMin,
              numericalScaleMax,
              // Only auto-set scale when column or type changes, not when user manually adjusts
              shouldAutoScale ? setNumericalScaleMin : undefined,
              shouldAutoScale ? setNumericalScaleMax : undefined,
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
        // Apply current sizeScale to dotSize uniform
        updateDotSize(
          pointCloudRef.current,
          baseDotSizeRef.current * sizeScale,
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
    geneScaleMin,
    geneScaleMax,
    numericalScaleMin,
    numericalScaleMax,
    mode,
    clusterVersion,
    columnTypeOverrides,
    pointCloudVersion,
  ]);

  // Effect 3: Update dotSize uniform when slider changes (instant, no per-point loop)
  useEffect(() => {
    if (pointCloudRef.current) {
      updateDotSize(pointCloudRef.current, baseDotSizeRef.current * sizeScale);
    }
  }, [sizeScale]);

  return (
    <>
      <div
        ref={containerRef}
        className="absolute inset-0 w-full h-full"
        style={{ margin: 0, padding: 0 }}
      />

      {/* Visualization legends panel (includes scale bar) */}
      <VisualizationLegends />

      {/* Spatial scale bar - only for non-normalized (raw coordinate) datasets */}
      {dataset && !dataset.normalized && (
        <SpatialScaleBar
          cameraRef={
            cameraRef as React.RefObject<THREE.PerspectiveCamera | null>
          }
          rendererRef={rendererRef}
          controlsRef={controlsRef}
        />
      )}
    </>
  );
}
