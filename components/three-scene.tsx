"use client";

import type { StandardizedDataset } from "@/lib/StandardizedDataset";
import type { SingleMoleculeDataset } from "@/lib/SingleMoleculeDataset";
import type { MoleculeShape } from "@/lib/stores/createSingleMoleculeVisualizationStore";

import { useEffect, useRef, useState } from "react";
import * as THREE from "three";
import { toast } from "react-toastify";
import gsap from "gsap";

import { getClusterValue } from "@/lib/StandardizedDataset";
import { useSingleMoleculeVisualizationStore } from "@/lib/stores/singleMoleculeVisualizationStore";
import { initializeScene } from "@/lib/webgl/scene-manager";
import {
  createPointCloudFromBuffers,
  createSmPointCloud,
  updatePointCloudAttributes,
  updatePointCloudColor,
  updatePointCloudShape,
  updateDotSize,
  SM_DOT_SIZE_FACTOR,
} from "@/lib/webgl/point-cloud";
import {
  updateGeneVisualization,
  updateCelltypeVisualization,
  updateNumericalCelltypeVisualization,
  updateCombinedVisualization,
  updateCoexpressionVisualization,
  type AdvancedVizSettings,
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
import { SingleMoleculeLegends } from "@/components/single-molecule-legends";
import { getEffectiveColumnType } from "@/lib/utils/column-type-utils";
import { ExportBoxOverlay } from "@/components/export-box-overlay";
import {
  markSceneReady,
  markRenderComplete,
  markGeneRendered,
} from "@/lib/utils/test-hooks";
import { SpatialScaleBar } from "@/components/spatial-scale-bar";
import { VISUALIZATION_CONFIG } from "@/lib/config/visualization.config";
import {
  applyTransformsToPositions,
  computeSampleCentroids,
} from "@/lib/utils/sample-transforms";
import { useTourPlayer } from "@/lib/hooks/useTourPlayer";
import { useTourStore } from "@/lib/stores/tourStore";

interface ThreeSceneProps {
  dataset?: StandardizedDataset | null;
}

export function ThreeScene({ dataset }: ThreeSceneProps) {
  const containerRef = useRef<HTMLDivElement>(null);
  const pointCloudRef = useRef<THREE.Points | null>(null);
  const [pointCloudVersion, setPointCloudVersion] = useState(0);
  const sceneRef = useRef<THREE.Scene | null>(null);
  const sceneGroupRef = useRef<THREE.Group | null>(null);
  // Untransformed scene-space positions (xyz interleaved). Used as the base
  // when applying per-sample transforms so edits compose from identity.
  const basePositionsRef = useRef<Float32Array | null>(null);
  const sceneScaleRef = useRef<number>(1);
  const geneToastIdRef = useRef<string | number | null>(null);

  // Raycaster and interaction state
  const raycasterRef = useRef<THREE.Raycaster>(new THREE.Raycaster());
  // World-space lookAt captured at scene init — restored by the camera reset.
  const initialLookAtRef = useRef<THREE.Vector3 | null>(null);
  // gsap.ticker callback for the continuous O-orbit. Set while orbiting.
  const orbitTickerRef = useRef<
    ((time: number, deltaTime: number) => void) | null
  >(null);
  // Mirror orbit settings into refs so the ticker reads current values
  // without having to restart the orbit on every slider drag.
  const orbitSpeedRef = useRef(1.0);
  const orbitYBobRef = useRef(0.3);
  // Distance-from-target filter state — pushed to shader uniforms each frame.
  // The radius is in absolute world units (microns for raw-coord datasets).
  const targetFilterEnabledRef = useRef(false);
  const targetFilterRadiusRef = useRef(5000);
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
  const pinnedClustersRef = useRef<
    Array<{
      column: string;
      cluster: any;
      palette: Record<string, string>;
      isNumerical: boolean;
    }>
  >([]);
  const baseDotSizeRef = useRef<number>(5);
  const cameraDistanceRef = useRef<number>(1);

  // SM overlay state (for __all__ mapping). Gene selection lives in the
  // global SM viz store so the picker + legends drive what's rendered.
  const smDatasetRef = useRef<SingleMoleculeDataset | null>(null);
  const smPointCloudsRef = useRef<Map<string, THREE.Points>>(new Map());
  const [smDataset, setSmDataset] = useState<SingleMoleculeDataset | null>(
    null,
  );
  const [smLoading, setSmLoading] = useState(false);
  // Last cmd/ctrl+click world coordinates; rendered as a small overlay chip
  // so the user can read them directly without opening devtools.
  const [lastClickCoords, setLastClickCoords] = useState<{
    x: number;
    y: number;
    z: number;
  } | null>(null);
  const smSelectedGenes = useSingleMoleculeVisualizationStore(
    (s) => s.selectedGenes,
  );
  const smSelectedGenesLegend = useSingleMoleculeVisualizationStore(
    (s) => s.selectedGenesLegend,
  );
  const smGlobalScale = useSingleMoleculeVisualizationStore(
    (s) => s.globalScale,
  );
  const smShowAssigned = useSingleMoleculeVisualizationStore(
    (s) => s.showAssigned,
  );
  const smShowUnassigned = useSingleMoleculeVisualizationStore(
    (s) => s.showUnassigned,
  );

  // Render-order helper: later in selectedGenesLegend → higher renderOrder.
  // Assigned cloud sits on top of its own unassigned. All values stay below
  // the SC point cloud (which uses the default renderOrder of 0).
  const smRenderOrderFor = (legendIndex: number, isAssigned: boolean): number =>
    -1000 + legendIndex * 2 + (isAssigned ? 1 : 0);
  // Mirror globalScale into a ref so Effect 5 can read the current value when
  // creating new point clouds without depending on globalScale (which would
  // re-create all clouds on every slider tick).
  const smGlobalScaleRef = useRef(smGlobalScale);

  useEffect(() => {
    smGlobalScaleRef.current = smGlobalScale;
  }, [smGlobalScale]);

  // Get visualization settings from store
  const {
    mode,
    selectedGene,
    selectedColumn,
    selectedCelltypes,
    geneEverywhere,
    hiddenCelltypes,
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
    cameraResetSignal,
    orbitSpeed,
    orbitYBob,
    targetFilterEnabled,
    targetFilterRadius,
    columnTypeOverrides,
    viewMode,
    targetPx,
    selectedSizeMultiplier,
    greyedOutSizeMultiplier,
    greyedOutAlpha,
    expressionAlphaMin,
    expressionAlphaMax,
    pointSizeMultiplierMin,
    pointSizeMultiplierMax,
    sceneRotation,
    flipX,
    flipY,
    pinnedTooltipColumns,
    colormap,
    transformColumn,
    sampleTransforms,
    transformVersion,
    coexpressEnabled,
    selectedGene2,
    coexpressSwapped,
    gene2ScaleMin,
    gene2ScaleMax,
    setGene2ScaleMin,
    setGene2ScaleMax,
    exportBoxEnabled,
    exportBoxWidthMm,
    exportBoxHeightMm,
    exportBoxCenterPx,
    setExportBoxCenterPx,
  } = usePanelVisualizationStore();

  // Mirror orbit settings into refs (read by gsap.ticker each frame so the
  // sliders update a running orbit live without restarting it).
  useEffect(() => {
    orbitSpeedRef.current = orbitSpeed;
  }, [orbitSpeed]);
  useEffect(() => {
    orbitYBobRef.current = orbitYBob;
  }, [orbitYBob]);
  useEffect(() => {
    targetFilterEnabledRef.current = targetFilterEnabled;
  }, [targetFilterEnabled]);
  useEffect(() => {
    targetFilterRadiusRef.current = targetFilterRadius;
  }, [targetFilterRadius]);

  // Tour player — drives camera + SM gene selection from the loaded tour JSON.
  // No-op until the user uploads a tour file and presses Play.
  const setTargetFilterEnabledFn = usePanelVisualizationStore(
    (s) => s.setTargetFilterEnabled,
  );
  const setTargetFilterRadiusFn = usePanelVisualizationStore(
    (s) => s.setTargetFilterRadius,
  );

  useTourPlayer({
    cameraRef,
    controlsRef,
    orbitTickerRef,
    isNormalized: dataset?.normalized !== false,
    targetFilterEnabledRef,
    targetFilterRadiusRef,
    setTargetFilterEnabled: setTargetFilterEnabledFn,
    setTargetFilterRadius: setTargetFilterRadiusFn,
    scPointCloudRef: pointCloudRef,
  });

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

        // "__all__" means every cell links to the same SM dataset — load overlay
        if (mappingConfig.linkColumn === "__all__") {
          const smUrl = mappingConfig.links["__all__"];

          if (!smUrl) return;

          // Load SM dataset for overlay
          setSmLoading(true);
          try {
            const { SingleMoleculeDataset: SMDataset } = await import(
              "@/lib/SingleMoleculeDataset"
            );
            const smDs = await SMDataset.fromCustomS3(smUrl);

            smDatasetRef.current = smDs;
            setSmDataset(smDs);

            // Push into the (global) SM dataset store so the SM gene picker
            // and legends see this dataset on the SC overlay page.
            const { useSingleMoleculeStore } = await import(
              "@/lib/stores/singleMoleculeStore"
            );

            useSingleMoleculeStore.getState().addDataset(smDs);
          } catch (error) {
            console.warn("Failed to load SM overlay dataset:", error);
          } finally {
            setSmLoading(false);
          }

          return;
        }

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
  const previousGene2Ref = useRef<string | null>(null);
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

    // Build pinned column rows
    const pinnedRows = pinnedClustersRef.current
      .map((p) => {
        const val = getClusterValue(p.cluster, index);

        if (p.isNumerical) {
          return `<div class="text-xs text-default-400">${p.column}: ${val}</div>`;
        }
        const color = p.palette[String(val)] || "#808080";

        return `<div class="flex items-center text-xs">
        <div style="width: 10px; height: 10px; border-radius: 50%; background-color: ${color}; margin-right: 5px;"></div>
        <span class="text-default-400">${p.column}: ${val}</span>
      </div>`;
      })
      .join("");

    let tooltipContent = "";

    if (currentMode.includes("gene") && currentGene) {
      // Gene mode
      if (isNumerical) {
        tooltipContent = `
          <div class="flex flex-col gap-1">
            <div>${currentColumn}: ${clusterValue}</div>
            ${pinnedRows}
            <div>${currentGene}: ${geneValue?.toFixed(2) ?? "N/A"}</div>
          </div>
        `;
      } else {
        const clusterColor =
          colorPaletteRef.current[String(clusterValue)] || "#808080";

        tooltipContent = `
          <div class="flex flex-col gap-1">
            <div class="flex items-center">
              <div style="width: 12px; height: 12px; border-radius: 50%; background-color: ${clusterColor}; margin-right: 6px;"></div>
              <span>${clusterValue}</span>
            </div>
            ${pinnedRows}
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
        tooltipContent = `
          <div class="flex flex-col gap-1">
            <div>${currentColumn}: ${clusterValue}</div>
            ${pinnedRows}
          </div>
        `;
      } else {
        const clusterColor =
          colorPaletteRef.current[String(clusterValue)] || "#808080";

        tooltipContent = `
          <div class="flex flex-col gap-1">
            <div class="flex items-center">
              <div style="width: 12px; height: 12px; border-radius: 50%; background-color: ${clusterColor}; margin-right: 6px;"></div>
              <span>${clusterValue}</span>
            </div>
            ${pinnedRows}
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

      // Stash the initial lookAt so the camera reset can snap back to it.
      initialLookAtRef.current = center.clone();

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

        // "__all__" linkColumn means every cell maps to the same SM dataset,
        // which is already rendered as an overlay on this scene — no need to
        // open a split panel on right-click.
        if (linkConfig.linkColumn === "__all__") return;

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

      // Cmd/Ctrl+click: pan the camera (and its orbit/pan pivot) so the
      // clicked SC cell lands at the screen center. Translation only —
      // camera↔target distance and angle are preserved. Works in both 2D
      // (OrbitControls, top-down) and 3D (TrackballControls) because both
      // expose the same .target/.position API. ~350ms gsap tween; a second
      // click during the tween overrides cleanly.
      const handleClickRecenter = (event: MouseEvent) => {
        if (!(event.metaKey || event.ctrlKey)) return;
        const idx = hoveredPointRef.current;
        const pc = pointCloudRef.current;
        const ctrls = controlsRef.current;
        const cam = cameraRef.current;

        if (idx == null || !pc || !ctrls || !cam) return;

        // Stop a running orbit before panning — otherwise the orbit ticker
        // keeps overwriting the camera position during the pan tween.
        if (orbitTickerRef.current) {
          gsap.ticker.remove(orbitTickerRef.current);
          orbitTickerRef.current = null;
        }

        const posAttr = pc.geometry.attributes
          .position as THREE.BufferAttribute;
        const local = new THREE.Vector3(
          posAttr.getX(idx),
          posAttr.getY(idx),
          posAttr.getZ(idx),
        );

        pc.updateMatrixWorld();
        const world = local.applyMatrix4(pc.matrixWorld);

        setLastClickCoords({ x: world.x, y: world.y, z: world.z });

        const targetStart = ctrls.target.clone();
        const positionStart = cam.position.clone();
        const delta = world.clone().sub(targetStart);
        const positionEnd = positionStart.clone().add(delta);
        const tweenObj = { t: 0 };

        gsap.to(tweenObj, {
          t: 1,
          duration: 0.35,
          ease: "power2.out",
          overwrite: true,
          onUpdate: () => {
            ctrls.target.lerpVectors(targetStart, world, tweenObj.t);
            cam.position.lerpVectors(positionStart, positionEnd, tweenObj.t);
            ctrls.update();
          },
        });
      };

      // Add event listeners
      renderer.domElement.addEventListener("mousemove", handleMouseMove);
      renderer.domElement.addEventListener("dblclick", handleDoubleClick);
      renderer.domElement.addEventListener("contextmenu", handleContextMenu);
      renderer.domElement.addEventListener("click", handleClickRecenter);

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

      // Stash base (untransformed) positions and the scene-space scale so the
      // per-sample transform effect can recompute live without re-deriving.
      basePositionsRef.current = new Float32Array(positions);
      sceneScaleRef.current = cs;

      // Create point cloud mesh with custom shaders
      // The shader computes: gl_PointSize = size * dotSize * proj[1][1] / -mvPosition.z
      // We want ~3px dots at the initial zoom level.
      // Back-calculate: dotSize = targetPx * distance / (baseSize * proj[1][1])
      // where distance ≈ size*1.5, proj[1][1] ≈ 1.3 (75° FOV), baseSize = POINT_BASE_SIZE
      const proj11 = 1.0 / Math.tan((75 * Math.PI) / 180 / 2); // ~1.3

      cameraDistanceRef.current = distance;
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

      // Wrap point cloud in a group for scene transforms (rotation, flip)
      const group = new THREE.Group();

      group.add(pointCloud);
      scene.add(group);
      sceneGroupRef.current = group;
      setPointCloudVersion((v) => v + 1); // Trigger visualization update

      // Start animation
      animate();
      markSceneReady();

      // Cleanup on unmount
      return () => {
        // Remove event listeners
        renderer.domElement.removeEventListener("mousemove", handleMouseMove);
        renderer.domElement.removeEventListener("dblclick", handleDoubleClick);
        renderer.domElement.removeEventListener(
          "contextmenu",
          handleContextMenu,
        );
        renderer.domElement.removeEventListener("click", handleClickRecenter);

        // Hide tooltip
        hideTooltip();

        if (sceneGroupRef.current) {
          scene.remove(sceneGroupRef.current);
          sceneGroupRef.current = null;
        }
        pointCloud.geometry.dispose();
        (pointCloud.material as any).dispose();
        pointCloudRef.current = null;
        sceneRef.current = null;
        basePositionsRef.current = null;
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

  // Camera reset: snap controls.target back to the initial lookAt captured at
  // scene init. Triggered by cameraResetSignal (Reset View button) or R key.
  useEffect(() => {
    if (cameraResetSignal === 0) return;
    const ctrls = controlsRef.current;
    const initial = initialLookAtRef.current;

    if (!ctrls || !initial) return;
    ctrls.target.copy(initial);
    ctrls.update();
  }, [cameraResetSignal]);

  // R key resets the camera (skipped while typing in an input).
  useEffect(() => {
    const onKey = (e: KeyboardEvent) => {
      if (e.key !== "r" && e.key !== "R") return;
      if (e.metaKey || e.ctrlKey || e.altKey) return;
      const el = document.activeElement as HTMLElement | null;

      if (
        el &&
        (el.tagName === "INPUT" ||
          el.tagName === "TEXTAREA" ||
          el.isContentEditable)
      ) {
        return;
      }
      const ctrls = controlsRef.current;
      const initial = initialLookAtRef.current;

      if (!ctrls || !initial) return;
      ctrls.target.copy(initial);
      ctrls.update();
    };

    window.addEventListener("keydown", onKey);

    return () => window.removeEventListener("keydown", onKey);
  }, []);

  // O key toggles a continuous orbit around controls.target in 3D mode.
  // The camera rotates around world Y at ~4.5s per revolution and also bobs
  // up/down (Y offset) on a slower phase, so the view sweeps the cell from
  // multiple angles rather than tracing a flat ring. Pressing O again stops.
  useEffect(() => {
    const stopOrbit = () => {
      if (orbitTickerRef.current) {
        gsap.ticker.remove(orbitTickerRef.current);
        orbitTickerRef.current = null;
      }
    };

    const onKey = (e: KeyboardEvent) => {
      if (e.key !== "o" && e.key !== "O") return;
      if (e.metaKey || e.ctrlKey || e.altKey) return;
      const el = document.activeElement as HTMLElement | null;

      if (
        el &&
        (el.tagName === "INPUT" ||
          el.tagName === "TEXTAREA" ||
          el.isContentEditable)
      ) {
        return;
      }
      if (viewMode !== "3D") return;

      // Toggle off if already orbiting.
      if (orbitTickerRef.current) {
        stopOrbit();

        return;
      }

      const ctrls = controlsRef.current;
      const cam = cameraRef.current;

      if (!ctrls || !cam) return;

      const target = ctrls.target.clone();
      const startRadius = cam.position.clone().sub(target);
      const radiusXZ = Math.sqrt(
        startRadius.x * startRadius.x + startRadius.z * startRadius.z,
      );
      const baseY = startRadius.y;
      // Base period at speed 1.0 is 4.5s/rev. Live-updates via the ref.
      const BASE_PERIOD_MS = 4500;
      let angle = 0;
      // The Y bob runs at half the azimuth frequency so XZ and Y waves drift
      // out of sync — gives the "shifting figure-8" look.
      let elev = 0;

      const ticker = (_time: number, deltaTime: number) => {
        const speed = Math.max(0.01, orbitSpeedRef.current);

        angle += (deltaTime / (BASE_PERIOD_MS / speed)) * Math.PI * 2;
        elev += (deltaTime / ((BASE_PERIOD_MS * 2) / speed)) * Math.PI * 2;
        const sin = Math.sin(angle);
        const cos = Math.cos(angle);
        const yOsc = Math.sin(elev) * radiusXZ * orbitYBobRef.current;

        cam.position.set(
          target.x + startRadius.x * cos + startRadius.z * sin,
          target.y + baseY + yOsc,
          target.z + -startRadius.x * sin + startRadius.z * cos,
        );
        ctrls.update();
      };

      orbitTickerRef.current = ticker;
      gsap.ticker.add(ticker);
    };

    window.addEventListener("keydown", onKey);

    return () => {
      window.removeEventListener("keydown", onKey);
      stopOrbit();
    };
  }, [viewMode]);

  // T key toggles tour playback. Works even while the UI is hidden (H), so a
  // tour can be started/stopped for a clean screen recording. Requires a
  // loaded tour config and 3D view (same gate as the Play tour button).
  useEffect(() => {
    const onKey = (e: KeyboardEvent) => {
      if (e.key !== "t" && e.key !== "T") return;
      if (e.metaKey || e.ctrlKey || e.altKey) return;
      const el = document.activeElement as HTMLElement | null;

      if (
        el &&
        (el.tagName === "INPUT" ||
          el.tagName === "TEXTAREA" ||
          el.isContentEditable)
      ) {
        return;
      }
      const tour = useTourStore.getState();

      if (tour.isPlaying) {
        tour.stop();

        return;
      }
      if (!tour.config || viewMode !== "3D") return;
      tour.start();
    };

    window.addEventListener("keydown", onKey);

    return () => window.removeEventListener("keydown", onKey);
  }, [viewMode]);

  // Effect 1b: Apply per-sample transforms to the live position buffer.
  // Reads base positions stashed during scene init and writes the result
  // into geometry.attributes.position. Marks needsUpdate so the GPU re-uploads.
  const centroidsCacheRef = useRef<{
    column: string | null;
    centroids: Map<string, [number, number]> | null;
  }>({ column: null, centroids: null });

  // Invalidate centroid cache when dataset changes.
  useEffect(() => {
    centroidsCacheRef.current = { column: null, centroids: null };
  }, [dataset]);

  useEffect(() => {
    const base = basePositionsRef.current;
    const pc = pointCloudRef.current;

    if (!base || !pc || !dataset) return;

    const positionAttr = pc.geometry.attributes
      .position as THREE.BufferAttribute;
    const out = positionAttr.array as Float32Array;

    if (!transformColumn || sampleTransforms.size === 0) {
      out.set(base);
      positionAttr.needsUpdate = true;

      return;
    }

    const cluster =
      dataset.clusters?.find((c) => c.column === transformColumn) ?? null;

    if (!cluster) {
      // Column not loaded yet — fall back to base positions.
      out.set(base);
      positionAttr.needsUpdate = true;

      return;
    }

    // Centroids depend only on (spatial, column); cache them across edits.
    let centroids = centroidsCacheRef.current.centroids;

    if (centroidsCacheRef.current.column !== transformColumn) {
      centroids = computeSampleCentroids(dataset.spatial, cluster);
      centroidsCacheRef.current = { column: transformColumn, centroids };
    }

    applyTransformsToPositions(
      base,
      out,
      cluster,
      centroids,
      sampleTransforms,
      sceneScaleRef.current,
    );
    positionAttr.needsUpdate = true;
  }, [
    dataset,
    pointCloudVersion,
    transformColumn,
    transformVersion,
    clusterVersion,
  ]);

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

      // Populate pinned columns for tooltip
      pinnedClustersRef.current = [];
      for (const col of pinnedTooltipColumns) {
        if (col === selectedColumn) continue; // active column already shown
        const cluster = dataset.clusters?.find((c) => c.column === col);

        if (!cluster) continue;
        const isNum =
          getEffectiveColumnType(col, dataset, columnTypeOverrides) ===
          "numerical";

        pinnedClustersRef.current.push({
          column: col,
          cluster,
          palette: cluster.palette || {},
          isNumerical: isNum,
        });
      }

      // Build advanced settings from store
      const adv: AdvancedVizSettings = {
        selectedSizeMultiplier,
        greyedOutSizeMultiplier,
        greyedOutAlpha,
        expressionAlphaMin,
        expressionAlphaMax,
        pointSizeMultiplierMin,
        pointSizeMultiplierMax,
      };

      // Determine which visualization to use based on mode array
      const hasGeneMode = mode.includes("gene");
      const hasCelltypeMode = mode.includes("celltype");

      // The 3D scene renders only the visible subset of the selection
      // (selectedCelltypes minus hiddenCelltypes, set via the per-badge eye
      // toggles). selectedCelltypes itself is untouched, so DEG / plots keep
      // the full selection.
      const effectiveCelltypes =
        hiddenCelltypes.size > 0
          ? new Set(
              [...selectedCelltypes].filter((ct) => !hiddenCelltypes.has(ct)),
            )
          : selectedCelltypes;

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

      // Check if 2nd gene has changed (for coexpression auto-scaling).
      const gene2Changed = previousGene2Ref.current !== selectedGene2;

      if (gene2Changed) {
        previousGene2Ref.current = selectedGene2;
      }

      const coexpressActive =
        coexpressEnabled && !!selectedGene && !!selectedGene2;

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

      if (coexpressActive) {
        // Two-gene coexpression: green + magenta additive blend.
        try {
          if (geneChanged || gene2Changed) {
            geneToastIdRef.current = toast.loading(
              `Loading "${selectedGene}" + "${selectedGene2}"…`,
            );
          }

          // Cache primary gene expression for tooltips (use gene 1).
          geneExpressionRef.current = await dataset.getGeneExpression(
            selectedGene as string,
          );

          const greyOutNonSelected =
            hasCelltypeMode && selectedCelltypes.size > 0 && !geneEverywhere;

          result = await updateCoexpressionVisualization(
            dataset,
            selectedGene as string,
            selectedGene2 as string,
            geneScaleMin,
            geneScaleMax,
            geneChanged ? setGeneScaleMin : undefined,
            geneChanged ? setGeneScaleMax : undefined,
            gene2ScaleMin,
            gene2ScaleMax,
            gene2Changed ? setGene2ScaleMin : undefined,
            gene2Changed ? setGene2ScaleMax : undefined,
            coexpressSwapped,
            alphaScale,
            sizeScale,
            selectedColumn,
            effectiveCelltypes,
            greyOutNonSelected,
            adv,
          );
        } finally {
          if (geneToastIdRef.current != null) {
            toast.dismiss(geneToastIdRef.current);
            geneToastIdRef.current = null;
          }
        }
      } else if (
        hasGeneMode &&
        hasCelltypeMode &&
        selectedGene &&
        selectedCelltypes.size > 0 &&
        !geneEverywhere
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
            effectiveCelltypes,
            alphaScale,
            sizeScale,
            geneScaleMin,
            geneScaleMax,
            geneChanged ? setGeneScaleMin : undefined,
            geneChanged ? setGeneScaleMax : undefined,
            adv,
            colormap,
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
            geneChanged ? setGeneScaleMin : undefined,
            geneChanged ? setGeneScaleMax : undefined,
            adv,
            colormap,
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
              shouldAutoScale ? setNumericalScaleMin : undefined,
              shouldAutoScale ? setNumericalScaleMax : undefined,
              adv,
              colormap,
            )
          : updateCelltypeVisualization(
              dataset,
              selectedColumn,
              effectiveCelltypes,
              colorPalette,
              alphaScale,
              sizeScale,
              adv,
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

      // E2E timing hook: mark this render complete (covers initial render and
      // every gene/celltype change). Stats are read straight off the dataset.
      markRenderComplete({
        pointCount: dataset.getPointCount(),
        geneCount: dataset.genes.length,
        dimensions: dataset.spatial.dimensions,
        dataType: "single_cell",
      });
      if (mode.includes("gene") && selectedGene) {
        markGeneRendered(selectedGene, dataset.getPointCount());
      }
    };

    updateVisualization();
  }, [
    dataset,
    selectedGene,
    selectedColumn,
    selectedCelltypes,
    geneEverywhere,
    hiddenCelltypes,
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
    selectedSizeMultiplier,
    greyedOutSizeMultiplier,
    greyedOutAlpha,
    expressionAlphaMin,
    expressionAlphaMax,
    pointSizeMultiplierMin,
    pointSizeMultiplierMax,
    pinnedTooltipColumns,
    colormap,
    coexpressEnabled,
    selectedGene2,
    coexpressSwapped,
    gene2ScaleMin,
    gene2ScaleMax,
  ]);

  // Effect 3: Update dotSize uniform when slider or targetPx changes (instant)
  useEffect(() => {
    if (pointCloudRef.current) {
      const proj11 = 1.0 / Math.tan((75 * Math.PI) / 180 / 2);
      const newBaseDotSize =
        (targetPx * cameraDistanceRef.current) /
        (VISUALIZATION_CONFIG.POINT_BASE_SIZE * proj11);

      baseDotSizeRef.current = newBaseDotSize;
      updateDotSize(pointCloudRef.current, newBaseDotSize * sizeScale);
    }
  }, [sizeScale, targetPx]);

  // Effect 4: Apply scene transforms (rotation, flip) — instant via group matrix
  useEffect(() => {
    const group = sceneGroupRef.current;

    if (!group) return;

    group.rotation.z = (sceneRotation * Math.PI) / 180;
    group.scale.x = flipX ? -1 : 1;
    group.scale.y = flipY ? -1 : 1;
  }, [sceneRotation, flipX, flipY]);

  // Effect 5: Manage SM overlay point clouds.
  // Keys are composite: `${gene}:a` for assigned and `${gene}:u` for unassigned.
  // Visibility per gene = (global flag) AND (per-gene flag) AND (dataset has it).
  useEffect(() => {
    const parent = sceneGroupRef.current ?? sceneRef.current;
    const smDs = smDatasetRef.current;

    if (!parent || !smDs || smSelectedGenes.size === 0) {
      smPointCloudsRef.current.forEach((pc) => {
        pc.parent?.remove(pc);
        pc.geometry.dispose();
        (pc.material as THREE.ShaderMaterial).dispose();
      });
      smPointCloudsRef.current.clear();

      return;
    }

    const currentClouds = smPointCloudsRef.current;

    // Drop orphaned clouds whose parent group was disposed (e.g. on a 2D↔3D
    // view-mode switch, which rebuilds the scene). They'll get recreated
    // below under the fresh sceneGroupRef.
    for (const [key, pc] of currentClouds) {
      if (pc.parent !== parent) {
        pc.geometry.dispose();
        (pc.material as THREE.ShaderMaterial).dispose();
        currentClouds.delete(key);
      }
    }

    // Remove clouds for deselected genes (both keys).
    for (const [key, pc] of currentClouds) {
      const gene = key.slice(0, -2);

      if (!smSelectedGenes.has(gene)) {
        pc.parent?.remove(pc);
        pc.geometry.dispose();
        (pc.material as THREE.ShaderMaterial).dispose();
        currentClouds.delete(key);
      }
    }

    // Snapshot legend order so each cloud's renderOrder is set deterministically
    // from the user's legend (later-inserted gene renders on top).
    const legendOrder = new Map<string, number>();

    {
      let i = 0;

      for (const gene of smSelectedGenesLegend) legendOrder.set(gene, i++);
    }

    const makeCloud = (
      coords: Float32Array,
      colorHex: string,
      localScale: number,
      shape: MoleculeShape,
      gene: string,
      isAssigned: boolean,
    ): THREE.Points => {
      const userSize =
        VISUALIZATION_CONFIG.SINGLE_MOLECULE_POINT_BASE_SIZE *
        smGlobalScaleRef.current *
        localScale;
      const pc = createSmPointCloud(
        coords,
        coords.length / 3,
        colorHex,
        userSize * SM_DOT_SIZE_FACTOR,
        shape,
      );

      pc.renderOrder = smRenderOrderFor(legendOrder.get(gene) ?? 0, isAssigned);

      return pc;
    };

    const removeKey = (key: string) => {
      const pc = currentClouds.get(key);

      if (!pc) return;
      pc.parent?.remove(pc);
      pc.geometry.dispose();
      (pc.material as THREE.ShaderMaterial).dispose();
      currentClouds.delete(key);
    };

    // Per-invocation cancellation. If deps change while async loads are in
    // flight (typical during auto-select when addGene fires three times back
    // to back), the previous invocation's awaits would otherwise resolve and
    // attach duplicate clouds to the scene. Cleanup flips this so post-await
    // code exits without doing work.
    let cancelled = false;

    const work = async () => {
      for (const [gene, viz] of smSelectedGenes) {
        if (cancelled) return;
        const aKey = `${gene}:a`;
        const uKey = `${gene}:u`;
        const wantAssigned = smShowAssigned && viz.showAssigned;
        const wantUnassigned =
          smDs.hasUnassigned && smShowUnassigned && viz.showUnassigned;
        const toastId = `loading-sm-${gene}`;
        const needsAssigned = wantAssigned && !currentClouds.has(aKey);
        const needsUnassigned = wantUnassigned && !currentClouds.has(uKey);

        if (needsAssigned || needsUnassigned) {
          toast.loading(`Loading ${gene}...`, {
            toastId,
            position: "bottom-left",
            autoClose: false,
          });
        }

        try {
          // Assigned
          if (needsAssigned) {
            const coords = await smDs.getCoordinatesByGene(gene);

            if (cancelled) return;
            if (coords && coords.length > 0 && !currentClouds.has(aKey)) {
              const pc = makeCloud(
                coords,
                viz.color,
                viz.localScale,
                viz.assignedShape,
                gene,
                true,
              );

              parent.add(pc);
              currentClouds.set(aKey, pc);
            }
          } else if (!wantAssigned && currentClouds.has(aKey)) {
            removeKey(aKey);
          }

          if (cancelled) return;

          // Unassigned
          if (needsUnassigned) {
            const uCoords = await smDs.getUnassignedCoordinatesByGene(gene);

            if (cancelled) return;
            if (uCoords && uCoords.length > 0 && !currentClouds.has(uKey)) {
              const pc = makeCloud(
                uCoords,
                viz.unassignedColor,
                viz.unassignedLocalScale,
                viz.unassignedShape,
                gene,
                false,
              );

              parent.add(pc);
              currentClouds.set(uKey, pc);
            }
          } else if (!wantUnassigned && currentClouds.has(uKey)) {
            removeKey(uKey);
          }
        } catch (error) {
          console.warn(`Failed to load SM gene ${gene}:`, error);
        } finally {
          if (needsAssigned || needsUnassigned) {
            toast.dismiss(toastId);
          }
        }
      }
    };

    work();

    return () => {
      cancelled = true;
    };
  }, [
    smDataset,
    smSelectedGenes,
    smShowAssigned,
    smShowUnassigned,
    pointCloudVersion,
  ]);

  // Effect 6: live-update SM material sizes when globalScale or per-gene
  // localScale changes — avoids tearing down/rebuilding the point clouds.
  useEffect(() => {
    for (const [key, pc] of smPointCloudsRef.current) {
      const gene = key.slice(0, -2);
      const isAssigned = key.endsWith(":a");
      const viz = smSelectedGenes.get(gene);

      if (!viz) continue;
      const localScale = isAssigned ? viz.localScale : viz.unassignedLocalScale;
      const userSize =
        VISUALIZATION_CONFIG.SINGLE_MOLECULE_POINT_BASE_SIZE *
        smGlobalScale *
        localScale;

      updateDotSize(pc, userSize * SM_DOT_SIZE_FACTOR);
    }
  }, [smGlobalScale, smSelectedGenes]);

  // Effect 7: live-update SM cloud colors + shapes from the SM viz store.
  // Runs whenever the selectedGenes Map changes (color/shape edits create a
  // new Map ref in the store).
  useEffect(() => {
    for (const [key, pc] of smPointCloudsRef.current) {
      const gene = key.slice(0, -2);
      const isAssigned = key.endsWith(":a");
      const viz = smSelectedGenes.get(gene);

      if (!viz) continue;
      const colorHex = isAssigned ? viz.color : viz.unassignedColor;
      const shape = isAssigned ? viz.assignedShape : viz.unassignedShape;

      updatePointCloudColor(pc, colorHex);
      updatePointCloudShape(pc, shape);
    }
  }, [smSelectedGenes]);

  // Distance-filter uniforms: pushed every frame because controls.target
  // mutates outside React (cmd+click pan, mouse drag, R reset, orbit). Cost
  // is a handful of Vector3 copies per frame — negligible.
  useEffect(() => {
    const pushFilterUniforms = () => {
      const ctrls = controlsRef.current;
      const pc = pointCloudRef.current;

      if (!ctrls) return;
      const enabled = targetFilterEnabledRef.current ? 1.0 : 0.0;
      const radius = targetFilterRadiusRef.current;
      const feather = Math.max(1.0, radius * 0.1);
      const center = ctrls.target;

      if (pc) {
        const mat = pc.material as THREE.ShaderMaterial;

        if (mat.uniforms.uTargetCenter) {
          (mat.uniforms.uTargetCenter.value as THREE.Vector3).copy(center);
          mat.uniforms.uTargetRadius.value = radius;
          mat.uniforms.uTargetFeather.value = feather;
          mat.uniforms.uTargetFilterEnabled.value = enabled;
        }
      }

      for (const [, smPc] of smPointCloudsRef.current) {
        const mat = smPc.material as THREE.ShaderMaterial;

        if (!mat.uniforms?.uTargetCenter) continue;
        (mat.uniforms.uTargetCenter.value as THREE.Vector3).copy(center);
        mat.uniforms.uTargetRadius.value = radius;
        mat.uniforms.uTargetFeather.value = feather;
        mat.uniforms.uTargetFilterEnabled.value = enabled;
      }
    };

    gsap.ticker.add(pushFilterUniforms);

    return () => gsap.ticker.remove(pushFilterUniforms);
  }, []);

  // Effect 8: live-update SM cloud render order from selectedGenesLegend.
  // Later-inserted genes get higher renderOrder → render on top. Deselect +
  // reactivate from the gene picker pushes the gene to the end of the Set,
  // so it comes to the front.
  useEffect(() => {
    const legendOrder = new Map<string, number>();
    let i = 0;

    for (const gene of smSelectedGenesLegend) legendOrder.set(gene, i++);

    for (const [key, pc] of smPointCloudsRef.current) {
      const gene = key.slice(0, -2);
      const isAssigned = key.endsWith(":a");
      const idx = legendOrder.get(gene) ?? 0;

      pc.renderOrder = smRenderOrderFor(idx, isAssigned);
    }
  }, [smSelectedGenesLegend]);

  return (
    <>
      <div
        ref={containerRef}
        className="absolute inset-0 w-full h-full"
        data-testid="sc-scene-canvas"
        style={{ margin: 0, padding: 0 }}
      />

      {/* SM loading indicator */}
      {smLoading && (
        <div className="absolute top-28 right-6 z-50 flex items-center gap-2 px-3 py-2 rounded-lg bg-black/60 backdrop-blur-sm text-white text-sm">
          <div className="w-4 h-4 border-2 border-white/30 border-t-white rounded-full animate-spin" />
          Loading molecules...
        </div>
      )}

      {/* Last cmd/ctrl+click coordinates. */}
      {lastClickCoords && (
        <div className="absolute bottom-4 left-4 z-50 px-3 py-1.5 rounded-lg bg-black/60 backdrop-blur-sm text-white text-xs font-mono">
          ({lastClickCoords.x.toFixed(2)}, {lastClickCoords.y.toFixed(2)},{" "}
          {lastClickCoords.z.toFixed(2)})
        </div>
      )}

      {/* SC + SM legends share a single absolute container so the SM gene
          pills stack below the SC celltype/gene legend instead of overlapping
          it at top-right. */}
      <div className="absolute right-6 top-24 z-10 flex flex-col items-end gap-4 max-w-xs">
        <VisualizationLegends embedded />
        <SingleMoleculeLegends embedded />
      </div>

      {/* Spatial scale bar — only for datasets with reliable raw coordinates */}
      {dataset && !dataset.metadata?.wasNormalized && (
        <SpatialScaleBar
          cameraRef={
            cameraRef as React.RefObject<THREE.PerspectiveCamera | null>
          }
          controlsRef={controlsRef}
          rendererRef={rendererRef}
        />
      )}

      {/* Export box overlay — gated by viewMode + raw-coords (same as scale bar). */}
      <ExportBoxOverlay
        cameraRef={cameraRef as React.RefObject<THREE.PerspectiveCamera | null>}
        canShow={
          viewMode === "2D" && !!dataset && !dataset.metadata?.wasNormalized
        }
        controlsRef={controlsRef}
        rendererRef={rendererRef}
        setCenterPx={setExportBoxCenterPx}
        state={{
          enabled: exportBoxEnabled,
          widthMm: exportBoxWidthMm,
          heightMm: exportBoxHeightMm,
          centerPx: exportBoxCenterPx,
        }}
      />
    </>
  );
}
