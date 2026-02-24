"use client";

import type { StandardizedDataset } from "@/lib/StandardizedDataset";
import type { PointData } from "@/lib/webgl/types";

import { motion, AnimatePresence } from "framer-motion";
import { useState, useEffect, useRef, useCallback } from "react";
import { Tooltip, Select, SelectItem } from "@heroui/react";
import { toast } from "react-toastify";
import { TbChartDots3 } from "react-icons/tb";
import { IoClose } from "react-icons/io5";
import * as THREE from "three";
import { OrbitControls } from "three/examples/jsm/controls/OrbitControls.js";

import {
  usePanelDatasetStore,
  usePanelVisualizationStore,
} from "@/lib/hooks/usePanelStores";
import {
  createPointCloud,
  updatePointCloudAttributes,
} from "@/lib/webgl/point-cloud";
import { normalizeCoordinates } from "@/lib/utils/coordinates";
import {
  updateGeneVisualization,
  updateCelltypeVisualization,
  updateNumericalCelltypeVisualization,
  updateCombinedVisualization,
} from "@/lib/webgl/visualization-utils";
import { VISUALIZATION_CONFIG } from "@/lib/config/visualization.config";

export default function UMAPPanel() {
  const [isOpen, setIsOpen] = useState(false);
  const [selectedEmbedding, setSelectedEmbedding] = useState<string>("");
  const [availableEmbeddings, setAvailableEmbeddings] = useState<string[]>([]);
  const [sceneReady, setSceneReady] = useState(false);
  const [pointCloudVersion, setPointCloudVersion] = useState(0);

  const containerRef = useRef<HTMLDivElement>(null);
  const sceneRef = useRef<THREE.Scene | null>(null);
  const cameraRef = useRef<THREE.OrthographicCamera | null>(null);
  const rendererRef = useRef<THREE.WebGLRenderer | null>(null);
  const controlsRef = useRef<OrbitControls | null>(null);
  const pointCloudRef = useRef<THREE.Points | null>(null);
  const animationIdRef = useRef<number | null>(null);

  // Raycaster and interaction state
  const raycasterRef = useRef<THREE.Raycaster>(new THREE.Raycaster());
  const mouseRef = useRef<THREE.Vector2>(new THREE.Vector2());
  const tooltipRef = useRef<HTMLDivElement | null>(null);
  const hoveredPointRef = useRef<number | null>(null);

  const dataset = usePanelDatasetStore((state) => state.getCurrentDataset()) as
    | StandardizedDataset
    | null
    | undefined;

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
    numericalScaleMin,
    numericalScaleMax,
    toggleCelltype,
  } = usePanelVisualizationStore();

  // Store current visualization data for tooltips
  const geneExpressionRef = useRef<number[] | null>(null);
  const colorPaletteRef = useRef<Record<string, string>>({});
  const clusterValuesRef = useRef<(string | number)[]>([]);
  const isNumericalClusterRef = useRef<boolean>(false);

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

  // Check if embeddings are available and auto-select
  useEffect(() => {
    if (!dataset) {
      setAvailableEmbeddings([]);
      setSelectedEmbedding("");

      return;
    }

    // Prefer allEmbeddingNames (lazy-load aware) over loaded keys
    const embeddingKeys =
      dataset.allEmbeddingNames && dataset.allEmbeddingNames.length > 0
        ? dataset.allEmbeddingNames
        : Object.keys(dataset.embeddings ?? {});

    setAvailableEmbeddings(embeddingKeys);

    // Auto-select: umap > pca > first available (panel stays closed by default)
    if (embeddingKeys.length > 0) {
      if (embeddingKeys.includes("umap")) {
        setSelectedEmbedding("umap");
      } else if (embeddingKeys.includes("pca")) {
        setSelectedEmbedding("pca");
      } else {
        setSelectedEmbedding(embeddingKeys[0]);
      }
    }
  }, [dataset]);

  // Helper function: Create tooltip element
  const createTooltip = (): HTMLDivElement => {
    const tooltip = document.createElement("div");

    tooltip.className =
      "absolute bg-black/80 text-white px-2.5 py-1.5 rounded text-sm font-sans pointer-events-none hidden z-[10000] shadow-lg min-w-[80px]";
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

    return `rgb(${r}, ${g}, ${b})`;
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
        // Categorical cluster: show with color circle from palette
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

  // Initialize Three.js scene
  useEffect(() => {
    if (!containerRef.current || !isOpen) return;

    // If scene already exists, just mark it ready
    if (sceneRef.current) {
      setSceneReady(true);

      return;
    }

    // Create scene
    const scene = new THREE.Scene();

    scene.background = new THREE.Color(0x000000);
    sceneRef.current = scene;

    // Create orthographic camera for 2D view
    const aspect = 1; // Square container
    const frustumSize = 20; // Larger view area to see more points
    const camera = new THREE.OrthographicCamera(
      (-frustumSize * aspect) / 2,
      (frustumSize * aspect) / 2,
      frustumSize / 2,
      -frustumSize / 2,
      0.1,
      1000,
    );

    camera.position.set(0, 0, 20);
    camera.lookAt(0, 0, 0);
    cameraRef.current = camera;

    // Create renderer
    const renderer = new THREE.WebGLRenderer({ antialias: true });
    const container = containerRef.current;

    renderer.setSize(container.clientWidth, container.clientHeight);
    renderer.setPixelRatio(window.devicePixelRatio);
    container.appendChild(renderer.domElement);
    rendererRef.current = renderer;

    // Create controls (2D - no rotation)
    const controls = new OrbitControls(camera, renderer.domElement);

    controls.enableRotate = false;
    controls.enableDamping = true;
    controls.dampingFactor = 0.05;
    controls.screenSpacePanning = true;
    controls.mouseButtons = {
      LEFT: THREE.MOUSE.PAN,
      MIDDLE: THREE.MOUSE.DOLLY,
      RIGHT: THREE.MOUSE.PAN,
    };
    controlsRef.current = controls;

    // Create tooltip
    tooltipRef.current = createTooltip();

    // Raycaster setup
    raycasterRef.current.params.Points = { threshold: 0.5 };

    // Animation loop
    const animate = () => {
      animationIdRef.current = requestAnimationFrame(animate);
      controls.update();
      renderer.render(scene, camera);
    };

    animate();

    // Mark scene as ready
    setSceneReady(true);

    // Handle resize
    const handleResize = () => {
      if (!containerRef.current || !cameraRef.current || !rendererRef.current)
        return;

      const width = containerRef.current.clientWidth;
      const height = containerRef.current.clientHeight;

      rendererRef.current.setSize(width, height);
    };

    window.addEventListener("resize", handleResize);

    // Cleanup
    return () => {
      if (animationIdRef.current) {
        cancelAnimationFrame(animationIdRef.current);
      }
      window.removeEventListener("resize", handleResize);
      if (tooltipRef.current) {
        document.body.removeChild(tooltipRef.current);
        tooltipRef.current = null;
      }
      controls.dispose();
      renderer.dispose();
      if (pointCloudRef.current) {
        scene.remove(pointCloudRef.current);
        pointCloudRef.current.geometry.dispose();
        if (pointCloudRef.current.material instanceof THREE.Material) {
          pointCloudRef.current.material.dispose();
        }
      }
      if (containerRef.current && renderer.domElement) {
        containerRef.current.removeChild(renderer.domElement);
      }
      sceneRef.current = null;
      cameraRef.current = null;
      rendererRef.current = null;
      controlsRef.current = null;
      pointCloudRef.current = null;
      setSceneReady(false);
    };
  }, [isOpen]);

  // Helper: create point cloud from embedding data array
  const buildPointCloud = useCallback(
    (embeddingData: number[][], scene: THREE.Scene) => {
      // Remove existing point cloud
      if (pointCloudRef.current) {
        scene.remove(pointCloudRef.current);
        pointCloudRef.current.geometry.dispose();
        if (pointCloudRef.current.material instanceof THREE.Material) {
          pointCloudRef.current.material.dispose();
        }
        pointCloudRef.current = null;
      }

      // Normalize coordinates to [-1, 1] range
      const normalized = normalizeCoordinates(embeddingData);

      if (!normalized) return;

      // Scale factor to fit in view (frustum size is 20)
      const scaleFactor = 8;

      // Create point data from embedding
      const pointData: PointData[] = normalized.normalized.map(
        (coord: number[]) => ({
          x: coord[0] * scaleFactor,
          y: coord[1] * scaleFactor,
          z: 0, // 2D embedding
          r: 0.7, // Default gray color
          g: 0.7,
          b: 0.7,
          size: 1.0,
          alpha: 1.0,
        }),
      );

      // Create point cloud with configurable dot size
      const pointCloud = createPointCloud(
        pointData,
        VISUALIZATION_CONFIG.UMAP_POINT_SIZE,
      );

      scene.add(pointCloud);
      pointCloudRef.current = pointCloud;
      setPointCloudVersion((v) => v + 1); // Trigger visualization update

      // Force a render after adding point cloud
      if (rendererRef.current && cameraRef.current) {
        rendererRef.current.render(scene, cameraRef.current);
      }
    },
    [],
  );

  // Create point cloud from embedding data (with on-demand loading)
  useEffect(() => {
    if (
      !dataset ||
      !selectedEmbedding ||
      !sceneRef.current ||
      !isOpen ||
      !sceneReady
    ) {
      return;
    }

    const scene = sceneRef.current;
    let cancelled = false;

    const loadAndBuild = async () => {
      // Check if embedding is already loaded
      let embeddingData = dataset.embeddings[selectedEmbedding];

      if (!embeddingData || embeddingData.length === 0) {
        if (!dataset.adapter) return;

        const toastId = toast.loading(
          `Loading ${selectedEmbedding.toUpperCase()} embedding...`,
        );

        try {
          let result: { name: string; data: number[][] } | null = null;

          if (dataset.adapter.mode === "local") {
            // Local datasets: load on main thread (File objects can't be sent to worker)
            result = await dataset.adapter.loadEmbedding(selectedEmbedding);
          } else {
            // S3 / custom S3: load in worker to avoid freezing UI
            const { getStandardizedDatasetWorker } = await import(
              "@/lib/workers/standardizedDatasetWorkerManager"
            );
            const worker = await getStandardizedDatasetWorker();

            result = await worker.loadEmbeddingFromS3(
              dataset.id,
              selectedEmbedding,
              dataset.metadata?.customS3BaseUrl,
            );
          }

          if (cancelled) {
            toast.dismiss(toastId);

            return;
          }

          if (result) {
            dataset.addEmbedding(result.name, result.data);
            embeddingData = result.data;
            toast.dismiss(toastId);
          } else {
            toast.update(toastId, {
              render: `Failed to load ${selectedEmbedding.toUpperCase()} embedding`,
              type: "error",
              isLoading: false,
              autoClose: 3000,
            });

            return;
          }
        } catch (error) {
          if (cancelled) {
            toast.dismiss(toastId);

            return;
          }
          toast.update(toastId, {
            render: `Error loading ${selectedEmbedding.toUpperCase()} embedding`,
            type: "error",
            isLoading: false,
            autoClose: 3000,
          });

          return;
        }
      }

      if (!cancelled && embeddingData && embeddingData.length > 0) {
        buildPointCloud(embeddingData, scene);
      }
    };

    loadAndBuild();

    return () => {
      cancelled = true;
    };
  }, [dataset, selectedEmbedding, isOpen, sceneReady, buildPointCloud]);

  // Update visualization colors based on mode (sync with main viewer)
  useEffect(() => {
    if (!pointCloudRef.current || !dataset || !selectedColumn || !sceneReady) {
      return;
    }

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

      // Determine which visualization mode to use
      const hasGeneMode = mode.includes("gene");
      const hasCelltypeMode = mode.includes("celltype");

      // Get cluster data
      const clusterData = dataset.clusters?.find(
        (c) => c.column === selectedColumn,
      );
      const isNumerical = clusterData?.type === "numerical";

      let result = null;

      if (
        hasGeneMode &&
        hasCelltypeMode &&
        selectedGene &&
        selectedCelltypes.size > 0
      ) {
        // Combined mode: gene expression on selected celltypes
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
        );
      } else if (hasGeneMode && selectedGene) {
        // Gene mode only
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
        );
      } else if (hasCelltypeMode) {
        // Clear gene expression ref when not in gene mode
        geneExpressionRef.current = null;
        // Celltype mode only
        result = isNumerical
          ? updateNumericalCelltypeVisualization(
              dataset,
              selectedColumn,
              alphaScale,
              sizeScale,
              numericalScaleMin,
              numericalScaleMax,
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

      // Update point cloud colors if we have a result
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
    mode,
    selectedGene,
    selectedColumn,
    selectedCelltypes,
    colorPalette,
    alphaScale,
    sizeScale,
    geneScaleMin,
    geneScaleMax,
    numericalScaleMin,
    numericalScaleMax,
    sceneReady,
    selectedEmbedding, // Re-run visualization when embedding changes
    pointCloudVersion,
  ]);

  // Mouse interaction handlers
  useEffect(() => {
    if (!containerRef.current || !isOpen) return;
    if (!rendererRef.current || !cameraRef.current || !pointCloudRef.current)
      return;

    const renderer = rendererRef.current;
    const camera = cameraRef.current;
    const pointCloud = pointCloudRef.current;

    let lastCheckTime = 0;
    const throttleDelay = 16; // ms (~60fps for smoother detection)

    const checkIntersections = () => {
      if (!pointCloud || !camera) return;

      raycasterRef.current.setFromCamera(mouseRef.current, camera);

      // Use a smaller threshold for better accuracy with orthographic camera
      // Threshold represents the maximum distance from point to ray in world units
      raycasterRef.current.params.Points!.threshold = 0.3;

      const intersects = raycasterRef.current.intersectObject(pointCloud);

      if (intersects.length > 0) {
        const intersect = intersects[0];
        const index = intersect.index!;

        hoveredPointRef.current = index;

        // Show tooltip
        showTooltip(intersect.point, index);
      } else {
        hoveredPointRef.current = null;
        hideTooltip();
      }
    };

    // Mouse move handler
    const handleMouseMove = (event: MouseEvent) => {
      const rect = renderer.domElement.getBoundingClientRect();

      mouseRef.current.x = ((event.clientX - rect.left) / rect.width) * 2 - 1;
      mouseRef.current.y = -((event.clientY - rect.top) / rect.height) * 2 + 1;

      // Throttle intersection checking
      const now = Date.now();

      if (now - lastCheckTime > throttleDelay) {
        lastCheckTime = now;
        checkIntersections();
      }
    };

    // Double-click handler
    const handleDoubleClick = () => {
      if (hoveredPointRef.current !== null && dataset?.clusters) {
        const index = hoveredPointRef.current;
        const clusterValue = clusterValuesRef.current[index];
        const clusterValueStr = String(clusterValue);

        // Only toggle celltype if it's not a numerical cluster
        if (!isNumericalClusterRef.current) {
          toggleCelltype(clusterValueStr);
        }
      }
    };

    // Add event listeners
    renderer.domElement.addEventListener("mousemove", handleMouseMove);
    renderer.domElement.addEventListener("dblclick", handleDoubleClick);

    return () => {
      renderer.domElement.removeEventListener("mousemove", handleMouseMove);
      renderer.domElement.removeEventListener("dblclick", handleDoubleClick);
    };
  }, [
    isOpen,
    dataset,
    toggleCelltype,
    getPointColor,
    showTooltip,
    hideTooltip,
  ]);

  // Check if embeddings are available (either loaded or available for lazy loading)
  const hasEmbeddings =
    (dataset?.allEmbeddingNames && dataset.allEmbeddingNames.length > 0) ||
    (dataset?.embeddings && Object.keys(dataset.embeddings).length > 0);

  return (
    <>
      {/* Main UMAP Panel */}
      <AnimatePresence>
        {isOpen && (
          <motion.div
            animate={{ opacity: 1, y: 0 }}
            className="!absolute bottom-4 right-4 md:bottom-8 md:right-8
                       w-[50vh] h-[50vh]
                       bg-black
                       rounded-3xl
                       shadow-2xl
                       border border-white/10
                       z-[9999]
                       "
            exit={{ opacity: 0, y: 100 }}
            initial={{ opacity: 0, y: 100 }}
            style={{ position: "absolute" }}
            transition={{ duration: 0.3, ease: "easeInOut" }}
          >
            {/* Floating Embedding Selector */}
            <div className="absolute top-3 left-3 z-50 flex items-center gap-2">
              <Select
                aria-label="Select embedding"
                className="w-32"
                classNames={{
                  trigger:
                    "bg-black/50 backdrop-blur-sm border border-white/20 h-8 min-h-8",
                  value: "text-white text-xs",
                  popoverContent: "bg-black/90 backdrop-blur-sm",
                }}
                selectedKeys={selectedEmbedding ? [selectedEmbedding] : []}
                size="sm"
                onChange={(e) => setSelectedEmbedding(e.target.value)}
              >
                {availableEmbeddings.map((key) => (
                  <SelectItem key={key}>{key.toUpperCase()}</SelectItem>
                ))}
              </Select>
            </div>

            {/* Floating Close Button */}
            <button
              aria-label="Close UMAP"
              className="absolute top-3 right-3 z-50
                         p-1.5 rounded-lg
                         bg-black/50 backdrop-blur-sm
                         border border-white/20
                         hover:bg-white/10 hover:border-white/30
                         transition-all duration-200"
              onClick={() => setIsOpen(false)}
            >
              <IoClose className="w-5 h-5 text-white/70 hover:text-white" />
            </button>

            {/* Three.js Scene Container */}
            <div ref={containerRef} className="w-full h-full relative" />
          </motion.div>
        )}
      </AnimatePresence>

      {/* Minimized Toggle Button */}
      <AnimatePresence>
        {!isOpen && (
          <motion.div
            animate={{ opacity: 1, scale: 1 }}
            className="!absolute bottom-4 right-4 md:bottom-8 md:right-8 z-[9999]"
            exit={{ opacity: 0, scale: 0.8 }}
            initial={{ opacity: 0, scale: 0.8 }}
            style={{ position: "absolute" }}
            transition={{ duration: 0.2 }}
          >
            <Tooltip
              className="bg-black/90 text-white text-xs px-3 py-1.5 rounded-lg"
              content={
                hasEmbeddings
                  ? "Open UMAP Visualization"
                  : "No embeddings available"
              }
              placement="left"
            >
              <button
                aria-label="Open UMAP"
                className={`p-4 bg-black rounded-full
                           border border-white/20
                           shadow-2xl
                           transition-all duration-200
                           ${
                             hasEmbeddings
                               ? "hover:bg-white/10 hover:border-white/30 hover:scale-110 cursor-pointer"
                               : "opacity-50 cursor-not-allowed"
                           }`}
                disabled={!hasEmbeddings}
                onClick={() => hasEmbeddings && setIsOpen(true)}
              >
                <TbChartDots3 className="w-6 h-6 text-white" />
              </button>
            </Tooltip>
          </motion.div>
        )}
      </AnimatePresence>
    </>
  );
}
