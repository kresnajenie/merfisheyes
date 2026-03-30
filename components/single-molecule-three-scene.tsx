"use client";

import { useEffect, useRef } from "react";
import * as THREE from "three";
import { Spinner } from "@heroui/react";
import { toast } from "react-toastify";

import { initializeScene } from "@/lib/webgl/scene-manager";
import {
  usePanelSingleMoleculeStore,
  usePanelSingleMoleculeVisualizationStore,
} from "@/lib/hooks/usePanelStores";
import { VISUALIZATION_CONFIG } from "@/lib/config/visualization.config";
import type { MoleculeShape } from "@/lib/stores/createSingleMoleculeVisualizationStore";


// Create solid circular sprite texture for points
function createCircleTexture(): THREE.Texture {
  const canvas = document.createElement("canvas");
  const size = 64;

  canvas.width = size;
  canvas.height = size;

  const context = canvas.getContext("2d");

  if (!context) throw new Error("Could not get 2D context");

  // Draw a solid white circle with sharp edges
  context.fillStyle = "white";
  context.beginPath();
  context.arc(size / 2, size / 2, size / 2, 0, Math.PI * 2);
  context.fill();

  const texture = new THREE.CanvasTexture(canvas);

  texture.needsUpdate = true;

  return texture;
}

// Create solid square sprite texture for unassigned points
function createSquareTexture(): THREE.Texture {
  const canvas = document.createElement("canvas");
  const size = 64;

  canvas.width = size;
  canvas.height = size;

  const context = canvas.getContext("2d");

  if (!context) throw new Error("Could not get 2D context");

  // Draw a solid white square
  context.fillStyle = "white";
  context.fillRect(4, 4, size - 8, size - 8);

  const texture = new THREE.CanvasTexture(canvas);

  texture.needsUpdate = true;

  return texture;
}

// Point cloud key helpers for assigned/unassigned separation
const ASSIGNED_KEY_SUFFIX = ":assigned";
const UNASSIGNED_KEY_SUFFIX = ":unassigned";

function assignedKey(gene: string): string {
  return gene + ASSIGNED_KEY_SUFFIX;
}

function unassignedKey(gene: string): string {
  return gene + UNASSIGNED_KEY_SUFFIX;
}

function geneFromKey(key: string): string {
  if (key.endsWith(ASSIGNED_KEY_SUFFIX)) {
    return key.slice(0, -ASSIGNED_KEY_SUFFIX.length);
  }
  if (key.endsWith(UNASSIGNED_KEY_SUFFIX)) {
    return key.slice(0, -UNASSIGNED_KEY_SUFFIX.length);
  }

  return key;
}

export function SingleMoleculeThreeScene() {
  const containerRef = useRef<HTMLDivElement>(null);
  const sceneRef = useRef<THREE.Scene | null>(null);
  const rendererRef = useRef<THREE.WebGLRenderer | null>(null);
  const cameraRef = useRef<THREE.Camera | null>(null);
  const controlsRef = useRef<any>(null);
  const pointCloudsRef = useRef<Map<string, THREE.Points>>(new Map());
  const circleTextureRef = useRef<THREE.Texture | null>(null);
  const squareTextureRef = useRef<THREE.Texture | null>(null);
  const lastDatasetIdRef = useRef<string | null>(null);
  const lastViewModeRef = useRef<string | null>(null);
  const baselineCameraDistanceRef = useRef<number | null>(null);
  const animationFrameIdRef = useRef<number | null>(null);
  const selectedGenesRef = useRef<Map<string, any>>(new Map());
  const globalScaleRef = useRef<number>(1);

  // Get dataset from store - using stable selector to prevent re-renders
  const currentDatasetId = usePanelSingleMoleculeStore(
    (state) => state.currentDatasetId,
  );
  const dataset = usePanelSingleMoleculeStore((state) =>
    currentDatasetId ? state.datasets.get(currentDatasetId) : null,
  );

  // Get visualization settings from store
  const { selectedGenes, globalScale, viewMode, showAssigned, showUnassigned } =
    usePanelSingleMoleculeVisualizationStore();

  // Keep refs in sync with store values
  selectedGenesRef.current = selectedGenes;
  globalScaleRef.current = globalScale;

  // Debug: Log when component re-renders
  console.log("[SingleMoleculeThreeScene] Component render", {
    datasetId: dataset?.id,
    datasetName: dataset?.name,
    selectedGenesCount: selectedGenes.size,
    lastDatasetId: lastDatasetIdRef.current,
    sceneExists: !!sceneRef.current,
  });

  // Effect 1: Scene initialization
  useEffect(() => {
    console.log("[SingleMoleculeThreeScene] Effect triggered", {
      hasContainer: !!containerRef.current,
      hasDataset: !!dataset,
      datasetId: dataset?.id,
      lastDatasetId: lastDatasetIdRef.current,
      hasScene: !!sceneRef.current,
    });

    if (!containerRef.current || !dataset) return;

    // Only re-initialize if dataset ID AND viewMode are unchanged AND scene exists
    // This prevents unnecessary re-init during React Strict Mode double-render
    if (
      lastDatasetIdRef.current === dataset.id &&
      lastViewModeRef.current === viewMode &&
      sceneRef.current
    ) {
      console.log(
        "[SingleMoleculeThreeScene] Dataset ID and viewMode unchanged and scene exists, skipping re-initialization",
        {
          datasetId: dataset.id,
          lastDatasetId: lastDatasetIdRef.current,
          viewMode,
          lastViewMode: lastViewModeRef.current,
        },
      );

      return;
    }

    // If scene was cleaned up but dataset ID is same, we still need to recreate
    if (lastDatasetIdRef.current === dataset.id && !sceneRef.current) {
      console.log(
        "[SingleMoleculeThreeScene] Scene was cleaned up but dataset ID unchanged, will recreate scene",
        {
          datasetId: dataset.id,
        },
      );
    }

    console.log("[SingleMoleculeThreeScene] Initializing scene");
    console.log("Dataset:", dataset.name);
    console.log("Dataset ID:", dataset.id);
    console.log("Previous Dataset ID:", lastDatasetIdRef.current);
    console.log("View Mode:", viewMode);
    console.log("Previous View Mode:", lastViewModeRef.current);
    console.log("Dimensions:", dataset.dimensions);
    console.log("Total molecules:", dataset.getMoleculeCount());
    console.log("Unique genes:", dataset.uniqueGenes.length);

    // Update last dataset ID and viewMode
    lastDatasetIdRef.current = dataset.id;
    lastViewModeRef.current = viewMode;

    // Clear any existing canvas before creating new one
    if (containerRef.current) {
      while (containerRef.current.firstChild) {
        containerRef.current.removeChild(containerRef.current.firstChild);
      }
    }

    // Initialize Three.js scene with viewMode preference (not dataset dimensions)
    const { scene, camera, renderer, controls } = initializeScene(
      containerRef.current,
      { is2D: viewMode === "2D" },
    );

    sceneRef.current = scene;
    rendererRef.current = renderer;
    cameraRef.current = camera;
    controlsRef.current = controls;

    // Create textures for points (reuse across all point clouds)
    if (!circleTextureRef.current) {
      circleTextureRef.current = createCircleTexture();
    }
    if (!squareTextureRef.current) {
      squareTextureRef.current = createSquareTexture();
    }

    // Set baseline camera distance (initial zoomed-out distance)
    if (baselineCameraDistanceRef.current === null) {
      baselineCameraDistanceRef.current = camera.position.distanceTo(
        new THREE.Vector3(0, 0, 0),
      );
    }

    // Custom animation loop with dynamic point size scaling
    const customAnimate = () => {
      animationFrameIdRef.current = requestAnimationFrame(customAnimate);

      controls.update();

      // Calculate zoom factor (k) and update point sizes
      const currentDistance = camera.position.distanceTo(
        new THREE.Vector3(0, 0, 0),
      );
      const zoomFactor = baselineCameraDistanceRef.current! / currentDistance;

      // Power-law scaling: s(k) = clamp(s0 * k^alpha, sMin, sMax)
      const alpha = 0.8; // Power-law exponent (0.7-0.9 range)
      const s0 = VISUALIZATION_CONFIG.SINGLE_MOLECULE_POINT_BASE_SIZE;
      const sMin = s0; // Minimum = base size (zoomed out)
      const sMax = 200; // Maximum = 200 (zoomed in)

      const scaledSize = Math.pow(zoomFactor, alpha) * s0;
      const clampedSize = Math.max(sMin, Math.min(sMax, scaledSize));

      // Update all point cloud sizes
      pointCloudsRef.current.forEach((pointCloud, key) => {
        const gene = geneFromKey(key);
        const geneViz = selectedGenesRef.current.get(gene);

        if (geneViz) {
          const material = pointCloud.material as THREE.PointsMaterial;
          const isUnassigned = key.endsWith(UNASSIGNED_KEY_SUFFIX);
          const scale = isUnassigned
            ? geneViz.unassignedLocalScale
            : geneViz.localScale;

          material.size = scale * globalScaleRef.current * clampedSize;
        }
      });

      renderer.render(scene, camera);
    };

    // Start custom animation loop
    customAnimate();

    // Cleanup
    return () => {
      console.log("[SingleMoleculeThreeScene] Cleaning up scene", {
        lastDatasetId: lastDatasetIdRef.current,
        reason: "Effect cleanup triggered - dataset dependency changed",
      });

      // Cancel animation frame
      if (animationFrameIdRef.current !== null) {
        cancelAnimationFrame(animationFrameIdRef.current);
        animationFrameIdRef.current = null;
      }

      // Clear all point clouds
      console.log(
        "[SingleMoleculeThreeScene] Clearing point clouds, count:",
        pointCloudsRef.current.size,
      );
      pointCloudsRef.current.forEach((pointCloud, key) => {
        console.log(`  Disposing point cloud: ${key}`);
        pointCloud.geometry.dispose();
        if (Array.isArray(pointCloud.material)) {
          pointCloud.material.forEach((m) => m.dispose());
        } else {
          pointCloud.material.dispose();
        }
      });
      pointCloudsRef.current.clear();
      console.log(
        "[SingleMoleculeThreeScene] Point clouds cleared, new count:",
        pointCloudsRef.current.size,
      );

      // Properly dispose of Three.js resources
      if (controls) controls.dispose();
      if (renderer) {
        // Force lose WebGL context before disposing
        const gl = renderer.getContext();

        if (gl) {
          const loseContextExt = gl.getExtension("WEBGL_lose_context");

          if (loseContextExt) loseContextExt.loseContext();
        }
        renderer.dispose();

        // Remove canvas from DOM
        if (containerRef.current?.contains(renderer.domElement)) {
          containerRef.current.removeChild(renderer.domElement);
        }
      }

      rendererRef.current = null;
      sceneRef.current = null;
      // DON'T reset lastDatasetIdRef here - keep it so we can skip re-init on React Strict Mode double-render
      console.log("[SingleMoleculeThreeScene] Cleanup complete");
    };
  }, [dataset, viewMode]); // Re-initialize scene when viewMode changes to get correct controls (TrackballControls for 3D, OrbitControls for 2D)

  // Effect 3: Update point clouds based on selected genes
  useEffect(() => {
    if (!sceneRef.current || !dataset) {
      console.log("[SingleMoleculeThreeScene] Skipping point cloud update:", {
        hasScene: !!sceneRef.current,
        hasDataset: !!dataset,
      });

      return;
    }

    console.log("[SingleMoleculeThreeScene] === POINT CLOUD UPDATE START ===");
    console.log("Selected genes:", Array.from(selectedGenes.keys()));
    console.log("Current point clouds count:", pointCloudsRef.current.size);
    console.log(
      "Current point clouds genes:",
      Array.from(pointCloudsRef.current.keys()),
    );

    const scene = sceneRef.current;
    const currentPointClouds = pointCloudsRef.current;

    // Create abort flag for this effect run
    let isCancelled = false;

    // Remove point clouds for unselected genes
    for (const [key, pointCloud] of currentPointClouds.entries()) {
      const gene = geneFromKey(key);

      if (!selectedGenes.has(gene)) {
        console.log(`Removing point cloud: ${key}`);
        scene.remove(pointCloud);
        pointCloud.geometry.dispose();
        if (Array.isArray(pointCloud.material)) {
          pointCloud.material.forEach((m) => m.dispose());
        } else {
          pointCloud.material.dispose();
        }
        currentPointClouds.delete(key);
      }
    }

    // Helper: create a point cloud from coordinates
    const createPointCloud = (
      coords: number[],
      geneViz: { color: string; localScale: number },
      texture: THREE.Texture | null,
      renderOrder: number = 0,
    ): THREE.Points => {
      const moleculeCount = coords.length / 3;
      const positions: number[] = [];

      for (let i = 0; i < coords.length; i += 3) {
        positions.push(
          coords[i] * 100,
          coords[i + 1] * 100,
          coords[i + 2] * 100,
        );
      }

      const geometry = new THREE.BufferGeometry();

      geometry.setAttribute(
        "position",
        new THREE.Float32BufferAttribute(positions, 3),
      );

      const color = new THREE.Color(geneViz.color);
      const colors = new Float32Array(moleculeCount * 3);

      for (let i = 0; i < moleculeCount; i++) {
        colors[i * 3] = color.r;
        colors[i * 3 + 1] = color.g;
        colors[i * 3 + 2] = color.b;
      }
      geometry.setAttribute("color", new THREE.BufferAttribute(colors, 3));

      const material = new THREE.PointsMaterial({
        size:
          geneViz.localScale *
          globalScale *
          VISUALIZATION_CONFIG.SINGLE_MOLECULE_POINT_BASE_SIZE,
        vertexColors: true,
        transparent: true,
        opacity: 1.0,
        sizeAttenuation: false,
        map: texture,
        alphaTest: 0.5,
      });

      const pointCloud = new THREE.Points(geometry, material);

      pointCloud.renderOrder = renderOrder;

      return pointCloud;
    };

    // Helper: update an existing point cloud's color and size
    const updatePointCloudAppearance = (
      pointCloud: THREE.Points,
      geneViz: { color: string; localScale: number },
    ) => {
      const material = pointCloud.material as THREE.PointsMaterial;

      material.size =
        geneViz.localScale *
        globalScale *
        VISUALIZATION_CONFIG.SINGLE_MOLECULE_POINT_BASE_SIZE;

      const color = new THREE.Color(geneViz.color);
      const colorAttr = pointCloud.geometry.getAttribute("color");
      const colors = colorAttr.array as Float32Array;
      const moleculeCount = colors.length / 3;

      for (let i = 0; i < moleculeCount; i++) {
        colors[i * 3] = color.r;
        colors[i * 3 + 1] = color.g;
        colors[i * 3 + 2] = color.b;
      }
      colorAttr.needsUpdate = true;
    };

    // Helper: remove and dispose a point cloud by key
    const removePointCloud = (key: string) => {
      const pc = currentPointClouds.get(key);

      if (pc) {
        scene.remove(pc);
        pc.geometry.dispose();
        if (Array.isArray(pc.material)) {
          pc.material.forEach((m) => m.dispose());
        } else {
          pc.material.dispose();
        }
        currentPointClouds.delete(key);
      }
    };

    // Add/update point clouds for selected genes (async)
    const updatePointClouds = async () => {
      for (const [gene, geneViz] of selectedGenes.entries()) {
        if (isCancelled) {
          console.log(
            `[SingleMoleculeThreeScene] Update cancelled for gene: ${gene}`,
          );

          return;
        }

        const aKey = assignedKey(gene);
        const uKey = unassignedKey(gene);
        const assignedExists = currentPointClouds.has(aKey);
        const toastId = `loading-gene-${gene}`;

        try {
          // Only show loading toast for new genes
          if (!assignedExists) {
            toast.loading(`Loading ${gene}...`, {
              toastId,
              position: "bottom-left",
              autoClose: false,
            });
          }

          // Get assigned coordinates
          const coords = await dataset.getCoordinatesByGene(gene);

          if (isCancelled) {
            toast.dismiss(toastId);

            return;
          }

          const moleculeCount = coords.length / 3;

          console.log(`Creating/updating point cloud for gene: ${gene}`);
          console.log(`  Assigned molecules: ${moleculeCount.toLocaleString()}`);
          console.log(`  Color: ${geneViz.color}`);

          // Helper to get texture for a shape
          const getTexture = (shape: MoleculeShape) =>
            shape === "square"
              ? squareTextureRef.current
              : circleTextureRef.current;

          // --- Assigned point cloud ---
          if (geneViz.showAssigned) {
            let pointCloud = currentPointClouds.get(aKey);

            if (!pointCloud) {
              pointCloud = createPointCloud(
                coords,
                geneViz,
                getTexture(geneViz.assignedShape),
                0,
              );

              if (isCancelled) {
                pointCloud.geometry.dispose();
                if (pointCloud.material instanceof THREE.Material) {
                  pointCloud.material.dispose();
                }
                toast.dismiss(toastId);

                return;
              }

              scene.add(pointCloud);
              currentPointClouds.set(aKey, pointCloud);

              console.log(
                `  ✅ Assigned point cloud created with ${moleculeCount} molecules`,
              );
            } else {
              updatePointCloudAppearance(pointCloud, geneViz);
              // Update texture if shape changed
              const mat = pointCloud.material as THREE.PointsMaterial;

              mat.map = getTexture(geneViz.assignedShape);
              mat.needsUpdate = true;
              console.log(`  ✅ Assigned point cloud updated`);
            }
          } else {
            removePointCloud(aKey);
          }

          // --- Unassigned point cloud ---
          const shouldShowUnassigned =
            dataset.hasUnassigned && showUnassigned && geneViz.showUnassigned;

          if (shouldShowUnassigned) {
            const unassignedExists = currentPointClouds.has(uKey);
            const uViz = {
              color: geneViz.unassignedColor,
              localScale: geneViz.unassignedLocalScale,
            };

            if (!unassignedExists) {
              // Load unassigned coordinates
              const uCoords =
                await dataset.getUnassignedCoordinatesByGene(gene);

              if (isCancelled) {
                toast.dismiss(toastId);

                return;
              }

              if (uCoords.length > 0) {
                const uPointCloud = createPointCloud(
                  uCoords,
                  uViz,
                  getTexture(geneViz.unassignedShape),
                  -1, // Render behind assigned
                );

                if (isCancelled) {
                  uPointCloud.geometry.dispose();
                  if (uPointCloud.material instanceof THREE.Material) {
                    uPointCloud.material.dispose();
                  }
                  toast.dismiss(toastId);

                  return;
                }

                scene.add(uPointCloud);
                currentPointClouds.set(uKey, uPointCloud);

                const uMoleculeCount = uCoords.length / 3;

                console.log(
                  `  ✅ Unassigned point cloud created with ${uMoleculeCount} molecules`,
                );
              }
            } else {
              // Update existing unassigned cloud color/size/texture
              const uPC = currentPointClouds.get(uKey)!;

              updatePointCloudAppearance(uPC, uViz);
              const uMat = uPC.material as THREE.PointsMaterial;

              uMat.map = getTexture(geneViz.unassignedShape);
              uMat.needsUpdate = true;
            }
          } else {
            // Remove unassigned cloud if toggle is off or dataset has none
            removePointCloud(uKey);
          }

          // Update toast
          if (!assignedExists) {
            toast.update(toastId, {
              render: `${gene}: ${moleculeCount.toLocaleString()} molecules loaded`,
              type: "success",
              isLoading: false,
              autoClose: 3000,
              position: "bottom-left",
            });
          } else {
            toast.dismiss(toastId);
          }
        } catch (error) {
          console.error(`Error creating point cloud for gene ${gene}:`, error);

          if (assignedExists) {
            toast.error(`Failed to update ${gene}`, {
              toastId,
              position: "bottom-left",
              autoClose: 5000,
            });
          } else {
            toast.update(toastId, {
              render: `Failed to load ${gene}`,
              type: "error",
              isLoading: false,
              autoClose: 5000,
              position: "bottom-left",
            });
          }
        }
      }
    };

    // Call async function
    updatePointClouds().then(() => {
      console.log(
        "[SingleMoleculeThreeScene] === POINT CLOUD UPDATE COMPLETE ===",
      );
      console.log("Final point clouds count:", pointCloudsRef.current.size);
      console.log(
        "Final point clouds keys:",
        Array.from(pointCloudsRef.current.keys()),
      );
    });

    // Cleanup: cancel the async operation if effect is cleaned up
    return () => {
      console.log("[SingleMoleculeThreeScene] Cancelling point cloud updates");
      isCancelled = true;
    };
  }, [dataset, selectedGenes, globalScale, viewMode, showAssigned, showUnassigned]);

  if (!dataset) {
    return (
      <div className="w-full h-full flex items-center justify-center bg-black">
        <div className="text-center">
          <Spinner color="primary" size="lg" />
          <p className="mt-4 text-white">No dataset loaded</p>
        </div>
      </div>
    );
  }

  return (
    <>
      <div
        ref={containerRef}
        className="absolute inset-0 w-full h-full"
        style={{ margin: 0, padding: 0 }}
      />

    </>
  );
}
