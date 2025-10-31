"use client";

import { useEffect, useRef } from "react";
import * as THREE from "three";
import { Spinner } from "@heroui/react";
import { toast } from "react-toastify";
import { OrbitControls } from "three/examples/jsm/controls/OrbitControls";
import { TrackballControls } from "three/examples/jsm/controls/TrackballControls";

import { initializeScene } from "@/lib/webgl/scene-manager";
import { useSingleMoleculeStore } from "@/lib/stores/singleMoleculeStore";
import { useSingleMoleculeVisualizationStore } from "@/lib/stores/singleMoleculeVisualizationStore";
import { VISUALIZATION_CONFIG } from "@/lib/config/visualization.config";

// Create circular sprite texture for points
function createCircleTexture(): THREE.Texture {
  const canvas = document.createElement("canvas");
  const size = 64;
  canvas.width = size;
  canvas.height = size;

  const context = canvas.getContext("2d");
  if (!context) throw new Error("Could not get 2D context");

  const gradient = context.createRadialGradient(
    size / 2,
    size / 2,
    0,
    size / 2,
    size / 2,
    size / 2,
  );
  gradient.addColorStop(0, "rgba(255,255,255,1)");
  gradient.addColorStop(0.5, "rgba(255,255,255,1)");
  gradient.addColorStop(1, "rgba(255,255,255,0)");

  context.fillStyle = gradient;
  context.fillRect(0, 0, size, size);

  const texture = new THREE.CanvasTexture(canvas);
  texture.needsUpdate = true;

  return texture;
}

export function SingleMoleculeThreeScene() {
  const containerRef = useRef<HTMLDivElement>(null);
  const sceneRef = useRef<THREE.Scene | null>(null);
  const rendererRef = useRef<THREE.WebGLRenderer | null>(null);
  const cameraRef = useRef<THREE.Camera | null>(null);
  const controlsRef = useRef<any>(null);
  const pointCloudsRef = useRef<Map<string, THREE.Points>>(new Map());
  const circleTextureRef = useRef<THREE.Texture | null>(null);
  const lastDatasetIdRef = useRef<string | null>(null);

  // Get dataset from store - using stable selector to prevent re-renders
  const currentDatasetId = useSingleMoleculeStore(
    (state) => state.currentDatasetId,
  );
  const dataset = useSingleMoleculeStore((state) =>
    currentDatasetId ? state.datasets.get(currentDatasetId) : null,
  );

  // Get visualization settings from store
  const { selectedGenes, globalScale, viewMode } =
    useSingleMoleculeVisualizationStore();

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

    // Only re-initialize if dataset ID actually changed AND scene exists
    // This prevents unnecessary re-init during React Strict Mode double-render
    if (lastDatasetIdRef.current === dataset.id && sceneRef.current) {
      console.log(
        "[SingleMoleculeThreeScene] Dataset ID unchanged and scene exists, skipping re-initialization",
        {
          datasetId: dataset.id,
          lastDatasetId: lastDatasetIdRef.current,
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
    console.log("Dimensions:", dataset.dimensions);
    console.log("Total molecules:", dataset.getMoleculeCount());
    console.log("Unique genes:", dataset.uniqueGenes.length);

    // Update last dataset ID
    lastDatasetIdRef.current = dataset.id;

    // Clear any existing canvas before creating new one
    if (containerRef.current) {
      while (containerRef.current.firstChild) {
        containerRef.current.removeChild(containerRef.current.firstChild);
      }
    }

    // Initialize Three.js scene with viewMode preference (not dataset dimensions)
    const { scene, camera, renderer, controls, animate } = initializeScene(
      containerRef.current,
      { is2D: viewMode === "2D" },
    );

    sceneRef.current = scene;
    rendererRef.current = renderer;
    cameraRef.current = camera;
    controlsRef.current = controls;

    // Create circle texture for points (reuse across all point clouds)
    if (!circleTextureRef.current) {
      circleTextureRef.current = createCircleTexture();
    }

    // Start animation loop
    animate();

    // Cleanup
    return () => {
      console.log("[SingleMoleculeThreeScene] Cleaning up scene", {
        lastDatasetId: lastDatasetIdRef.current,
        reason: "Effect cleanup triggered - dataset dependency changed",
      });

      // Properly dispose of Three.js resources
      if (controls) controls.dispose();
      if (renderer) {
        // Force lose WebGL context before disposing
        const gl = renderer.getContext();
        if (gl) {
          const loseContextExt = gl.getExtension('WEBGL_lose_context');
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
  }, [dataset]); // Removed viewMode to prevent unnecessary re-initialization

  // Effect 2: Update camera when viewMode changes
  // Note: Controls are created in initializeScene, we just adjust camera here
  useEffect(() => {
    if (!controlsRef.current || !cameraRef.current) return;

    console.log("[SingleMoleculeThreeScene] Adjusting camera for mode:", viewMode);

    const camera = cameraRef.current;
    const controls = controlsRef.current;

    if (viewMode === "2D") {
      console.log("[SingleMoleculeThreeScene] Setting 2D camera position");

      // Reset camera to top-down view
      camera.position.set(0, 0, 200);
      camera.lookAt(0, 0, 0);

      // Ensure controls target is at origin
      controls.target.set(0, 0, 0);
      controls.update();

      console.log("[SingleMoleculeThreeScene] 2D camera configured:", {
        cameraPosition: camera.position.toArray(),
        controlsTarget: controls.target.toArray(),
      });
    } else {
      console.log("[SingleMoleculeThreeScene] Setting 3D camera position");

      // Reset camera to 3D perspective view
      camera.position.set(150, 150, 150);
      camera.lookAt(0, 0, 0);

      // Ensure controls target is at origin
      controls.target.set(0, 0, 0);
      controls.update();

      console.log("[SingleMoleculeThreeScene] 3D camera configured");
    }
  }, [viewMode]);

  // Effect 3: Update point clouds based on selected genes
  useEffect(() => {
    if (!sceneRef.current || !dataset) return;

    console.log("[SingleMoleculeThreeScene] Updating point clouds");
    console.log("Selected genes:", Array.from(selectedGenes.keys()));

    const scene = sceneRef.current;
    const currentPointClouds = pointCloudsRef.current;

    // Create abort flag for this effect run
    let isCancelled = false;

    // Remove point clouds for unselected genes
    for (const [gene, pointCloud] of currentPointClouds.entries()) {
      if (!selectedGenes.has(gene)) {
        console.log(`Removing point cloud for gene: ${gene}`);
        scene.remove(pointCloud);
        pointCloud.geometry.dispose();
        if (Array.isArray(pointCloud.material)) {
          pointCloud.material.forEach((m) => m.dispose());
        } else {
          pointCloud.material.dispose();
        }
        currentPointClouds.delete(gene);
      }
    }

    // Add/update point clouds for selected genes (async)
    const updatePointClouds = async () => {
      for (const [gene, geneViz] of selectedGenes.entries()) {
        // Check if cancelled before processing each gene
        if (isCancelled) {
          console.log(`[SingleMoleculeThreeScene] Update cancelled for gene: ${gene}`);
          return;
        }

        // Create a unique toast ID for this gene
        const toastId = `loading-gene-${gene}`;

        try {
          // Show loading toast
          toast.loading(`Loading ${gene}...`, {
            toastId,
            position: "bottom-left",
            autoClose: false,
          });

          // Get coordinates for this gene (async for lazy loading from S3)
          const coords = await dataset.getCoordinatesByGene(gene);

          // Check again after async operation
          if (isCancelled) {
            console.log(`[SingleMoleculeThreeScene] Update cancelled after loading gene: ${gene}`);
            toast.dismiss(toastId);
            return;
          }

          const moleculeCount = coords.length / 3; // Each molecule has x, y, z

          console.log(`Creating/updating point cloud for gene: ${gene}`);
          console.log(`  Molecules: ${moleculeCount.toLocaleString()}`);
          console.log(`  Color: ${geneViz.color}`);
          console.log(`  Local scale: ${geneViz.localScale}`);

          // Check if point cloud already exists
          let pointCloud = currentPointClouds.get(gene);

          if (!pointCloud) {
            // Create new point cloud
            const positions: number[] = [];

            // Extract coordinates (already normalized to [-1, 1]) and scale by 100
            for (let i = 0; i < coords.length; i += 3) {
              positions.push(
                coords[i] * 100,
                coords[i + 1] * 100,
                coords[i + 2] * 100,
              );
            }

            // Create point cloud with single color
            const geometry = new THREE.BufferGeometry();

            geometry.setAttribute(
              "position",
              new THREE.Float32BufferAttribute(positions, 3),
            );

            // Parse HSL color and convert to RGB
            const color = new THREE.Color(geneViz.color);
            const colors = new Float32Array(moleculeCount * 3);

            for (let i = 0; i < moleculeCount; i++) {
              colors[i * 3] = color.r;
              colors[i * 3 + 1] = color.g;
              colors[i * 3 + 2] = color.b;
            }
            geometry.setAttribute(
              "color",
              new THREE.BufferAttribute(colors, 3),
            );

            // Create material with circular texture
            const material = new THREE.PointsMaterial({
              size: geneViz.localScale * globalScale * VISUALIZATION_CONFIG.SINGLE_MOLECULE_POINT_BASE_SIZE,
              vertexColors: true,
              transparent: true,
              opacity: 1.0,
              sizeAttenuation: false,
              map: circleTextureRef.current,
              alphaTest: 0.5,
            });

            pointCloud = new THREE.Points(geometry, material);

            // Final check before adding to scene
            if (isCancelled) {
              console.log(`[SingleMoleculeThreeScene] Cancelled before adding point cloud for: ${gene}`);
              pointCloud.geometry.dispose();
              pointCloud.material.dispose();
              toast.dismiss(toastId);
              return;
            }

            scene.add(pointCloud);
            currentPointClouds.set(gene, pointCloud);

            console.log(
              `  ✅ Point cloud created with ${moleculeCount} molecules`,
            );

            // Update toast to success
            toast.update(toastId, {
              render: `${gene}: ${moleculeCount.toLocaleString()} molecules loaded`,
              type: "success",
              isLoading: false,
              autoClose: 3000,
              position: "bottom-left",
            });
          } else {
            // Update existing point cloud color and size
            const material = pointCloud.material as THREE.PointsMaterial;

            material.size = geneViz.localScale * globalScale * VISUALIZATION_CONFIG.SINGLE_MOLECULE_POINT_BASE_SIZE;

            // Update colors
            const color = new THREE.Color(geneViz.color);
            const colorAttr = pointCloud.geometry.getAttribute("color");
            const colors = colorAttr.array as Float32Array;

            for (let i = 0; i < moleculeCount; i++) {
              colors[i * 3] = color.r;
              colors[i * 3 + 1] = color.g;
              colors[i * 3 + 2] = color.b;
            }
            colorAttr.needsUpdate = true;

            console.log(`  ✅ Point cloud updated`);

            // Dismiss loading toast for already-cached genes
            toast.dismiss(toastId);
          }
        } catch (error) {
          console.error(`Error creating point cloud for gene ${gene}:`, error);

          // Update toast to error
          toast.update(toastId, {
            render: `Failed to load ${gene}`,
            type: "error",
            isLoading: false,
            autoClose: 5000,
            position: "bottom-left",
          });
        }
      }
    };

    // Call async function
    updatePointClouds();

    // Cleanup: cancel the async operation if effect is cleaned up
    return () => {
      console.log("[SingleMoleculeThreeScene] Cancelling point cloud updates");
      isCancelled = true;
    };
  }, [dataset, selectedGenes, globalScale]);

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
    <div
      ref={containerRef}
      className="fixed inset-0 w-full h-full"
      style={{ margin: 0, padding: 0 }}
    />
  );
}
