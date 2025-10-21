"use client";

import { useEffect, useRef } from "react";
import * as THREE from "three";
import { Spinner } from "@heroui/react";
import { initializeScene } from "@/lib/webgl/scene-manager";
import { createPointCloud } from "@/lib/webgl/point-cloud";
import { useSingleMoleculeStore } from "@/lib/stores/singleMoleculeStore";
import { useSingleMoleculeVisualizationStore } from "@/lib/stores/singleMoleculeVisualizationStore";

export function SingleMoleculeThreeScene() {
  const containerRef = useRef<HTMLDivElement>(null);
  const sceneRef = useRef<THREE.Scene | null>(null);
  const rendererRef = useRef<THREE.WebGLRenderer | null>(null);
  const cameraRef = useRef<THREE.Camera | null>(null);
  const controlsRef = useRef<any>(null);
  const pointCloudsRef = useRef<Map<string, THREE.Points>>(new Map());

  // Get dataset from store
  const dataset = useSingleMoleculeStore((state) => {
    const id = state.currentDatasetId;
    return id ? state.datasets.get(id) : null;
  });

  // Get visualization settings from store
  const { selectedGenes, globalScale, viewMode } = useSingleMoleculeVisualizationStore();

  // Effect 1: Scene initialization
  useEffect(() => {
    if (!containerRef.current || !dataset) return;

    console.log("[SingleMoleculeThreeScene] Initializing scene");
    console.log("Dataset:", dataset.name);
    console.log("Dimensions:", dataset.dimensions);
    console.log("Total molecules:", dataset.getMoleculeCount());
    console.log("Unique genes:", dataset.uniqueGenes.length);

    // Clear any existing canvas before creating new one
    if (containerRef.current) {
      while (containerRef.current.firstChild) {
        containerRef.current.removeChild(containerRef.current.firstChild);
      }
    }

    // Initialize Three.js scene
    const { scene, camera, renderer, controls, animate } = initializeScene(
      containerRef.current,
      dataset.dimensions
    );

    sceneRef.current = scene;
    rendererRef.current = renderer;
    cameraRef.current = camera;
    controlsRef.current = controls;

    // For 2D: disable rotation, enable panning with left click
    if (dataset.dimensions === 2 || viewMode === "2D") {
      controls.enableRotate = false;
      controls.mouseButtons = {
        LEFT: THREE.MOUSE.PAN,
        MIDDLE: THREE.MOUSE.DOLLY,
        RIGHT: THREE.MOUSE.PAN,
      };
    }

    // Start animation loop
    animate();

    // Cleanup
    return () => {
      console.log("[SingleMoleculeThreeScene] Cleaning up scene");
      controls.dispose();
      renderer.dispose();
      if (containerRef.current?.contains(renderer.domElement)) {
        containerRef.current.removeChild(renderer.domElement);
      }
      rendererRef.current = null;
      sceneRef.current = null;
    };
  }, [dataset, viewMode]);

  // Effect 2: Update camera controls when viewMode changes
  useEffect(() => {
    if (!controlsRef.current || !cameraRef.current) return;

    const controls = controlsRef.current;
    const camera = cameraRef.current;

    if (viewMode === "2D") {
      console.log("[SingleMoleculeThreeScene] Switching to 2D mode");

      // Reset camera to top-down view
      camera.position.set(0, 0, 200);
      camera.lookAt(0, 0, 0);

      // Reset controls target to origin
      controls.target.set(0, 0, 0);

      controls.enableRotate = false;
      controls.mouseButtons = {
        LEFT: THREE.MOUSE.PAN,
        MIDDLE: THREE.MOUSE.DOLLY,
        RIGHT: THREE.MOUSE.PAN,
      };

      controls.update();
    } else {
      console.log("[SingleMoleculeThreeScene] Switching to 3D mode");

      // Reset camera to 3D perspective view
      camera.position.set(150, 150, 150);
      camera.lookAt(0, 0, 0);

      // Reset controls target to origin
      controls.target.set(0, 0, 0);

      controls.enableRotate = true;
      controls.mouseButtons = {
        LEFT: THREE.MOUSE.ROTATE,
        MIDDLE: THREE.MOUSE.DOLLY,
        RIGHT: THREE.MOUSE.PAN,
      };

      controls.update();
    }
  }, [viewMode]);

  // Effect 3: Update point clouds based on selected genes
  useEffect(() => {
    if (!sceneRef.current || !dataset) return;

    console.log("[SingleMoleculeThreeScene] Updating point clouds");
    console.log("Selected genes:", Array.from(selectedGenes.keys()));

    const scene = sceneRef.current;
    const currentPointClouds = pointCloudsRef.current;

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
        try {
          // Get coordinates for this gene (async for lazy loading from S3)
          const coords = await dataset.getCoordinatesByGene(gene);
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
            positions.push(coords[i] * 100, coords[i + 1] * 100, coords[i + 2] * 100);
          }

          // Create point cloud with single color
          const geometry = new THREE.BufferGeometry();
          geometry.setAttribute(
            "position",
            new THREE.Float32BufferAttribute(positions, 3)
          );

          // Parse HSL color and convert to RGB
          const color = new THREE.Color(geneViz.color);
          const colors = new Float32Array(moleculeCount * 3);
          for (let i = 0; i < moleculeCount; i++) {
            colors[i * 3] = color.r;
            colors[i * 3 + 1] = color.g;
            colors[i * 3 + 2] = color.b;
          }
          geometry.setAttribute("color", new THREE.BufferAttribute(colors, 3));

          // Create material (increased base size)
          const material = new THREE.PointsMaterial({
            size: geneViz.localScale * globalScale * 2.0,
            vertexColors: true,
            transparent: true,
            opacity: 1.0,
            sizeAttenuation: false,
          });

          pointCloud = new THREE.Points(geometry, material);
          scene.add(pointCloud);
          currentPointClouds.set(gene, pointCloud);

          console.log(`  ✅ Point cloud created with ${moleculeCount} molecules`);
        } else {
          // Update existing point cloud color and size
          const material = pointCloud.material as THREE.PointsMaterial;
          material.size = geneViz.localScale * globalScale * 2.0;

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
        }
        } catch (error) {
          console.error(`Error creating point cloud for gene ${gene}:`, error);
        }
      }
    };

    // Call async function
    updatePointClouds();
  }, [dataset, selectedGenes, globalScale]);

  if (!dataset) {
    return (
      <div className="w-full h-full flex items-center justify-center bg-black">
        <div className="text-center">
          <Spinner size="lg" color="primary" />
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
