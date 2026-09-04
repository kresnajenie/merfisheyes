"use client";

import { useEffect, useRef, useState } from "react";
import * as THREE from "three";
import { Spinner } from "@heroui/react";
import { toast } from "react-toastify";

import { initializeScene } from "@/lib/webgl/scene-manager";
import {
  appendToSmPointCloud,
  createSmPointCloud,
  createStreamingSmPointCloud,
  updateDotSize,
  updatePointCloudColor,
  updatePointCloudShape,
  SM_DOT_SIZE_FACTOR,
} from "@/lib/webgl/point-cloud";
import {
  usePanelSingleMoleculeStore,
  usePanelSingleMoleculeVisualizationStore,
} from "@/lib/hooks/usePanelStores";
import { VISUALIZATION_CONFIG } from "@/lib/config/visualization.config";
import {
  markSceneReady,
  markRenderComplete,
  markGeneRendered,
} from "@/lib/utils/test-hooks";
import type { MoleculeShape } from "@/lib/stores/createSingleMoleculeVisualizationStore";
import { ExportBoxOverlay } from "@/components/export-box-overlay";
import { WebGLErrorOverlay } from "@/components/webgl-error-overlay";
import { SpatialScaleBar } from "@/components/spatial-scale-bar";


// Point cloud key helpers for assigned/unassigned separation
const ASSIGNED_KEY_SUFFIX = ":assigned";
const UNASSIGNED_KEY_SUFFIX = ":unassigned";

// Blue background for in-progress gene-loading toasts, so they read as
// activity rather than a neutral notification. Closing one cancels the load.
const LOADING_TOAST_STYLE = { background: "#1d4ed8", color: "#ffffff" };

function formatMB(bytes: number): string {
  return (bytes / 1024 ** 2).toFixed(1);
}

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
  const [webglFailed, setWebglFailed] = useState(false);
  const sceneRef = useRef<THREE.Scene | null>(null);
  const rendererRef = useRef<THREE.WebGLRenderer | null>(null);
  const cameraRef = useRef<THREE.Camera | null>(null);
  const controlsRef = useRef<any>(null);
  const pointCloudsRef = useRef<Map<string, THREE.Points>>(new Map());
  const lastDatasetIdRef = useRef<string | null>(null);
  const lastViewModeRef = useRef<string | null>(null);
  const baselineCameraDistanceRef = useRef<number | null>(null);
  const animationFrameIdRef = useRef<number | null>(null);
  const selectedGenesRef = useRef<Map<string, any>>(new Map());
  const globalScaleRef = useRef<number>(1);
  const hasAutoFittedRef = useRef<boolean>(false);
  // Gene streams in flight, keyed by point-cloud key. They deliberately
  // OUTLIVE the effect run that started them: selecting another gene, or
  // toggling one on/off, re-runs the effect, and a 30-second download must
  // not restart (nor block the genes processed after it).
  const activeStreamsRef = useRef<Map<string, AbortController>>(new Map());
  // Whole-file gene downloads in flight, keyed by gene. Unlike streams these
  // belong to a single effect run: cleanup aborts them, and closing a gene's
  // loading toast aborts + deselects it.
  const wholeLoadControllersRef = useRef<Map<string, AbortController>>(
    new Map(),
  );
  // Scene rebuilds (view-mode switch) dispose every cloud, so streams filling
  // them have to stop. Null until the first point-cloud effect run.
  const streamViewModeRef = useRef<string | null>(null);
  const sceneGroupRef = useRef<THREE.Group | null>(null);  // outer: rotation/flip
  const innerGroupRef = useRef<THREE.Group | null>(null);  // inner: offset by -center
  const dataCenterRef = useRef<THREE.Vector3>(new THREE.Vector3());

  // Get dataset from store - using stable selector to prevent re-renders
  const currentDatasetId = usePanelSingleMoleculeStore(
    (state) => state.currentDatasetId,
  );
  const dataset = usePanelSingleMoleculeStore((state) =>
    currentDatasetId ? state.datasets.get(currentDatasetId) : null,
  );

  // Get visualization settings from store
  const {
    selectedGenes, globalScale, viewMode, showAssigned, showUnassigned,
    sceneRotation, flipX, flipY,
    exportBoxEnabled, exportBoxWidthMm, exportBoxHeightMm,
    exportBoxCenterPx, setExportBoxCenterPx, removeGene,
  } = usePanelSingleMoleculeVisualizationStore();

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
    // Reset auto-fit flag so camera refits for new dataset/viewMode
    hasAutoFittedRef.current = false;

    // Clear any existing canvas before creating new one
    if (containerRef.current) {
      while (containerRef.current.firstChild) {
        containerRef.current.removeChild(containerRef.current.firstChild);
      }
    }

    // Initialize Three.js scene with viewMode preference (not dataset dimensions)
    // WebGL context creation can fail (GPU blocklist, stale browser update,
    // software GL) — surface that instead of leaving a silently blank canvas.
    let init: ReturnType<typeof initializeScene> | null = null;

    try {
      init = initializeScene(containerRef.current, {
        is2D: viewMode === "2D",
      });
      setWebglFailed(false);
    } catch (e) {
      console.error("[SingleMoleculeThreeScene] WebGL init failed:", e);
      setWebglFailed(true);
    }
    if (!init) return;
    const { scene, camera, renderer, controls } = init;

    sceneRef.current = scene;
    rendererRef.current = renderer;
    cameraRef.current = camera;
    controlsRef.current = controls;

    // Nested groups for scene transforms:
    // outerGroup (position=center, rotation/flip) → innerGroup (position=-center, holds point clouds)
    const outerGroup = new THREE.Group();
    const innerGroup = new THREE.Group();
    outerGroup.add(innerGroup);
    scene.add(outerGroup);
    sceneGroupRef.current = outerGroup;
    innerGroupRef.current = innerGroup;

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

      // Update point sizes via the shader's dotSize uniform. The shader does
      // perspective division itself (dotSize * proj11 / -z), so zoom scaling
      // happens automatically without per-frame work beyond uniform writes.
      const s0 = VISUALIZATION_CONFIG.SINGLE_MOLECULE_POINT_BASE_SIZE;

      pointCloudsRef.current.forEach((pointCloud, key) => {
        const gene = geneFromKey(key);
        const geneViz = selectedGenesRef.current.get(gene);

        if (geneViz) {
          const isUnassigned = key.endsWith(UNASSIGNED_KEY_SUFFIX);
          const scale = isUnassigned
            ? geneViz.unassignedLocalScale
            : geneViz.localScale;

          updateDotSize(
            pointCloud,
            scale * globalScaleRef.current * s0 * SM_DOT_SIZE_FACTOR,
          );
        }
      });

      renderer.render(scene, camera);
    };

    // Start custom animation loop
    customAnimate();
    markSceneReady();

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
    const activeStreams = activeStreamsRef.current;

    /** Stop a gene's stream (deselected, hidden, or its scene went away). */
    const abortStream = (key: string) => {
      const controller = activeStreams.get(key);

      if (controller) {
        controller.abort();
        activeStreams.delete(key);
      }
    };

    /**
     * Fit the camera to the current point-cloud bounds, once. Called as soon
     * as the FIRST data exists — the first whole-loaded gene, or the first
     * streamed batch — so the scene shows immediately instead of staying
     * black until every default gene finishes. Guards against empty/degenerate
     * bounds (all clouds still streaming with zero points): fitting an empty
     * THREE.Box3 produces a NaN camera that never recovers, which is why
     * streamed default genes used to never show up.
     */
    const tryAutoFitCamera = () => {
      if (
        hasAutoFittedRef.current ||
        pointCloudsRef.current.size === 0 ||
        !cameraRef.current ||
        !controlsRef.current
      ) {
        return;
      }
      const box = new THREE.Box3();

      pointCloudsRef.current.forEach((pc) => {
        box.expandByObject(pc);
      });
      if (box.isEmpty()) return;

      const center = new THREE.Vector3();
      const size = new THREE.Vector3();

      box.getCenter(center);
      box.getSize(size);

      const maxDim = Math.max(size.x, size.y, size.z);

      if (!Number.isFinite(maxDim) || maxDim <= 0) return;

      hasAutoFittedRef.current = true;
      const camera = cameraRef.current as THREE.PerspectiveCamera;
      const fov = camera.fov * (Math.PI / 180);
      const cameraDistance = ((maxDim / 2) / Math.tan(fov / 2)) * 1.5;

      camera.position.set(center.x, center.y, center.z + cameraDistance);
      camera.lookAt(center);

      // Set group pivot for rotation around data center
      dataCenterRef.current.copy(center);
      if (sceneGroupRef.current && innerGroupRef.current) {
        sceneGroupRef.current.position.copy(center);
        innerGroupRef.current.position.set(-center.x, -center.y, -center.z);
      }
      camera.near = cameraDistance * 0.001;
      camera.far = cameraDistance * 10;
      camera.updateProjectionMatrix();

      // Update controls target to center of data
      const controls = controlsRef.current;

      if (controls.target) {
        controls.target.copy(center);
      }
      controls.update();

      // Update baseline camera distance for zoom-based point sizing
      baselineCameraDistanceRef.current = cameraDistance;

      console.log(
        `[SingleMoleculeThreeScene] Camera auto-fitted: center=(${center.x.toFixed(1)}, ${center.y.toFixed(1)}, ${center.z.toFixed(1)}), distance=${cameraDistance.toFixed(1)}`,
      );
    };

    // A view-mode switch rebuilds the scene, so every cloud being streamed
    // into is gone.
    if (streamViewModeRef.current !== viewMode) {
      const rebuilt = streamViewModeRef.current !== null;

      streamViewModeRef.current = viewMode;
      if (rebuilt) {
        for (const key of Array.from(activeStreams.keys())) abortStream(key);
      }
    }

    // Remove point clouds for unselected genes
    for (const [key, pointCloud] of currentPointClouds.entries()) {
      const gene = geneFromKey(key);

      if (!selectedGenes.has(gene)) {
        console.log(`Removing point cloud: ${key}`);
        abortStream(key);
        if (innerGroupRef.current) innerGroupRef.current.remove(pointCloud);
        else scene.remove(pointCloud);
        pointCloud.geometry.dispose();
        if (Array.isArray(pointCloud.material)) {
          pointCloud.material.forEach((m) => m.dispose());
        } else {
          pointCloud.material.dispose();
        }
        currentPointClouds.delete(key);
      }
    }

    // Helper: create a point cloud from coordinates (Float32Array passed directly).
    // Uses the shared shader (sub-pixel cull + distance-from-target fade) via
    // createSmPointCloud so SM molecules render the same way as single-cell points.
    const createPointCloud = (
      coords: Float32Array,
      geneViz: { color: string; localScale: number },
      shape: MoleculeShape,
      renderOrder: number = 0,
    ): THREE.Points => {
      const moleculeCount = coords.length / 3;
      const userSize =
        geneViz.localScale *
        globalScale *
        VISUALIZATION_CONFIG.SINGLE_MOLECULE_POINT_BASE_SIZE;
      const pointCloud = createSmPointCloud(
        coords,
        moleculeCount,
        geneViz.color,
        userSize * SM_DOT_SIZE_FACTOR,
        shape,
      );

      pointCloud.renderOrder = renderOrder;

      return pointCloud;
    };

    /**
     * Molecules to stream for this gene, or 0 when it should load whole:
     * streaming pays off past SINGLE_MOLECULE_STREAM_THRESHOLD, needs the
     * count up front (to size the buffer) and a dataset that supports it
     * (S3-backed; a locally parsed dataset already has the data in memory).
     */
    const streamableMolecules = async (
      gene: string,
      kind: "assigned" | "unassigned",
    ): Promise<number> => {
      const stream =
        kind === "assigned"
          ? dataset.streamCoordinatesByGene
          : dataset.streamUnassignedCoordinatesByGene;

      if (typeof stream !== "function" || !dataset.resolveMoleculeCount) {
        return 0;
      }

      const expected = await dataset.resolveMoleculeCount(gene, kind);

      return expected !== null &&
        expected >= VISUALIZATION_CONFIG.SINGLE_MOLECULE_STREAM_THRESHOLD
        ? expected
        : 0;
    };

    /**
     * Start filling a cloud for a large gene in the BACKGROUND and return
     * immediately, so the rest of the genes render while the download runs.
     *
     * The stream owns its lifetime through `activeStreams`: it survives the
     * effect runs caused by selecting or toggling other genes, and stops only
     * when its own gene is deselected/hidden or the scene is rebuilt. That is
     * why nothing here closes over `isCancelled` — that flag belongs to one
     * effect run, and this work deliberately outlives it.
     */
    const startStreamingPointCloud = (
      key: string,
      gene: string,
      kind: "assigned" | "unassigned",
      viz: { color: string; localScale: number },
      shape: MoleculeShape,
      renderOrder: number,
      expected: number,
      toastId: string,
    ): void => {
      const userSize =
        viz.localScale *
        globalScale *
        VISUALIZATION_CONFIG.SINGLE_MOLECULE_POINT_BASE_SIZE;
      const cloud = createStreamingSmPointCloud(
        expected,
        viz.color,
        userSize * SM_DOT_SIZE_FACTOR,
        shape,
      );

      cloud.renderOrder = renderOrder;

      if (!innerGroupRef.current) {
        cloud.geometry.dispose();
        (cloud.material as THREE.ShaderMaterial).dispose();

        return;
      }
      innerGroupRef.current.add(cloud);
      currentPointClouds.set(key, cloud);

      const controller = new AbortController();

      activeStreams.set(key, controller);

      const label = kind === "assigned" ? gene : `${gene} (unassigned)`;
      const stream =
        kind === "assigned"
          ? dataset.streamCoordinatesByGene.bind(dataset)
          : dataset.streamUnassignedCoordinatesByGene.bind(dataset);

      void stream(gene, {
        signal: controller.signal,
        onBatch: ({ coords, loadedMolecules }) => {
          // Batches already in flight can arrive after cancellation; they must
          // not repaint the scene or revive the dismissed toast.
          if (controller.signal.aborted) return;
          appendToSmPointCloud(cloud, coords);
          // First data anywhere in the scene: point the camera at it so the
          // stream is visible as it fills.
          tryAutoFitCamera();
          toast.update(toastId, {
            render: `${label}: ${loadedMolecules.toLocaleString()} / ${expected.toLocaleString()} molecules`,
            isLoading: true,
            autoClose: false,
            style: LOADING_TOAST_STYLE,
          });
        },
      })
        .then((coords) => {
          if (controller.signal.aborted) return;
          const total = coords.length / 3;

          console.log(
            `  ✅ ${label} streamed with ${total.toLocaleString()} molecules`,
          );
          toast.update(toastId, {
            render: `${label}: ${total.toLocaleString()} molecules loaded`,
            type: "success",
            isLoading: false,
            autoClose: 3000,
            position: "bottom-left",
            style: undefined,
          });
          if (kind === "assigned") markGeneRendered(gene, total);
        })
        .catch((error: unknown) => {
          if ((error as Error)?.name === "AbortError") {
            // Deselected mid-stream: drop the partial cloud and its toast.
            removePointCloud(key);
            toast.dismiss(toastId);

            return;
          }
          console.error(`[SM] Streaming ${label} failed:`, error);
          removePointCloud(key);
          toast.update(toastId, {
            render: `Failed to load ${label}`,
            type: "error",
            isLoading: false,
            autoClose: 4000,
            style: undefined,
          });
        })
        .finally(() => {
          if (activeStreams.get(key) === controller) activeStreams.delete(key);
        });
    };

    // Helper: update an existing point cloud's color and size in place.
    const updatePointCloudAppearance = (
      pointCloud: THREE.Points,
      geneViz: { color: string; localScale: number },
    ) => {
      const userSize =
        geneViz.localScale *
        globalScale *
        VISUALIZATION_CONFIG.SINGLE_MOLECULE_POINT_BASE_SIZE;

      updateDotSize(pointCloud, userSize * SM_DOT_SIZE_FACTOR);
      updatePointCloudColor(pointCloud, geneViz.color);
    };

    // Helper: remove and dispose a point cloud by key
    const removePointCloud = (key: string) => {
      const pc = currentPointClouds.get(key);

      if (pc) {
        if (innerGroupRef.current) innerGroupRef.current.remove(pc);
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
        // One abort handle for this gene's whole-file downloads (assigned +
        // unassigned). Registered for the iteration so effect cleanup and the
        // toast's close button can truly cancel the fetch.
        const wholeController = new AbortController();

        wholeLoadControllersRef.current.set(gene, wholeController);

        try {
          // Only show loading toast for new genes. Closing the toast is a
          // real cancel: it aborts the in-flight download (stream or whole
          // file) and deselects the gene so nothing restarts it.
          if (!assignedExists) {
            toast.loading(`Loading ${gene}...`, {
              toastId,
              position: "bottom-left",
              autoClose: false,
              style: LOADING_TOAST_STYLE,
              onClose: (reason) => {
                // Only the close button (reason === true) cancels; autoClose
                // and programmatic dismissals pass no reason.
                if (reason !== true) return;
                const whole = wholeLoadControllersRef.current.get(gene);

                if (
                  !whole &&
                  !activeStreams.has(assignedKey(gene)) &&
                  !activeStreams.has(unassignedKey(gene))
                ) {
                  return; // load finished — nothing to cancel
                }
                abortStream(assignedKey(gene));
                abortStream(unassignedKey(gene));
                whole?.abort();
                wholeLoadControllersRef.current.delete(gene);
                removeGene(gene);
              },
            });
          }

          // A stream started by an earlier effect run is still filling this
          // cloud — never re-fetch the gene, just apply appearance changes.
          const assignedStreaming = activeStreams.has(aKey);

          // Big genes stream (points appear while the file downloads); small
          // ones load whole, which stays simpler and is fast enough.
          const expectedAssigned =
            !assignedStreaming && geneViz.showAssigned && !assignedExists
              ? await streamableMolecules(gene, "assigned")
              : 0;
          let moleculeCount = 0;
          // True once a background stream owns this gene's toast/progress.
          let streamStarted = false;

          if (assignedStreaming) {
            const streamingCloud = currentPointClouds.get(aKey);

            if (!geneViz.showAssigned) {
              abortStream(aKey);
              removePointCloud(aKey);
            } else if (streamingCloud) {
              updatePointCloudAppearance(streamingCloud, geneViz);
              updatePointCloudShape(streamingCloud, geneViz.assignedShape);
            }
          } else if (expectedAssigned) {
            // Fills in the background; the loop moves on to the next gene.
            startStreamingPointCloud(
              aKey,
              gene,
              "assigned",
              geneViz,
              geneViz.assignedShape,
              0,
              expectedAssigned,
              toastId,
            );
            streamStarted = true;
          } else {
            // Get assigned coordinates. Abortable (effect cleanup / toast
            // close) with byte progress in the toast, so a slow download is
            // distinguishable from a hung one.
            const coords = await dataset.getCoordinatesByGene(gene, {
              signal: wholeController.signal,
              onProgress: (loadedBytes, totalBytes) => {
                toast.update(toastId, {
                  render: totalBytes
                    ? `${gene}: ${formatMB(loadedBytes)} / ${formatMB(totalBytes)} MB`
                    : `${gene}: ${formatMB(loadedBytes)} MB`,
                  isLoading: true,
                  autoClose: false,
                  style: LOADING_TOAST_STYLE,
                });
              },
            });

            if (isCancelled) {
              toast.dismiss(toastId);

              return;
            }

            moleculeCount = coords.length / 3;

            console.log(`Creating/updating point cloud for gene: ${gene}`);
            console.log(`  Assigned molecules: ${moleculeCount.toLocaleString()}`);
            console.log(`  Color: ${geneViz.color}`);

            // --- Assigned point cloud ---
            if (geneViz.showAssigned) {
              let pointCloud = currentPointClouds.get(aKey);

            if (!pointCloud) {
              pointCloud = createPointCloud(
                coords,
                geneViz,
                geneViz.assignedShape,
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

                if (!innerGroupRef.current) {
                  console.warn(`[SM] innerGroupRef is null when adding ${aKey}, skipping`);
                  pointCloud.geometry.dispose();
                  (pointCloud.material as THREE.ShaderMaterial).dispose();
                  return;
                }
                innerGroupRef.current.add(pointCloud);
                currentPointClouds.set(aKey, pointCloud);
                // First data in the scene: show it now rather than after the
                // remaining genes finish downloading.
                tryAutoFitCamera();

                console.log(
                  `  ✅ Assigned point cloud created with ${moleculeCount} molecules`,
                );
              } else {
                updatePointCloudAppearance(pointCloud, geneViz);
                updatePointCloudShape(pointCloud, geneViz.assignedShape);
                console.log(`  ✅ Assigned point cloud updated`);
              }
            } else {
              abortStream(aKey);
              removePointCloud(aKey);
            }
          }

          // --- Unassigned point cloud ---
          const shouldShowUnassigned =
            dataset.hasUnassigned && showUnassigned && geneViz.showUnassigned;

            if (!streamStarted && !assignedStreaming) {
              console.log(
                `  ✅ Point cloud created with ${moleculeCount} molecules`,
              );
              markGeneRendered(gene, moleculeCount);
            }
          if (shouldShowUnassigned) {
            const unassignedExists = currentPointClouds.has(uKey);
            const uViz = {
              color: geneViz.unassignedColor,
              localScale: geneViz.unassignedLocalScale,
            };

            const unassignedStreaming = activeStreams.has(uKey);
            const expectedUnassigned =
              !unassignedStreaming && !unassignedExists
                ? await streamableMolecules(gene, "unassigned")
                : 0;

            if (unassignedStreaming) {
              const streamingCloud = currentPointClouds.get(uKey);

              if (streamingCloud) {
                updatePointCloudAppearance(streamingCloud, uViz);
                updatePointCloudShape(streamingCloud, geneViz.unassignedShape);
              }
              streamStarted = true;
            } else if (expectedUnassigned) {
              startStreamingPointCloud(
                uKey,
                gene,
                "unassigned",
                uViz,
                geneViz.unassignedShape,
                -1, // Render behind assigned
                expectedUnassigned,
                toastId,
              );
              streamStarted = true;
            } else if (!unassignedExists) {
              // Load unassigned coordinates
              const uCoords = await dataset.getUnassignedCoordinatesByGene(
                gene,
                { signal: wholeController.signal },
              );

              if (isCancelled) {
                toast.dismiss(toastId);

                return;
              }

              if (uCoords.length > 0) {
                const uPointCloud = createPointCloud(
                  uCoords,
                  uViz,
                  geneViz.unassignedShape,
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

                if (!innerGroupRef.current) {
                  console.warn(`[SM] innerGroupRef is null when adding unassigned ${uKey}, skipping`);
                  uPointCloud.geometry.dispose();
                  (uPointCloud.material as THREE.ShaderMaterial).dispose();
                  return;
                }
                innerGroupRef.current.add(uPointCloud);
                currentPointClouds.set(uKey, uPointCloud);

                const uMoleculeCount = uCoords.length / 3;

                console.log(
                  `  ✅ Unassigned point cloud created with ${uMoleculeCount} molecules`,
                );
              }
            } else {
              // Update existing unassigned cloud color/size/shape in place.
              const uPC = currentPointClouds.get(uKey)!;

              updatePointCloudAppearance(uPC, uViz);
              updatePointCloudShape(uPC, geneViz.unassignedShape);
            }
          } else {
            // Remove unassigned cloud if toggle is off or dataset has none
            abortStream(uKey);
            removePointCloud(uKey);
          }

          // Update toast — a background stream finishes its own.
          if (streamStarted || assignedStreaming) {
            // left to the stream
          } else if (!assignedExists) {
            toast.update(toastId, {
              render: `${gene}: ${moleculeCount.toLocaleString()} molecules loaded`,
              type: "success",
              isLoading: false,
              autoClose: 3000,
              position: "bottom-left",
              style: undefined,
            });
          } else {
            toast.dismiss(toastId);
          }
        } catch (error) {
          // Cancelled (toast closed / effect cleanup): not a failure — drop
          // the toast and move on to the next gene.
          if ((error as Error)?.name === "AbortError") {
            toast.dismiss(toastId);
            continue;
          }
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
              style: undefined,
            });
          }
        } finally {
          if (wholeLoadControllersRef.current.get(gene) === wholeController) {
            wholeLoadControllersRef.current.delete(gene);
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
      // E2E timing hook: all selected genes' point clouds are now built.
      markRenderComplete({
        geneCount: dataset.uniqueGenes.length,
        dimensions: dataset.dimensions,
        dataType: "single_molecule",
        pointClouds: pointCloudsRef.current.size,
      });

      // Fallback fit for whatever data exists once the loop finishes (no-op
      // when the first gene/batch already fitted the camera).
      tryAutoFitCamera();
    });

    // Cleanup: cancel the async operation if effect is cleaned up
    return () => {
      console.log("[SingleMoleculeThreeScene] Cancelling point cloud updates");
      // Only this run's sequential work stops here; background streams keep
      // filling their clouds and are aborted when their gene goes away.
      isCancelled = true;
      // Whole-file downloads belong to this run — abort them for real (the
      // next run re-requests anything still selected).
      for (const controller of wholeLoadControllersRef.current.values()) {
        controller.abort();
      }
      wholeLoadControllersRef.current.clear();
    };
  }, [dataset, selectedGenes, globalScale, viewMode, showAssigned, showUnassigned]);

  // Unmount: nothing can consume a stream any more.
  useEffect(() => {
    const streams = activeStreamsRef.current;

    return () => {
      streams.forEach((controller) => controller.abort());
      streams.clear();
    };
  }, []);

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

  // Apply scene transforms (rotation, flip) — pivot around data center (via group)
  useEffect(() => {
    const group = sceneGroupRef.current;
    if (!group) return;

    group.rotation.z = (sceneRotation * Math.PI) / 180;
    group.scale.x = flipX ? -1 : 1;
    group.scale.y = flipY ? -1 : 1;

    console.log(`[SM] Transform: rot=${sceneRotation}° flipX=${flipX} flipY=${flipY}, outerGroup children=${group.children.length}, innerGroup children=${innerGroupRef.current?.children.length}, scene children=${sceneRef.current?.children.length}`);
  }, [sceneRotation, flipX, flipY]);

  return (
    <>
      <div
        ref={containerRef}
        data-testid="sm-scene-canvas"
        className="absolute inset-0 w-full h-full"
        style={{ margin: 0, padding: 0 }}
      />
      {webglFailed && <WebGLErrorOverlay />}
      <SpatialScaleBar
        cameraRef={cameraRef as React.RefObject<THREE.PerspectiveCamera | null>}
        rendererRef={rendererRef}
        controlsRef={controlsRef}
      />
      <ExportBoxOverlay
        canShow={viewMode === "2D"}
        cameraRef={cameraRef as React.RefObject<THREE.PerspectiveCamera | null>}
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
