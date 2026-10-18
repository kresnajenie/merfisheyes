"use client";

import type { MoleculeLayer } from "@/lib/webgl/MoleculeLayer";
import type { LayerContext, OpenedElements } from "@/lib/scene/registry";

import { useEffect, useMemo, useRef, useState } from "react";
import * as THREE from "three";
import { Spinner } from "@heroui/react";

import { initializeScene } from "@/lib/webgl/scene-manager";
import { LayerHost } from "@/lib/scene/LayerHost";
import { moleculeElements } from "@/lib/scene/registry";
import { useLayers } from "@/lib/scene/useLayers";
import { LayerBindings, LayerControls } from "@/components/layer-controls";
import {
  usePanelSingleMoleculeStore,
  usePanelSingleMoleculeVisualizationStore,
} from "@/lib/hooks/usePanelStores";
import {
  markSceneReady,
  markRenderComplete,
  markGeneRendered,
} from "@/lib/utils/test-hooks";
import { ExportBoxOverlay } from "@/components/export-box-overlay";
import { WebGLErrorOverlay } from "@/components/webgl-error-overlay";
import { SpatialScaleBar } from "@/components/spatial-scale-bar";

/**
 * A single-molecule dataset on its own: the scene, camera and controls, with
 * the genes drawn by a `MoleculeLayer` — the same layer the cell scene uses
 * for its molecule overlay, so loading, streaming and cancelling behave the
 * same everywhere.
 */
export function SingleMoleculeThreeScene({
  elements,
}: {
  /** What is drawn under the molecules, opened by the viewer (an image). */
  elements?: OpenedElements;
} = {}) {
  const containerRef = useRef<HTMLDivElement>(null);
  const [webglFailed, setWebglFailed] = useState(false);
  const sceneRef = useRef<THREE.Scene | null>(null);
  const rendererRef = useRef<THREE.WebGLRenderer | null>(null);
  const cameraRef = useRef<THREE.Camera | null>(null);
  const controlsRef = useRef<any>(null);
  const lastDatasetIdRef = useRef<string | null>(null);
  const lastViewModeRef = useRef<string | null>(null);
  const baselineCameraDistanceRef = useRef<number | null>(null);
  const hasAutoFittedRef = useRef<boolean>(false);
  const invalidateRef = useRef<() => void>(() => {});
  const sceneGroupRef = useRef<THREE.Group | null>(null); // outer: rotation/flip
  const innerGroupRef = useRef<THREE.Group | null>(null); // inner: offset by -center
  const dataCenterRef = useRef<THREE.Vector3>(new THREE.Vector3());

  // Get dataset from store - using stable selector to prevent re-renders
  const currentDatasetId = usePanelSingleMoleculeStore(
    (state) => state.currentDatasetId,
  );
  const dataset = usePanelSingleMoleculeStore((state) =>
    currentDatasetId ? state.datasets.get(currentDatasetId) : null,
  );

  const {
    viewMode,
    sceneRotation,
    flipX,
    flipY,
    exportBoxEnabled,
    exportBoxWidthMm,
    exportBoxHeightMm,
    exportBoxCenterPx,
    setExportBoxCenterPx,
    removeGene,
  } = usePanelSingleMoleculeVisualizationStore();

  /**
   * Fit the camera to the molecules drawn so far, once. Called as soon as
   * the FIRST data exists — the first whole-loaded gene, or the first
   * streamed batch — so the scene shows immediately instead of staying
   * black until every default gene finishes. Guards against empty bounds
   * (all clouds still streaming with zero points): fitting an empty
   * THREE.Box3 produces a NaN camera that never recovers.
   */
  const tryAutoFitCamera = (moleculeLayer: MoleculeLayer) => {
    if (
      hasAutoFittedRef.current ||
      !cameraRef.current ||
      !controlsRef.current
    ) {
      return;
    }
    const box = moleculeLayer.bounds;

    if (!box) return;

    const center = new THREE.Vector3();
    const size = new THREE.Vector3();

    box.getCenter(center);
    box.getSize(size);

    const maxDim = Math.max(size.x, size.y, size.z);

    if (!Number.isFinite(maxDim) || maxDim <= 0) return;

    hasAutoFittedRef.current = true;
    const camera = cameraRef.current as THREE.PerspectiveCamera;
    const fov = camera.fov * (Math.PI / 180);
    const cameraDistance = (maxDim / 2 / Math.tan(fov / 2)) * 1.5;

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

  // The molecules and whatever the viewer opened under them, as layers of
  // the contract. They outlive scene rebuilds (a view-mode switch);
  // `onSceneReady` attaches them to each new scene.
  const handedImage = elements?.image;
  const layerElements = useMemo<OpenedElements>(
    () => ({ image: handedImage, ...moleculeElements(dataset ?? null) }),
    [handedImage, dataset],
  );
  const layerCtx: LayerContext = {
    molecules: {
      // Fit the camera as soon as the first data exists
      onFirstData: tryAutoFitCamera,
      // Closing a gene's loading toast cancels its download; take it out of
      // the selection so nothing restarts it
      onCancelGene: removeGene,
      onGeneRendered: markGeneRendered,
      // Whole loads finish before this; streams keep filling in the background
      onSynced: (moleculeLayer) => {
        if (!dataset) return;
        // E2E timing hook: all selected genes' point clouds are now built.
        markRenderComplete({
          geneCount: dataset.uniqueGenes.length,
          dimensions: dataset.dimensions,
          dataType: "single_molecule",
          pointClouds: moleculeLayer.pointClouds.size,
        });
        // Fallback fit for whatever data exists once the loop finishes (no-op
        // when the first gene/batch already fitted the camera).
        tryAutoFitCamera(moleculeLayer);
      },
    },
  };
  const { onSceneReady, layers } = useLayers(layerElements, layerCtx);

  // Effect 1: Scene initialization
  useEffect(() => {
    if (!containerRef.current || !dataset) return;

    // Only re-initialize if dataset ID AND viewMode are unchanged AND scene exists
    // This prevents unnecessary re-init during React Strict Mode double-render
    if (
      lastDatasetIdRef.current === dataset.id &&
      lastViewModeRef.current === viewMode &&
      sceneRef.current
    ) {
      return;
    }

    console.log(
      `[SingleMoleculeThreeScene] Initializing scene for ${dataset.name} (${viewMode})`,
    );

    // Update last dataset ID and viewMode
    lastDatasetIdRef.current = dataset.id;
    lastViewModeRef.current = viewMode;
    // Reset auto-fit flag so camera refits for new dataset/viewMode
    hasAutoFittedRef.current = false;

    // Clear any existing canvas before creating new one
    while (containerRef.current.firstChild) {
      containerRef.current.removeChild(containerRef.current.firstChild);
    }

    // Initialize Three.js scene with viewMode preference (not dataset dimensions)
    // WebGL context creation can fail (GPU blocklist, stale browser update,
    // software GL) — surface that instead of leaving a silently blank canvas.
    let init: ReturnType<typeof initializeScene> | null = null;
    let host: LayerHost | null = null;

    try {
      init = initializeScene(containerRef.current, {
        is2D: viewMode === "2D",
        onFrame: (now) => host?.tick(now),
        // On demand: the layer and the controls ask for frames; 10 fps is
        // the safety net for anything that forgets to
        idleFps: 10,
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

    // The layers join through the host, which ticks them from the scene
    // loop. The genes: point sizes follow the camera inside the shader
    // (dotSize * proj11 / -z), so nothing has to be updated per frame.
    host = new LayerHost(innerGroup, camera, renderer, init.invalidate);
    invalidateRef.current = init.invalidate;

    init.animate();
    markSceneReady();
    const removeLayers = onSceneReady(host.handle);
    const { dispose } = init;

    // Cleanup
    return () => {
      removeLayers();
      host?.dispose();

      // Force lose WebGL context before disposing
      const gl = renderer.getContext();

      if (gl) {
        const loseContextExt = gl.getExtension("WEBGL_lose_context");

        if (loseContextExt) loseContextExt.loseContext();
      }
      // The loop, the controls, the renderer, the resize observer, the canvas
      dispose();

      rendererRef.current = null;
      sceneRef.current = null;
      // DON'T reset lastDatasetIdRef here - keep it so we can skip re-init on React Strict Mode double-render
    };
    // Re-initialize scene when viewMode changes to get the right controls
    // (TrackballControls for 3D, OrbitControls for 2D).
  }, [dataset, viewMode]);

  // Apply scene transforms (rotation, flip) — pivot around data center (via group)
  useEffect(() => {
    const group = sceneGroupRef.current;

    if (!group) return;

    group.rotation.z = (sceneRotation * Math.PI) / 180;
    group.scale.x = flipX ? -1 : 1;
    group.scale.y = flipY ? -1 : 1;
    invalidateRef.current();
  }, [sceneRotation, flipX, flipY]);

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
        data-testid="sm-scene-canvas"
        style={{ margin: 0, padding: 0 }}
      />
      {webglFailed && <WebGLErrorOverlay />}
      <LayerControls slot="panel" />
      <LayerBindings ctx={layerCtx} />
      <SpatialScaleBar
        cameraRef={cameraRef as React.RefObject<THREE.PerspectiveCamera | null>}
        controlsRef={controlsRef}
        rendererRef={rendererRef}
      />
      <ExportBoxOverlay
        cameraRef={cameraRef as React.RefObject<THREE.PerspectiveCamera | null>}
        canShow={viewMode === "2D"}
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
