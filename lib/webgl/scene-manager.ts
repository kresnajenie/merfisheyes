import * as THREE from "three";
import { OrbitControls } from "three/examples/jsm/controls/OrbitControls.js";
import { TrackballControls } from "three/examples/jsm/controls/TrackballControls.js";

import { hooksEnabled, markFrame, markMemory } from "@/lib/utils/test-hooks";

export interface SceneSetup {
  scene: THREE.Scene;
  camera: THREE.PerspectiveCamera;
  renderer: THREE.WebGLRenderer;
  controls: OrbitControls | TrackballControls;
  animate: () => void;
  dispose: () => void;
  /**
   * Request one more frame. Only meaningful with `idleFps`; call it after
   * mutating anything the scene draws (uniforms, textures, visibility).
   */
  invalidate: () => void;
}

export interface SceneOptions {
  is2D?: boolean;
  cameraPosition?: THREE.Vector3;
  lookAtPosition?: THREE.Vector3;
  near?: number;
  far?: number;
  /**
   * Called every animation frame before the scene is rendered — where a
   * scene ticks its layers (lib/scene/LayerHost.ts).
   */
  onFrame?: (now: number) => void;
  /**
   * Redraw on demand instead of every frame, falling back to this rate when
   * nothing has changed.
   *
   * The default (undefined) keeps the historical behaviour: redraw forever, so
   * a stationary multi-million-point cloud spends the whole GPU budget
   * reproducing an identical image. With this set, a frame is drawn when the
   * controls report movement or something calls `invalidate()`; the slow tick
   * is a safety net for any mutation that forgets to.
   */
  idleFps?: number;
}

/**
 * Initializes the Three.js scene, camera, renderer, and controls
 */
export function initializeScene(
  container: HTMLElement,
  options: SceneOptions = {},
): SceneSetup {
  const {
    is2D = false,
    cameraPosition,
    lookAtPosition,
    near,
    far,
    idleFps,
    onFrame,
  } = options;

  // Scene setup
  const scene = new THREE.Scene();

  scene.background = new THREE.Color(0x000000);

  // Camera setup
  const camera = new THREE.PerspectiveCamera(
    75,
    container.clientWidth / container.clientHeight,
    near ?? 0.1,
    far ?? 10000,
  );

  // Set camera position
  if (cameraPosition) {
    camera.position.copy(cameraPosition);
  } else {
    camera.position.set(0, 0, 500);
  }

  // Set camera look-at
  if (lookAtPosition) {
    camera.lookAt(lookAtPosition);
  }

  // Renderer setup. preserveDrawingBuffer keeps the WebGL framebuffer
  // around after compositing so canvas.toBlob() / toDataURL() can read
  // pixels (used by the export-box-overlay screenshot tool).
  const renderer = new THREE.WebGLRenderer({
    antialias: true,
    preserveDrawingBuffer: true,
    powerPreference: "high-performance",
  });

  renderer.setSize(container.clientWidth, container.clientHeight);
  // Capped: a 3× display would cost nine fragments per CSS pixel for no
  // visible gain on point clouds; the shaders size points in CSS pixels
  renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
  container.appendChild(renderer.domElement);

  // Controls setup based on dimensionality
  let controls: OrbitControls | TrackballControls;

  if (is2D) {
    // 2D: Use OrbitControls with rotation disabled and left-click panning
    const orbitControls = new OrbitControls(camera, renderer.domElement);

    orbitControls.enableRotate = false;
    orbitControls.mouseButtons = {
      LEFT: THREE.MOUSE.PAN,
      MIDDLE: THREE.MOUSE.DOLLY,
      RIGHT: THREE.MOUSE.ROTATE,
    };
    orbitControls.enableDamping = true;
    orbitControls.dampingFactor = 0.25;
    orbitControls.enableZoom = true;
    orbitControls.enablePan = true;

    if (lookAtPosition) {
      orbitControls.target.copy(lookAtPosition);
      orbitControls.update();
    }

    controls = orbitControls;
  } else {
    // 3D: Use TrackballControls
    const trackballControls = new TrackballControls(
      camera,
      renderer.domElement,
    );

    trackballControls.rotateSpeed = 1.5;
    trackballControls.zoomSpeed = 1.2;
    trackballControls.panSpeed = 0.8;
    trackballControls.dynamicDampingFactor = 0.15;

    if (lookAtPosition) {
      trackballControls.target.copy(lookAtPosition);
      trackballControls.update();
    }

    controls = trackballControls;
  }

  // Animation loop
  let animationId: number;
  let dirty = true;
  let lastDraw = 0;
  const invalidate = () => {
    dirty = true;
  };

  // Damping and inertia move the camera without any new input, so the controls
  // themselves are what tell us a frame is still worth drawing.
  controls.addEventListener("change", invalidate);

  // Frame instrumentation for the E2E budgets (docs/specs/BUDGETS.md): the
  // interval between drawn frames and the CPU time inside render. A no-op
  // outside NEXT_PUBLIC_E2E=1 — the check is one boolean per frame.
  const instrumented = hooksEnabled();
  let lastTickAt = 0;
  let lastMemoryAt = 0;
  // Every loop tick reports its interval (main-thread jank); a tick that
  // draws reports its render time too.
  const tick = (now: number, drew: boolean, renderMs: number) => {
    if (!instrumented) return;
    markFrame(lastTickAt ? now - lastTickAt : 0, drew ? renderMs : null);
    lastTickAt = now;
  };
  const draw = () => {
    if (!instrumented) {
      renderer.render(scene, camera);

      return;
    }
    const t0 = performance.now();

    renderer.render(scene, camera);
    const t1 = performance.now();

    tick(t0, true, t1 - t0);
    if (t1 - lastMemoryAt > 1000) {
      lastMemoryAt = t1;
      markMemory({
        geometries: renderer.info.memory.geometries,
        textures: renderer.info.memory.textures,
        drawCalls: renderer.info.render.calls,
        points: renderer.info.render.points,
        triangles: renderer.info.render.triangles,
      });
    }
  };

  const animate = () => {
    animationId = requestAnimationFrame(animate);
    controls.update();
    const now = performance.now();

    onFrame?.(now);

    if (idleFps === undefined) {
      draw();

      return;
    }

    if (!dirty && now - lastDraw < 1000 / idleFps) {
      tick(now, false, 0);

      return;
    }

    dirty = false;
    lastDraw = now;
    draw();
  };

  // Handle container resize via ResizeObserver (supports split screen panels)
  const handleResize = () => {
    const width = container.clientWidth;
    const height = container.clientHeight;

    if (width === 0 || height === 0) return;
    camera.aspect = width / height;
    camera.updateProjectionMatrix();
    renderer.setSize(width, height);
    invalidate();
  };

  const resizeObserver = new ResizeObserver(handleResize);

  resizeObserver.observe(container);

  // Cleanup function
  const dispose = () => {
    controls.removeEventListener("change", invalidate);
    resizeObserver.disconnect();
    cancelAnimationFrame(animationId);
    controls.dispose();
    renderer.dispose();
    if (container.contains(renderer.domElement)) {
      container.removeChild(renderer.domElement);
    }
  };

  return {
    scene,
    camera,
    renderer,
    controls,
    animate,
    dispose,
    invalidate,
  };
}
