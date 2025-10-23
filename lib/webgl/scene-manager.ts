import * as THREE from "three";
import { OrbitControls } from "three/examples/jsm/controls/OrbitControls.js";
import { TrackballControls } from "three/examples/jsm/controls/TrackballControls.js";

export interface SceneSetup {
  scene: THREE.Scene;
  camera: THREE.PerspectiveCamera;
  renderer: THREE.WebGLRenderer;
  controls: OrbitControls | TrackballControls;
  animate: () => void;
  dispose: () => void;
}

export interface SceneOptions {
  is2D?: boolean;
  cameraPosition?: THREE.Vector3;
  lookAtPosition?: THREE.Vector3;
}

/**
 * Initializes the Three.js scene, camera, renderer, and controls
 */
export function initializeScene(
  container: HTMLElement,
  options: SceneOptions = {},
): SceneSetup {
  const { is2D = false, cameraPosition, lookAtPosition } = options;

  // Scene setup
  const scene = new THREE.Scene();

  scene.background = new THREE.Color(0x000000);

  // Camera setup
  const camera = new THREE.PerspectiveCamera(
    75,
    container.clientWidth / container.clientHeight,
    0.1,
    10000,
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

  // Renderer setup
  const renderer = new THREE.WebGLRenderer({ antialias: true });

  renderer.setSize(container.clientWidth, container.clientHeight);
  renderer.setPixelRatio(window.devicePixelRatio);
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
  const animate = () => {
    animationId = requestAnimationFrame(animate);
    controls.update();
    renderer.render(scene, camera);
  };

  // Handle window resize
  const handleResize = () => {
    camera.aspect = container.clientWidth / container.clientHeight;
    camera.updateProjectionMatrix();
    renderer.setSize(container.clientWidth, container.clientHeight);
  };

  window.addEventListener("resize", handleResize);

  // Cleanup function
  const dispose = () => {
    window.removeEventListener("resize", handleResize);
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
  };
}
