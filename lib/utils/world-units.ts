import type * as THREE from "three";

// World units in this app are treated as microns (the spatial scale bar
// reads them directly as μm; viewers that store normalized coords hide the
// bar via `wasNormalized` rather than rescaling).
//
// Compute world-microns per CSS pixel for the given perspective camera and
// orbit-style controls. Same formula as `spatial-scale-bar.tsx`.
export function micronsPerPixel(
  camera: THREE.PerspectiveCamera,
  controls: { target?: THREE.Vector3 } | null | undefined,
  canvasHeightPx: number,
): number {
  if (canvasHeightPx <= 0) return 0;
  const target = controls?.target;
  const cameraDistance = target
    ? camera.position.distanceTo(target)
    : camera.position.length();
  const fovRad = (camera.fov / 2) * (Math.PI / 180);
  return (2 * cameraDistance * Math.tan(fovRad)) / canvasHeightPx;
}
