import type { ViewRect } from "./contract";

import * as THREE from "three";

const corner = new THREE.Vector3();
const origin = new THREE.Vector3();
const toLocal = new THREE.Matrix4();

/**
 * The part of the plane z = 0 in `group`'s frame that the camera sees, and
 * the scene units per screen pixel there. Null when a corner ray does not
 * reach the plane (the camera looks away from it).
 *
 * One implementation for every layer; the scene calls it once per frame.
 */
export function viewRect(
  group: THREE.Object3D,
  camera: THREE.Camera,
  renderer: THREE.WebGLRenderer,
): ViewRect | null {
  group.updateWorldMatrix(true, false);
  toLocal.copy(group.matrixWorld).invert();
  origin.copy(camera.position).applyMatrix4(toLocal);
  let minX = Infinity;
  let minY = Infinity;
  let maxX = -Infinity;
  let maxY = -Infinity;

  for (const [nx, ny] of [
    [-1, -1],
    [1, -1],
    [1, 1],
    [-1, 1],
  ]) {
    corner.set(nx, ny, 0.5).unproject(camera).applyMatrix4(toLocal);
    const t = -origin.z / (corner.z - origin.z);

    if (!(t > 0) || !isFinite(t)) return null;
    const x = origin.x + t * (corner.x - origin.x);
    const y = origin.y + t * (corner.y - origin.y);

    minX = Math.min(minX, x);
    maxX = Math.max(maxX, x);
    minY = Math.min(minY, y);
    maxY = Math.max(maxY, y);
  }

  const canvas = renderer.domElement;

  return {
    minX,
    minY,
    maxX,
    maxY,
    unitsPerPixel: Math.min(
      (maxX - minX) / Math.max(1, canvas.width),
      (maxY - minY) / Math.max(1, canvas.height),
    ),
  };
}
