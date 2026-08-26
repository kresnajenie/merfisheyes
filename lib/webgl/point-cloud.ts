import type { PointData, PointCloudConfig } from "./types";

import * as THREE from "three";

import { vertexShader, fragmentShader } from "./shaders";

/**
 * Generates random point data for visualization
 */
export function generateRandomPoints(config: PointCloudConfig): PointData[] {
  const points: PointData[] = [];
  const {
    count,
    xRange,
    yRange,
    zRange,
    defaultSize = 1.0,
    defaultAlpha = 0.5,
  } = config;

  for (let i = 0; i < count; i++) {
    points.push({
      x: Math.random() * (xRange[1] - xRange[0]) + xRange[0],
      y: Math.random() * (yRange[1] - yRange[0]) + yRange[0],
      z: Math.random() * (zRange[1] - zRange[0]) + zRange[0],
      r: Math.random(),
      g: Math.random(),
      b: Math.random(),
      size: defaultSize,
      alpha: defaultAlpha,
    });
  }

  return points;
}

/**
 * Creates a Three.js Points mesh with custom shader material from typed arrays.
 * Accepts positions directly as Float32Array — no intermediate PointData[] needed.
 */
export function createPointCloudFromBuffers(
  positions: Float32Array,
  count: number,
  dotSize: number = 5,
): THREE.Points {
  const geometry = new THREE.BufferGeometry();

  // Position buffer passed directly — zero copy
  geometry.setAttribute("position", new THREE.BufferAttribute(positions, 3));

  // Initialize color/size/alpha buffers with defaults
  const colors = new Float32Array(count * 3); // all zeros (black)
  const sizes = new Float32Array(count).fill(1.0);
  const alphas = new Float32Array(count); // all zeros

  geometry.setAttribute("color", new THREE.BufferAttribute(colors, 3));
  geometry.setAttribute("size", new THREE.BufferAttribute(sizes, 1));
  geometry.setAttribute("alpha", new THREE.BufferAttribute(alphas, 1));

  const material = new THREE.ShaderMaterial({
    uniforms: {
      dotSize: { value: dotSize },
      uTargetCenter: { value: new THREE.Vector3(0, 0, 0) },
      uTargetRadius: { value: 1.0 },
      uTargetFeather: { value: 0.1 },
      uTargetFilterEnabled: { value: 0.0 },
      uShape: { value: 0.0 }, // circle
    },
    vertexShader,
    fragmentShader,
    transparent: true,
  });

  return new THREE.Points(geometry, material);
}

/**
 * Multiplier that converts the legacy "PointsMaterial.size" SM scale (a
 * world-ish unit picked when SM molecules were rendered with
 * `sizeAttenuation: true`) into the `dotSize` uniform the shared shader
 * expects. PointsMaterial scales by `viewportHeight/2`; the shared shader
 * scales by `projectionMatrix[1][1]` (≈1.303 at 75° FOV). For a ~800px
 * viewport these match with a factor of ~300. Exported so it can be tuned
 * if SM molecules render too small/large on very different viewport sizes.
 */
export const SM_DOT_SIZE_FACTOR = 300;

/**
 * Creates an opaque, single-color point cloud using the shared shader.
 * Built for single-molecule clouds where every point in a cloud shares one
 * color and one shape, but we still want the shared sub-pixel cull,
 * distance-from-target fade, and consistent sizing math.
 *
 * Opaque + alphaTest lets the GPU depth-test reject points hidden behind
 * closer ones — the big perf win vs THREE.PointsMaterial with transparent:
 * true (which forces alpha blending and disables early-z).
 */
export function createSmPointCloud(
  positions: Float32Array,
  count: number,
  colorHex: string,
  dotSize: number,
  shape: "circle" | "square" = "circle",
): THREE.Points {
  const geometry = new THREE.BufferGeometry();

  geometry.setAttribute("position", new THREE.BufferAttribute(positions, 3));

  // Fill the per-vertex color buffer with the single cloud color so the
  // shared shader's `attribute vec3 color` path works without modification.
  const c = new THREE.Color(colorHex);
  const colors = new Float32Array(count * 3);
  for (let i = 0; i < count; i++) {
    colors[i * 3] = c.r;
    colors[i * 3 + 1] = c.g;
    colors[i * 3 + 2] = c.b;
  }
  const sizes = new Float32Array(count).fill(1.0);
  const alphas = new Float32Array(count).fill(1.0);

  geometry.setAttribute("color", new THREE.BufferAttribute(colors, 3));
  geometry.setAttribute("size", new THREE.BufferAttribute(sizes, 1));
  geometry.setAttribute("alpha", new THREE.BufferAttribute(alphas, 1));

  const material = new THREE.ShaderMaterial({
    uniforms: {
      dotSize: { value: dotSize },
      uTargetCenter: { value: new THREE.Vector3(0, 0, 0) },
      uTargetRadius: { value: 1.0 },
      uTargetFeather: { value: 0.1 },
      uTargetFilterEnabled: { value: 0.0 },
      uShape: { value: shape === "square" ? 1.0 : 0.0 },
    },
    vertexShader,
    fragmentShader,
    transparent: false,
    alphaTest: 0.5,
  });

  return new THREE.Points(geometry, material);
}

/**
 * An empty SM point cloud sized for `capacity` molecules, to be filled by
 * `appendToSmPointCloud` as a gene streams in.
 *
 * Same material and attributes as `createSmPointCloud`, but nothing is drawn
 * until molecules arrive: the draw range starts at 0 and grows, so the
 * unfilled tail of the buffer is never rendered. Bounds are maintained from
 * the filled portion only (`geometry.boundingBox` computed from all-zero
 * padding would wreck the camera auto-fit).
 */
export function createStreamingSmPointCloud(
  capacity: number,
  colorHex: string,
  dotSize: number,
  shape: "circle" | "square" = "circle",
): THREE.Points {
  const geometry = new THREE.BufferGeometry();
  const positions = new Float32Array(capacity * 3);
  const positionAttr = new THREE.BufferAttribute(positions, 3);

  positionAttr.setUsage(THREE.DynamicDrawUsage);
  geometry.setAttribute("position", positionAttr);

  const c = new THREE.Color(colorHex);
  const colors = new Float32Array(capacity * 3);

  for (let i = 0; i < capacity; i++) {
    colors[i * 3] = c.r;
    colors[i * 3 + 1] = c.g;
    colors[i * 3 + 2] = c.b;
  }
  geometry.setAttribute("color", new THREE.BufferAttribute(colors, 3));
  geometry.setAttribute(
    "size",
    new THREE.BufferAttribute(new Float32Array(capacity).fill(1.0), 1),
  );
  geometry.setAttribute(
    "alpha",
    new THREE.BufferAttribute(new Float32Array(capacity).fill(1.0), 1),
  );

  geometry.setDrawRange(0, 0);
  geometry.boundingBox = new THREE.Box3();
  geometry.boundingSphere = new THREE.Sphere();

  const material = new THREE.ShaderMaterial({
    uniforms: {
      dotSize: { value: dotSize },
      uTargetCenter: { value: new THREE.Vector3(0, 0, 0) },
      uTargetRadius: { value: 1.0 },
      uTargetFeather: { value: 0.1 },
      uTargetFilterEnabled: { value: 0.0 },
      uShape: { value: shape === "square" ? 1.0 : 0.0 },
    },
    vertexShader,
    fragmentShader,
    transparent: false,
    alphaTest: 0.5,
  });

  const points = new THREE.Points(geometry, material);

  // Only the filled prefix is drawn, so let three skip frustum culling until
  // the bounds are meaningful (they are updated on every append).
  points.userData.streamedMolecules = 0;
  points.userData.streamCapacity = capacity;

  return points;
}

/**
 * Copy one streamed batch into a cloud made by `createStreamingSmPointCloud`
 * and reveal it. Returns the new molecule count.
 *
 * Only the written slice is uploaded to the GPU (`addUpdateRange`), and the
 * bounding box grows over the batch's own min/max — both O(batch), not
 * O(total), so appending stays cheap as the cloud fills.
 */
export function appendToSmPointCloud(
  points: THREE.Points,
  coords: Float32Array,
): number {
  const geometry = points.geometry;
  const positionAttr = geometry.getAttribute("position") as THREE.BufferAttribute;
  const positions = positionAttr.array as Float32Array;
  const offsetMolecules = (points.userData.streamedMolecules as number) ?? 0;
  const capacity = (points.userData.streamCapacity as number) ?? positions.length / 3;

  const incoming = Math.min(coords.length / 3, capacity - offsetMolecules);

  if (incoming <= 0) return offsetMolecules;

  const start = offsetMolecules * 3;
  const length = incoming * 3;

  positions.set(coords.subarray(0, length), start);
  positionAttr.needsUpdate = true;
  positionAttr.addUpdateRange(start, length);

  const total = offsetMolecules + incoming;

  points.userData.streamedMolecules = total;
  geometry.setDrawRange(0, total);

  // Grow the bounds over the new points only.
  const box = geometry.boundingBox ?? new THREE.Box3();

  for (let i = 0; i < length; i += 3) {
    const x = coords[i];
    const y = coords[i + 1];
    const z = coords[i + 2];

    if (x < box.min.x) box.min.x = x;
    if (y < box.min.y) box.min.y = y;
    if (z < box.min.z) box.min.z = z;
    if (x > box.max.x) box.max.x = x;
    if (y > box.max.y) box.max.y = y;
    if (z > box.max.z) box.max.z = z;
  }
  geometry.boundingBox = box;
  if (!box.isEmpty()) {
    box.getBoundingSphere(geometry.boundingSphere ?? new THREE.Sphere());
  }

  return total;
}

/** Molecules rendered so far in a streaming cloud. */
export function streamedMoleculeCount(points: THREE.Points): number {
  return (points.userData.streamedMolecules as number) ?? 0;
}

/**
 * Repaints every vertex of an existing point cloud to a new uniform color.
 * Cheaper than recreating the cloud — just rewrites the color buffer.
 */
export function updatePointCloudColor(
  pointsMesh: THREE.Points,
  colorHex: string,
): void {
  const colorAttr = pointsMesh.geometry.getAttribute(
    "color",
  ) as THREE.BufferAttribute;
  const colors = colorAttr.array as Float32Array;
  const c = new THREE.Color(colorHex);
  const count = colors.length / 3;

  for (let i = 0; i < count; i++) {
    colors[i * 3] = c.r;
    colors[i * 3 + 1] = c.g;
    colors[i * 3 + 2] = c.b;
  }
  colorAttr.needsUpdate = true;
}

/**
 * Flips the shape uniform between circle and square.
 */
export function updatePointCloudShape(
  pointsMesh: THREE.Points,
  shape: "circle" | "square",
): void {
  const material = pointsMesh.material as THREE.ShaderMaterial;

  if (material.uniforms?.uShape) {
    material.uniforms.uShape.value = shape === "square" ? 1.0 : 0.0;
  }
}

/**
 * Creates a Three.js Points mesh with custom shader material from PointData objects.
 * @deprecated Use createPointCloudFromBuffers for better performance.
 */
export function createPointCloud(
  data: PointData[],
  dotSize: number = 5,
): THREE.Points {
  const count = data.length;

  // Create buffer geometry
  const geometry = new THREE.BufferGeometry();

  // Create typed arrays for attributes
  const positions = new Float32Array(count * 3);
  const colors = new Float32Array(count * 3);
  const sizes = new Float32Array(count);
  const alphas = new Float32Array(count);

  // Fill buffers from data
  data.forEach((point, i) => {
    positions[i * 3] = point.x;
    positions[i * 3 + 1] = point.y;
    positions[i * 3 + 2] = point.z;

    colors[i * 3] = point.r;
    colors[i * 3 + 1] = point.g;
    colors[i * 3 + 2] = point.b;

    sizes[i] = point.size ?? 1.0;
    alphas[i] = point.alpha ?? 0.5;
  });

  // Set buffer attributes
  geometry.setAttribute("position", new THREE.BufferAttribute(positions, 3));
  geometry.setAttribute("color", new THREE.BufferAttribute(colors, 3));
  geometry.setAttribute("size", new THREE.BufferAttribute(sizes, 1));
  geometry.setAttribute("alpha", new THREE.BufferAttribute(alphas, 1));

  // Create custom shader material
  const material = new THREE.ShaderMaterial({
    uniforms: {
      dotSize: { value: dotSize },
      uTargetCenter: { value: new THREE.Vector3(0, 0, 0) },
      uTargetRadius: { value: 1.0 },
      uTargetFeather: { value: 0.1 },
      uTargetFilterEnabled: { value: 0.0 },
      uShape: { value: 0.0 }, // circle
    },
    vertexShader,
    fragmentShader,
    transparent: true,
  });

  // Create and return points mesh
  const pointsMesh = new THREE.Points(geometry, material);

  return pointsMesh;
}

/**
 * Updates the dotSize uniform of a point cloud material
 */
export function updateDotSize(
  pointsMesh: THREE.Points,
  newDotSize: number,
): void {
  const material = pointsMesh.material as THREE.ShaderMaterial;

  if (material.uniforms && material.uniforms.dotSize) {
    material.uniforms.dotSize.value = newDotSize;
  }
}

/**
 * Updates point cloud buffer attributes (colors, sizes, alphas) without recreating geometry
 * This is much more performant than recreating the entire point cloud
 */
export function updatePointCloudAttributes(
  pointsMesh: THREE.Points,
  colors?: Float32Array,
  sizes?: Float32Array,
  alphas?: Float32Array,
): void {
  const geometry = pointsMesh.geometry;

  if (colors) {
    geometry.setAttribute("color", new THREE.BufferAttribute(colors, 3));
  }

  if (sizes) {
    geometry.setAttribute("size", new THREE.BufferAttribute(sizes, 1));
  }

  if (alphas) {
    geometry.setAttribute("alpha", new THREE.BufferAttribute(alphas, 1));
  }
}
