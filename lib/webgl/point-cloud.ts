import type { PointData, PointCloudConfig } from "./types";

import * as THREE from "three";

import {
  fragmentShader,
  moleculeFragmentShader,
  moleculeVertexShader,
  vertexShader,
} from "./shaders";

/** The renderer's pixel ratio, as scene-manager caps it. */
export function pixelRatio(): number {
  return typeof window === "undefined"
    ? 1
    : Math.min(window.devicePixelRatio || 1, 2);
}

/** The distance-filter uniforms every point material shares. */
function filterUniforms() {
  return {
    uPixelRatio: { value: pixelRatio() },
    uTargetCenter: { value: new THREE.Vector3(0, 0, 0) },
    uTargetRadius: { value: 1.0 },
    uTargetFeather: { value: 0.1 },
    uTargetFilterEnabled: { value: 0.0 },
  };
}

/**
 * The material of a molecule cloud: one colour, size, alpha and shape for
 * the whole cloud, as uniforms, so restyling it uploads nothing. Opaque:
 * the depth test rejects molecules hidden behind closer ones.
 */
export function createMoleculeMaterial(
  colorHex: string,
  dotSize: number,
  shape: "circle" | "square" = "circle",
): THREE.ShaderMaterial {
  return new THREE.ShaderMaterial({
    uniforms: {
      ...filterUniforms(),
      dotSize: { value: dotSize },
      uColor: { value: new THREE.Color(colorHex) },
      uAlpha: { value: 1.0 },
      uShape: { value: shape === "square" ? 1.0 : 0.0 },
    },
    vertexShader: moleculeVertexShader,
    fragmentShader: moleculeFragmentShader,
    transparent: false,
  });
}

/**
 * A molecule cloud's opacity, as its `uAlpha` uniform. Below 1 the cloud
 * blends and stops writing depth, so overlapping molecules build up where
 * they are dense instead of the nearest one hiding the rest; at 1 it is
 * opaque again. Switching mode recompiles once; moving the slider within a
 * mode uploads nothing.
 */
export function setMoleculeOpacity(
  material: THREE.ShaderMaterial,
  opacity: number,
): void {
  material.uniforms.uAlpha.value = opacity;
  const blended = opacity < 1;

  if (material.transparent !== blended) {
    material.transparent = blended;
    material.depthWrite = !blended;
    material.needsUpdate = true;
  }
}

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
      ...filterUniforms(),
      dotSize: { value: dotSize },
      uShape: { value: 0.0 }, // circle
      uWorldSized: { value: 0.0 },
      uHalfHeight: { value: 400.0 },
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
 * Creates an opaque, single-color point cloud: positions only (12 bytes a
 * molecule), colour and shape as uniforms of its material, with the shared
 * sub-pixel cull, distance-from-target fade and sizing math.
 */
export function createSmPointCloud(
  positions: Float32Array,
  count: number,
  colorHex: string,
  dotSize: number,
  shape: "circle" | "square" = "circle",
): THREE.Points {
  const geometry = new THREE.BufferGeometry();

  geometry.setAttribute(
    "position",
    new THREE.BufferAttribute(positions.subarray(0, count * 3), 3),
  );

  return new THREE.Points(
    geometry,
    createMoleculeMaterial(colorHex, dotSize, shape),
  );
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

  geometry.setDrawRange(0, 0);
  geometry.boundingBox = new THREE.Box3();
  geometry.boundingSphere = new THREE.Sphere();

  const points = new THREE.Points(
    geometry,
    createMoleculeMaterial(colorHex, dotSize, shape),
  );

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
  const positionAttr = geometry.getAttribute(
    "position",
  ) as THREE.BufferAttribute;
  const positions = positionAttr.array as Float32Array;
  const offsetMolecules = (points.userData.streamedMolecules as number) ?? 0;
  const capacity =
    (points.userData.streamCapacity as number) ?? positions.length / 3;

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
 * Recolours a molecule cloud: one uniform, nothing uploaded. A cloud with
 * a per-vertex colour attribute (the cells) is repainted instead.
 */
export function updatePointCloudColor(
  pointsMesh: THREE.Points,
  colorHex: string,
): void {
  const material = pointsMesh.material as THREE.ShaderMaterial;

  if (material.uniforms?.uColor) {
    (material.uniforms.uColor.value as THREE.Color).set(colorHex);

    return;
  }
  const colorAttr = pointsMesh.geometry.getAttribute(
    "color",
  ) as THREE.BufferAttribute;

  if (!colorAttr) return;
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
      ...filterUniforms(),
      dotSize: { value: dotSize },
      uShape: { value: 0.0 }, // circle
      uWorldSized: { value: 0.0 },
      uHalfHeight: { value: 400.0 },
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
/**
 * Draw a cell cloud at the cells' real size and shape (a Visium spot's
 * diameter, a bin's side, in world units), or back to sized dots (null).
 * The dot size is then set by `updateDotSize` in world units.
 */
export function setCellFootprint(
  pointsMesh: THREE.Points,
  shape: "circle" | "square" | null,
): void {
  const uniforms = (pointsMesh.material as THREE.ShaderMaterial).uniforms;

  uniforms.uWorldSized.value = shape ? 1.0 : 0.0;
  uniforms.uShape.value = shape === "square" ? 1.0 : 0.0;
}

/** The canvas's half height in CSS pixels, for a real-size cloud. */
export function updateHalfHeight(pointsMesh: THREE.Points, px: number): void {
  const uniforms = (pointsMesh.material as THREE.ShaderMaterial).uniforms;

  if (uniforms.uHalfHeight && uniforms.uHalfHeight.value !== px) {
    uniforms.uHalfHeight.value = px;
  }
}

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
 * Updates point cloud buffer attributes (colors, sizes, alphas) without
 * recreating geometry. The new values are copied into the attribute's own
 * array and re-uploaded in place, so a tick allocates no GPU buffer and
 * leaves none for the garbage collector; only a size change replaces one.
 */
export function updatePointCloudAttributes(
  pointsMesh: THREE.Points,
  colors?: Float32Array,
  sizes?: Float32Array,
  alphas?: Float32Array,
): void {
  const geometry = pointsMesh.geometry;
  const put = (name: string, values: Float32Array, itemSize: number) => {
    const attr = geometry.getAttribute(name) as
      | THREE.BufferAttribute
      | undefined;

    if (attr && attr.array.length === values.length) {
      (attr.array as Float32Array).set(values);
      attr.needsUpdate = true;
    } else {
      geometry.setAttribute(name, new THREE.BufferAttribute(values, itemSize));
    }
  };

  if (colors) put("color", colors, 3);
  if (sizes) put("size", sizes, 1);
  if (alphas) put("alpha", alphas, 1);
}
