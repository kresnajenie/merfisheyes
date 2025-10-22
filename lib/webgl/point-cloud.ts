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
 * Creates a Three.js Points mesh with custom shader material from point data
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
    const colorAttribute = geometry.getAttribute("color");

    if (colorAttribute) {
      geometry.setAttribute("color", new THREE.BufferAttribute(colors, 3));
    } else {
      geometry.setAttribute("color", new THREE.BufferAttribute(colors, 3));
    }
  }

  if (sizes) {
    const sizeAttribute = geometry.getAttribute("size");

    if (sizeAttribute) {
      geometry.setAttribute("size", new THREE.BufferAttribute(sizes, 1));
    } else {
      geometry.setAttribute("size", new THREE.BufferAttribute(sizes, 1));
    }
  }

  if (alphas) {
    const alphaAttribute = geometry.getAttribute("alpha");

    if (alphaAttribute) {
      geometry.setAttribute("alpha", new THREE.BufferAttribute(alphas, 1));
    } else {
      geometry.setAttribute("alpha", new THREE.BufferAttribute(alphas, 1));
    }
  }
}
