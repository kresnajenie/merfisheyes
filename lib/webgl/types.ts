/**
 * TypeScript types and interfaces for WebGL point cloud
 */

export interface PointData {
  x: number;
  y: number;
  z: number;
  r: number;
  g: number;
  b: number;
  size?: number;
  alpha?: number;
}

export interface PointCloudConfig {
  count: number;
  xRange: [number, number];
  yRange: [number, number];
  zRange: [number, number];
  defaultSize?: number;
  defaultAlpha?: number;
  dotSize?: number;
}
