/**
 * The contract between a scene and the layers drawn in it.
 *
 * A layer draws one element kind (an image, a points element, a shapes
 * element, legacy molecule files) inside the group a scene gives it. Every
 * layer has the same lifecycle — attach, tick, report, abort, dispose — so
 * a scene, a controls panel, the URL state, the budgets and the tests treat
 * them alike. How a layer is told WHAT to draw stays typed per layer (an
 * image's channels, a points layer's genes): see docs/specs/LAYER-CONTRACT.md.
 */
import type * as THREE from "three";

export type ElementKind = "image" | "points" | "shapes" | "molecules";

/**
 * The part of the data plane (z = 0 in the scene group's frame) on screen,
 * and how many scene units one screen pixel covers there. Computed once per
 * frame by the scene for every layer.
 */
export interface ViewRect {
  minX: number;
  minY: number;
  maxX: number;
  maxY: number;
  unitsPerPixel: number;
}

export interface FrameInfo {
  /** Monotonic frame counter of the scene. */
  frame: number;
  /** performance.now() at the start of the frame. */
  now: number;
  /** Null when the camera does not look onto the data plane. */
  view: ViewRect | null;
}

/** What a scene hands a layer, and what a viewer gets from `onSceneReady`. */
export interface SceneHandle {
  /** Holds the points; a layer's root goes in here and shares its transforms. */
  group: THREE.Group;
  camera: THREE.PerspectiveCamera;
  renderer: THREE.WebGLRenderer;
  /** The view this frame; cached, so every layer sees the same rectangle. */
  view(): ViewRect | null;
  /** Ask for a render. Nothing redraws otherwise once the scene is still. */
  invalidate(): void;
  /** Attach a layer to this scene; the returned function detaches it. */
  attach(layer: Layer): () => void;
  /**
   * The cell dots' colour attribute (rgb per cell), as currently drawn.
   * Absent in a scene without cells.
   */
  cellColors?: () => THREE.BufferAttribute | null;
}

export type LayerPhase = "idle" | "loading" | "ready" | "error";

/** What every layer can say about itself, in the same words. */
export interface LayerStatus {
  phase: LayerPhase;
  /** 0..1 of the current view's wanted content that has arrived; null when unknown. */
  progress: number | null;
  /** On screen now, in the layer's unit (molecules, tiles, polygons). */
  drawn: number;
  /** What the current view asks for, in the same unit. */
  wanted: number;
  /** Kept in memory for the layer, in the same unit. */
  resident: number;
  /** Network reads the layer counts itself, when it can. */
  bytes: number;
  requests: number;
  message?: string;
  error?: string;
}

export interface LayerBudget {
  /** Bytes one view may read before it is considered settled. */
  bytesPerView: number;
  /** Bytes the layer may keep resident (GPU and CPU). */
  residentBytes: number;
}

export interface Layer {
  readonly kind: ElementKind;
  /**
   * Join a scene. May be called again after a scene is rebuilt (a view-mode
   * switch); the layer re-creates what it draws from what it kept.
   */
  attach(scene: SceneHandle): void;
  /** Leave the scene, keeping what was loaded for a later attach. */
  detach(): void;
  /** Called by the scene once per frame, before it renders. */
  tick(frame: FrameInfo): void;
  status(): LayerStatus;
  /** Push: called whenever `status()` would answer differently. Never polled. */
  subscribe(listener: () => void): () => void;
  setBudget(budget: LayerBudget): void;
  /** Cancel every download in flight; what is drawn stays. */
  abort(): void;
  /** Free everything. The layer is unusable afterwards. */
  dispose(): void;
}
