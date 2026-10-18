import type { FrameInfo, Layer, SceneHandle, ViewRect } from "./contract";

import * as THREE from "three";

import { viewRect } from "./view-rect";

/**
 * The scene's side of the layer contract: owns the attached layers, hands
 * each the same frame (one view rectangle, computed once), and is the one
 * place a layer joins or leaves a scene. A scene creates one after its
 * renderer, group and camera exist, calls `tick()` from its loop before
 * rendering, and `dispose()` when it is torn down.
 */
export class LayerHost {
  readonly handle: SceneHandle;

  private layers = new Set<Layer>();
  private frame = 0;
  private viewFrame = -1;
  private cachedView: ViewRect | null = null;

  constructor(
    group: THREE.Group,
    camera: THREE.PerspectiveCamera,
    renderer: THREE.WebGLRenderer,
    invalidate: () => void,
    cellColors?: () => THREE.BufferAttribute | null,
  ) {
    this.handle = {
      group,
      camera,
      renderer,
      invalidate,
      view: () => this.view(),
      attach: (layer) => this.attach(layer),
      cellColors,
    };
  }

  /** Layers attached, in attach order. */
  get attached(): ReadonlySet<Layer> {
    return this.layers;
  }

  attach(layer: Layer): () => void {
    // Always through the layer, so one that was disposed refuses
    layer.attach(this.handle);
    this.layers.add(layer);
    this.handle.invalidate();

    return () => this.detach(layer);
  }

  detach(layer: Layer): void {
    if (!this.layers.delete(layer)) return;
    layer.detach();
    this.handle.invalidate();
  }

  /** One frame: the same view for every layer, then each layer's tick. */
  tick(now = performance.now()): void {
    this.frame++;
    const info: FrameInfo = { frame: this.frame, now, view: this.view() };

    for (const layer of this.layers) layer.tick(info);
  }

  /** Detach every layer (they keep their data; the scene is going away). */
  dispose(): void {
    for (const layer of [...this.layers]) this.detach(layer);
  }

  private view(): ViewRect | null {
    if (this.viewFrame !== this.frame) {
      this.viewFrame = this.frame;
      this.cachedView = viewRect(
        this.handle.group,
        this.handle.camera,
        this.handle.renderer,
      );
    }

    return this.cachedView;
  }
}
