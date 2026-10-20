import type {
  ElementKind,
  FrameInfo,
  Layer,
  LayerBudget,
  LayerStatus,
  SceneHandle,
} from "./contract";

import * as THREE from "three";

/**
 * What every layer shares: a root group that joins the scene's group on
 * attach and leaves it on detach, push notifications for its controls, the
 * budget, and the scene handle while attached.
 *
 * Subclasses draw into `this.root`, override `tick`, `status`, `abort` and
 * `dispose` (calling super), and call `this.changed()` whenever what they
 * draw or report has changed: that asks the scene for a frame and tells the
 * subscribers, once per frame at most.
 */
export abstract class BaseLayer implements Layer {
  abstract readonly kind: ElementKind;

  protected readonly root = new THREE.Group();
  protected scene: SceneHandle | null = null;
  protected budget: LayerBudget = {
    bytesPerView: Infinity,
    residentBytes: Infinity,
  };
  protected disposed = false;

  private listeners = new Set<() => void>();
  private notifyScheduled = false;

  attach(scene: SceneHandle): void {
    if (this.disposed) throw new Error(`${this.kind} layer is disposed`);
    if (this.scene) this.detach();
    this.scene = scene;
    scene.group.add(this.root);
    this.onAttach(scene);
    this.changed();
  }

  detach(): void {
    if (!this.scene) return;
    this.onDetach();
    this.root.removeFromParent();
    this.scene = null;
  }

  /** The scene while attached; subclasses that need it assert this. */
  protected get attached(): SceneHandle {
    if (!this.scene) throw new Error(`${this.kind} layer is not attached`);

    return this.scene;
  }

  /** Hook: the layer joined a scene (create what it draws). */
  protected onAttach(_scene: SceneHandle): void {}
  /** Hook: the layer is leaving its scene (drop scene-bound objects). */
  protected onDetach(): void {}

  abstract tick(frame: FrameInfo): void;
  abstract status(): LayerStatus;

  subscribe(listener: () => void): () => void {
    this.listeners.add(listener);

    return () => {
      this.listeners.delete(listener);
    };
  }

  setBudget(budget: LayerBudget): void {
    this.budget = budget;
    this.changed();
  }

  abort(): void {}

  dispose(): void {
    this.detach();
    this.disposed = true;
    this.listeners.clear();
  }

  /**
   * Something drawn or reported changed: redraw, and tell the subscribers
   * (coalesced to one call per task, so a burst of changes costs one render
   * of the controls).
   */
  protected changed(): void {
    this.scene?.invalidate();
    if (this.notifyScheduled || this.listeners.size === 0) return;
    this.notifyScheduled = true;
    queueMicrotask(() => {
      this.notifyScheduled = false;
      for (const listener of this.listeners) listener();
    });
  }
}
