import type { CellPolygons } from "./shapes";
import type {
  FrameInfo,
  LayerStatus,
  SceneHandle,
  ViewRect,
} from "@/lib/scene/contract";

import * as THREE from "three";

import { BaseLayer } from "@/lib/scene/BaseLayer";
import {
  markLayerChanged,
  markLayerContent,
  markLayerSettled,
} from "@/lib/utils/test-hooks";

export type ShapesMode = "outline" | "filled";
export type ShapesStatus = "idle" | "loading" | "ready" | "error";

// Polygons are drawn once a view is this narrow (µm), when cells can be told
// apart; wider, they would be a blur of outlines.
export const SHAPES_SHOW_WIDTH = 1000;
// Bucket size of the spatial index over polygons (µm).
const GRID = 100;
// Polygons drawn at once, at most.
const MAX_SHAPES = 60_000;
const REBUILD_MS = 100;

const ids = new WeakMap<object, number>();
let nextId = 0;
const idOf = (o: object) => {
  if (!ids.has(o)) ids.set(o, nextId++);

  return ids.get(o)!;
};

/**
 * Draws the cells' polygons (a SpatialData shapes element) in the cells'
 * current colours, when zoomed in: as outlines over the cell dots, or filled
 * under them, at a chosen opacity. The polygons are loaded the first time
 * they would be drawn, and only those in view are built.
 */
export class ShapesLayer extends BaseLayer {
  readonly kind = "shapes" as const;
  /** Called when status or settings change, for controls to re-render. */
  onChange?: () => void;

  private _visible = true;
  private _mode: ShapesMode = "outline";
  private _opacity = 0.8;
  private _status: ShapesStatus = "idle";
  private _zoomedIn = false;
  private polygons: CellPolygons | null = null;
  private grid: {
    x0: number;
    y0: number;
    nx: number;
    ny: number;
    cells: Map<number, number[]>;
  } | null = null;
  private triangles = new Map<number, number[]>();
  private outline: THREE.LineSegments | null = null;
  private fill: THREE.Mesh | null = null;
  private lineMaterial = new THREE.LineBasicMaterial({
    vertexColors: true,
    transparent: true,
    depthTest: false,
    depthWrite: false,
  });
  private fillMaterial = new THREE.MeshBasicMaterial({
    vertexColors: true,
    transparent: true,
    side: THREE.DoubleSide,
    depthTest: false,
    depthWrite: false,
  });
  private lastKey = "";
  private lastBuild = 0;
  private shown = 0;
  /** The cell dots' colour attribute (rgb per cell, in table order). */
  private cellColors: () => THREE.BufferAttribute | null = () => null;

  constructor(private load: () => Promise<CellPolygons>) {
    super();
    this.applyOpacity();
  }

  protected onAttach(scene: SceneHandle) {
    this.cellColors = scene.cellColors ?? (() => null);
    this.lastKey = "";
  }

  status(): LayerStatus {
    return {
      phase: this._status,
      progress: this._status === "ready" ? 1 : null,
      drawn: this.root.visible ? this.shown : 0,
      wanted: this.shown,
      resident: this.polygons ? this.polygons.starts.length - 1 : 0,
      bytes: 0,
      requests: 0,
      message:
        this._status === "ready" && !this._zoomedIn
          ? `Zoom in to under ${SHAPES_SHOW_WIDTH / 1000} mm across to see them.`
          : undefined,
    };
  }

  get visible() {
    return this._visible;
  }
  set visible(v: boolean) {
    this._visible = v;
    this.settingsChanged();
  }
  get mode() {
    return this._mode;
  }
  set mode(m: ShapesMode) {
    this._mode = m;
    this.settingsChanged();
  }
  get opacity() {
    return this._opacity;
  }
  set opacity(o: number) {
    this._opacity = o;
    this.applyOpacity();
    this.settingsChanged();
  }
  /** idle until first drawn, then loading, ready or error. */
  get shapesStatus() {
    return this._status;
  }
  /** Whether the view is narrow enough for polygons to be drawn. */
  get zoomedIn() {
    return this._zoomedIn;
  }

  dispose() {
    super.dispose();
    this.outline?.geometry.dispose();
    this.fill?.geometry.dispose();
    this.lineMaterial.dispose();
    this.fillMaterial.dispose();
  }

  private settingsChanged() {
    this.lastKey = "";
    this.onChange?.();
    this.changed();
  }

  private applyOpacity() {
    this.lineMaterial.opacity = this._opacity;
    this.fillMaterial.opacity = this._opacity;
  }

  tick(frame: FrameInfo) {
    if (this.disposed || !this.scene) return;
    const rect = frame.view;
    const zoomedIn = !!rect && rect.maxX - rect.minX <= SHAPES_SHOW_WIDTH;

    if (zoomedIn !== this._zoomedIn) {
      this._zoomedIn = zoomedIn;
      this.onChange?.();
      this.changed();
    }
    this.root.visible = this._visible && zoomedIn && !!this.polygons;
    if (!this._visible || !zoomedIn || !rect) return;
    if (!this.polygons) {
      if (this._status === "idle") this.fetch();

      return;
    }
    const colors = this.cellColors();
    const key = [
      Math.round(rect.minX / 25),
      Math.round(rect.minY / 25),
      Math.round(rect.maxX / 25),
      Math.round(rect.maxY / 25),
      this._mode,
      colors ? `${idOf(colors)}:${colors.version}` : "-",
    ].join(",");
    const now = performance.now();

    if (key !== this.lastKey && now - this.lastBuild > REBUILD_MS) {
      this.lastKey = key;
      this.lastBuild = now;
      markLayerChanged("shapes");
      this.build(rect, colors);
      // A build is synchronous: the view is settled as soon as it returns
      markLayerContent("shapes", { resident: this.polygons.starts.length - 1 });
      markLayerSettled("shapes");
      this.changed();
    }
  }

  private fetch() {
    this._status = "loading";
    this.onChange?.();
    this.changed();
    this.load()
      .then((p) => {
        if (this.disposed) return;
        this.polygons = p;
        this.index(p);
        this._status = "ready";
        this.settingsChanged();
      })
      .catch((e) => {
        console.warn("Failed to load cell shapes:", e);
        this._status = "error";
        this.onChange?.();
        this.changed();
      });
  }

  /** Bucket polygons by their first vertex into a GRID-sized index. */
  private index({ coords, starts }: CellPolygons) {
    let x0 = Infinity;
    let y0 = Infinity;
    let x1 = -Infinity;
    let y1 = -Infinity;

    for (let i = 0; i < coords.length; i += 2) {
      x0 = Math.min(x0, coords[i]);
      x1 = Math.max(x1, coords[i]);
      y0 = Math.min(y0, coords[i + 1]);
      y1 = Math.max(y1, coords[i + 1]);
    }
    const nx = Math.max(1, Math.ceil((x1 - x0) / GRID) + 1);
    const ny = Math.max(1, Math.ceil((y1 - y0) / GRID) + 1);
    const cells = new Map<number, number[]>();

    for (let p = 0; p + 1 < starts.length; p++) {
      const v = starts[p] * 2;
      const key =
        Math.floor((coords[v + 1] - y0) / GRID) * nx +
        Math.floor((coords[v] - x0) / GRID);
      let list = cells.get(key);

      if (!list) cells.set(key, (list = []));
      list.push(p);
    }
    this.grid = { x0, y0, nx, ny, cells };
  }

  /** Polygons of the cells in and just around the view. */
  private inView(rect: Rect): number[] {
    const g = this.grid!;
    const margin = (rect.maxX - rect.minX) * 0.15 + GRID;
    const gx0 = Math.max(0, Math.floor((rect.minX - margin - g.x0) / GRID));
    const gx1 = Math.min(
      g.nx - 1,
      Math.floor((rect.maxX + margin - g.x0) / GRID),
    );
    const gy0 = Math.max(0, Math.floor((rect.minY - margin - g.y0) / GRID));
    const gy1 = Math.min(
      g.ny - 1,
      Math.floor((rect.maxY + margin - g.y0) / GRID),
    );
    const out: number[] = [];

    for (let gy = gy0; gy <= gy1; gy++) {
      for (let gx = gx0; gx <= gx1; gx++) {
        const list = g.cells.get(gy * g.nx + gx);

        if (!list) continue;
        for (const p of list) {
          out.push(p);
          if (out.length >= MAX_SHAPES) return out;
        }
      }
    }

    return out;
  }

  private build(rect: Rect, colors: THREE.BufferAttribute | null) {
    const { coords, starts, cell } = this.polygons!;
    const shown = this.inView(rect);

    this.shown = shown.length;
    const rgb = colors?.array as ArrayLike<number> | undefined;
    const colorOf = (p: number, out: number[]) => {
      const c = cell[p];

      if (rgb && c >= 0 && c * 3 + 2 < rgb.length) {
        out[0] = rgb[c * 3];
        out[1] = rgb[c * 3 + 1];
        out[2] = rgb[c * 3 + 2];
      } else {
        out[0] = out[1] = out[2] = 0.6;
      }
    };
    const col = [0, 0, 0];

    if (this._mode === "outline") {
      let edges = 0;

      for (const p of shown) edges += starts[p + 1] - starts[p];
      const pos = new Float32Array(edges * 6);
      const clr = new Float32Array(edges * 6);
      let k = 0;

      for (const p of shown) {
        const a = starts[p];
        const n = starts[p + 1] - a;

        colorOf(p, col);
        for (let i = 0; i < n; i++) {
          const v0 = (a + i) * 2;
          const v1 = (a + ((i + 1) % n)) * 2;

          pos[k] = coords[v0];
          pos[k + 1] = coords[v0 + 1];
          pos[k + 3] = coords[v1];
          pos[k + 4] = coords[v1 + 1];
          for (let j = 0; j < 6; j += 3) {
            clr[k + j] = col[0];
            clr[k + j + 1] = col[1];
            clr[k + j + 2] = col[2];
          }
          k += 6;
        }
      }
      this.replace(pos, clr, null);
    } else {
      let verts = 0;

      for (const p of shown) verts += starts[p + 1] - starts[p];
      const pos = new Float32Array(verts * 3);
      const clr = new Float32Array(verts * 3);
      const index: number[] = [];
      let base = 0;

      for (const p of shown) {
        const a = starts[p];
        const n = starts[p + 1] - a;

        colorOf(p, col);
        for (let i = 0; i < n; i++) {
          pos[(base + i) * 3] = coords[(a + i) * 2];
          pos[(base + i) * 3 + 1] = coords[(a + i) * 2 + 1];
          clr.set(col, (base + i) * 3);
        }
        for (const t of this.triangulate(p)) index.push(base + t);
        base += n;
      }
      this.replace(pos, clr, index);
    }
  }

  /** A polygon's triangles as local vertex indices, computed once. */
  private triangulate(p: number): number[] {
    let tri = this.triangles.get(p);

    if (!tri) {
      const { coords, starts } = this.polygons!;
      const contour: THREE.Vector2[] = [];

      for (let v = starts[p]; v < starts[p + 1]; v++) {
        contour.push(new THREE.Vector2(coords[v * 2], coords[v * 2 + 1]));
      }
      tri = THREE.ShapeUtils.triangulateShape(contour, []).flat();
      this.triangles.set(p, tri);
    }

    return tri;
  }

  private replace(
    pos: Float32Array,
    clr: Float32Array,
    index: number[] | null,
  ) {
    const geometry = new THREE.BufferGeometry();

    geometry.setAttribute("position", new THREE.BufferAttribute(pos, 3));
    geometry.setAttribute("color", new THREE.BufferAttribute(clr, 3));
    this.outline?.removeFromParent();
    this.outline?.geometry.dispose();
    this.outline = null;
    this.fill?.removeFromParent();
    this.fill?.geometry.dispose();
    this.fill = null;
    if (index) {
      geometry.setIndex(index);
      this.fill = new THREE.Mesh(geometry, this.fillMaterial);
      // Above the image, under the molecules and the cell dots
      this.fill.renderOrder = -999.5;
      this.fill.frustumCulled = false;
      this.root.add(this.fill);
    } else {
      this.outline = new THREE.LineSegments(geometry, this.lineMaterial);
      // Over the cell dots
      this.outline.renderOrder = 0.5;
      this.outline.frustumCulled = false;
      this.root.add(this.outline);
    }
  }
}

type Rect = Pick<ViewRect, "minX" | "minY" | "maxX" | "maxY">;
