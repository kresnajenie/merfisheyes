import type { FeatureTiles, PointsElement } from "./points";

import * as THREE from "three";

import {
  createSmPointCloud,
  updateDotSize,
  updatePointCloudColor,
  updatePointCloudShape,
} from "@/lib/webgl/point-cloud";

export interface GeneStyle {
  color: string;
  /** Dot size as createSmPointCloud takes it. */
  dotSize: number;
  shape: "circle" | "square";
  /** THREE renderOrder of the gene's points. */
  renderOrder: number;
}

// Points drawn at once, shared by the selected genes.
const POINT_BUDGET = 4_000_000;
// Points kept in memory for tiles that have scrolled out of view.
const CACHE_POINTS = 24_000_000;
const MAX_CONCURRENT_LOADS = 6;
const REPLAN_MS = 120;
// Zoomed out, molecules would shrink below a pixel and vanish; the layer
// samples each view anyway, so they're kept at least this big on screen.
const MIN_PX = 1;

/** A feature's tiles at a level, with running counts for range sums. */
interface TileRuns extends FeatureTiles {
  /** cum[i] = points in tiles [0, i). */
  cum: Float64Array;
}

interface Piece {
  key: string;
  feature: number;
  level: number;
  start: number;
  count: number;
  points: THREE.Points | null;
  loading: boolean;
  lastWanted: number;
}

interface GenePlan {
  level: number;
  wanted: Piece[];
  /** What is on screen: the last plan whose pieces had all arrived. */
  shown: Piece[];
}

/**
 * Draws the molecules of a multiscale points element (docs/points-zarr-spec.md)
 * for the genes selected, reading only what the camera looks at.
 *
 * Per gene it picks the finest level whose molecules in view fit its share
 * of a point budget — every molecule when zoomed in, a uniform sample when
 * zoomed out — and reads those tiles' rows. A view stays on screen until the
 * next one has fully arrived, so panning and zooming never blank.
 */
export class PointsLayer {
  visible = true;
  /** Called when the molecules drawn change (count, level). */
  onStatus?: (status: {
    points: number;
    levels: Record<string, number>;
  }) => void;

  private root = new THREE.Group();
  private genes = new Map<string, GeneStyle>();
  private plans = new Map<string, GenePlan>();
  private pieces = new Map<string, Piece>();
  private queue: Piece[] = [];
  private activeLoads = 0;
  private frame = 0;
  private rafId = 0;
  private disposed = false;
  private lastRect = "";
  private lastPlan = 0;
  private planning = false;
  private dirty = true;
  private featureIndex: Map<string, number>;
  private runCache = new Map<string, TileRuns>();
  /** Dot size that is MIN_PX on screen at the current distance. */
  private minDot = 0;

  constructor(
    private element: PointsElement,
    parent: THREE.Object3D,
    private camera: THREE.PerspectiveCamera,
    private renderer: THREE.WebGLRenderer,
  ) {
    this.featureIndex = new Map(element.features.map((f, i) => [f, i]));
    parent.add(this.root);
    this.tick();
  }

  /** The genes to draw and how. Genes not in the map are removed. */
  setGenes(genes: Map<string, GeneStyle>) {
    for (const gene of this.genes.keys()) {
      if (!genes.has(gene)) {
        const plan = this.plans.get(gene);

        plan?.shown.forEach((p) => p.points && (p.points.visible = false));
        this.plans.delete(gene);
      }
    }
    for (const [gene, style] of genes) {
      const old = this.genes.get(gene);

      if (
        old &&
        (old.color !== style.color ||
          old.dotSize !== style.dotSize ||
          old.shape !== style.shape ||
          old.renderOrder !== style.renderOrder)
      ) {
        this.restyle(gene, style);
      }
    }
    this.genes = new Map(genes);
    this.dirty = true;
  }

  dispose() {
    this.disposed = true;
    cancelAnimationFrame(this.rafId);
    this.root.removeFromParent();
    for (const p of this.pieces.values()) this.drop(p);
    this.pieces.clear();
  }

  private tick = () => {
    if (this.disposed) return;
    this.rafId = requestAnimationFrame(this.tick);
    this.frame++;
    this.root.visible = this.visible;
    if (!this.visible || this.genes.size === 0) return;
    this.keepVisible();
    const rect = this.visibleRect();
    const sig = rect
      ? [rect.minX, rect.minY, rect.maxX, rect.maxY]
          .map((v) => Math.round(v))
          .join(",")
      : "all";
    const now = performance.now();

    if (
      (sig !== this.lastRect || this.dirty) &&
      !this.planning &&
      now - this.lastPlan > REPLAN_MS
    ) {
      this.lastRect = sig;
      this.dirty = false;
      this.lastPlan = now;
      this.planning = true;
      this.plan(rect).finally(() => (this.planning = false));
    }
  };

  /** Raise dot sizes to at least MIN_PX on screen, as the camera moves. */
  private keepVisible() {
    this.root.updateWorldMatrix(true, false);
    const local = this.camera.position
      .clone()
      .applyMatrix4(this.root.matrixWorld.clone().invert());
    const proj11 = this.camera.projectionMatrix.elements[5];
    const minDot = (MIN_PX * Math.abs(local.z)) / proj11;

    if (Math.abs(minDot - this.minDot) <= this.minDot * 0.02) return;
    this.minDot = minDot;
    for (const p of this.pieces.values()) {
      const style = this.genes.get(this.element.features[p.feature]);

      if (p.points && style) updateDotSize(p.points, this.dotSize(style));
    }
  }

  private dotSize(style: GeneStyle) {
    return Math.max(style.dotSize, this.minDot);
  }

  /** Choose each gene's level and tiles for this view, and queue what's missing. */
  private async plan(rect: Rect | null) {
    const genes = Array.from(this.genes.keys());
    const share = POINT_BUDGET / Math.max(1, genes.length);
    const levels: Record<string, number> = {};

    for (const gene of genes) {
      const f = this.featureIndex.get(gene);

      if (f === undefined) continue;
      const perLevel = await Promise.all(
        this.element.levels.map((_, l) => this.element.featureTiles(l, f)),
      );

      if (this.disposed || !this.genes.has(gene)) return;
      const runs = perLevel.map((ft, l) => this.runs(f, l, ft));
      // Finest level whose points in view fit this gene's share
      let level = runs.length - 1;

      for (let l = 0; l < runs.length; l++) {
        if (this.countInView(runs[l], l, rect) <= share) {
          level = l;
          break;
        }
      }
      levels[gene] = level;
      const wanted = this.piecesInView(f, level, runs[level], rect);

      for (const p of wanted) {
        p.lastWanted = this.frame;
        if (!p.points && !p.loading) {
          p.loading = true;
          this.queue.push(p);
        }
      }
      const plan = this.plans.get(gene) ?? { level, wanted: [], shown: [] };

      plan.level = level;
      plan.wanted = wanted;
      this.plans.set(gene, plan);
      this.settle(gene);
    }
    this.pump();
    this.evict();
    this.report(levels);
  }

  private runs(feature: number, level: number, ft: FeatureTiles): TileRuns {
    const key = `${feature}/${level}`;
    let r = this.runCache.get(key);

    if (!r) {
      const cum = new Float64Array(ft.counts.length + 1);

      for (let i = 0; i < ft.counts.length; i++)
        cum[i + 1] = cum[i] + ft.counts[i];
      r = { ...ft, cum };
      this.runCache.set(key, r);
    }

    return r;
  }

  /** The tile range a view covers at a level (all tiles, without a view). */
  private tileRange(level: number, rect: Rect | null) {
    const { origin, tileSize, nx, ny } = this.element.levels[level].grid;

    if (!rect) return { tx0: 0, tx1: nx - 1, ty0: 0, ty1: ny - 1, nx };

    return {
      tx0: Math.max(0, Math.floor((rect.minX - origin[0]) / tileSize)),
      tx1: Math.min(nx - 1, Math.floor((rect.maxX - origin[0]) / tileSize)),
      ty0: Math.max(0, Math.floor((rect.minY - origin[1]) / tileSize)),
      ty1: Math.min(ny - 1, Math.floor((rect.maxY - origin[1]) / tileSize)),
      nx,
    };
  }

  /** Entries [i0, i1) of a feature's (ascending) tiles with ids in [lo, hi]. */
  private span(r: TileRuns, lo: number, hi: number): [number, number] {
    return [lowerBound(r.tiles, lo), lowerBound(r.tiles, hi + 1)];
  }

  private countInView(r: TileRuns, level: number, rect: Rect | null): number {
    const { tx0, tx1, ty0, ty1, nx } = this.tileRange(level, rect);
    let n = 0;

    if (tx1 < tx0) return 0;
    for (let ty = ty0; ty <= ty1; ty++) {
      const [i0, i1] = this.span(r, ty * nx + tx0, ty * nx + tx1);

      n += r.cum[i1] - r.cum[i0];
    }

    return n;
  }

  /**
   * One piece per tile row in view: the feature's rows there are contiguous.
   * Run ends are snapped to a few tiles so a small pan reuses what's loaded.
   */
  private piecesInView(
    f: number,
    level: number,
    r: TileRuns,
    rect: Rect | null,
  ): Piece[] {
    const { tx0, tx1, ty0, ty1, nx } = this.tileRange(level, rect);
    const width = tx1 - tx0 + 1;
    const snap =
      2 ** Math.max(0, Math.floor(Math.log2(Math.max(1, width / 2))));
    const sx0 = Math.floor(tx0 / snap) * snap;
    const sx1 = Math.min(nx - 1, (Math.floor(tx1 / snap) + 1) * snap - 1);
    const out: Piece[] = [];

    if (tx1 < tx0) return out;
    for (let ty = Math.max(0, ty0 - 1); ty <= ty1 + 1; ty++) {
      const [i0, i1] = this.span(r, ty * nx + sx0, ty * nx + sx1);
      const count = r.cum[i1] - r.cum[i0];

      if (count > 0) {
        out.push(
          this.piece(f, level, `${ty}/${sx0}-${sx1}`, r.starts[i0], count),
        );
      }
    }

    return out;
  }

  private piece(
    feature: number,
    level: number,
    run: string,
    start: number,
    count: number,
  ) {
    const key = `${feature}/${level}/${run}`;
    let p = this.pieces.get(key);

    if (!p) {
      p = {
        key,
        feature,
        level,
        start,
        count,
        points: null,
        loading: false,
        lastWanted: 0,
      };
      this.pieces.set(key, p);
    }

    return p;
  }

  /** Show a gene's wanted pieces once all have arrived; until then keep the old ones. */
  private settle(gene: string) {
    const plan = this.plans.get(gene);

    if (!plan || !plan.wanted.every((p) => p.points)) return;
    const next = new Set(plan.wanted);

    for (const p of plan.shown)
      if (!next.has(p) && p.points) p.points.visible = false;
    for (const p of plan.wanted) p.points!.visible = true;
    plan.shown = plan.wanted;
  }

  private pump() {
    // Coarse tiles first, so a view fills in evenly
    this.queue.sort((a, b) => b.level - a.level);
    while (this.activeLoads < MAX_CONCURRENT_LOADS && this.queue.length) {
      const p = this.queue.shift()!;

      this.activeLoads++;
      this.load(p).finally(() => {
        this.activeLoads--;
        if (!this.disposed) this.pump();
      });
    }
  }

  private async load(p: Piece) {
    try {
      const xyz = await this.element.readXyz(
        p.level,
        p.start,
        p.start + p.count,
      );

      if (this.disposed) return;
      const gene = this.element.features[p.feature];
      const style = this.genes.get(gene) ?? {
        color: "#ffffff",
        dotSize: 1,
        shape: "circle",
        renderOrder: 0,
      };
      const pc = createSmPointCloud(
        xyz,
        p.count,
        style.color,
        this.dotSize(style),
        style.shape,
      );

      pc.visible = false;
      pc.renderOrder = style.renderOrder;
      pc.frustumCulled = true;
      this.root.add(pc);
      p.points = pc;
      this.settle(gene);
    } catch (e) {
      console.warn(`Failed to load molecules ${p.key}:`, e);
    } finally {
      p.loading = false;
    }
  }

  private restyle(gene: string, style: GeneStyle) {
    const f = this.featureIndex.get(gene);

    for (const p of this.pieces.values()) {
      if (p.feature !== f || !p.points) continue;
      updatePointCloudColor(p.points, style.color);
      updatePointCloudShape(p.points, style.shape);
      updateDotSize(p.points, this.dotSize(style));
      p.points.renderOrder = style.renderOrder;
    }
  }

  /** Free the least recently wanted pieces that aren't on screen. */
  private evict() {
    const inUse = new Set<Piece>();

    for (const plan of this.plans.values()) {
      plan.wanted.forEach((p) => inUse.add(p));
      plan.shown.forEach((p) => inUse.add(p));
    }
    let held = 0;
    const idle: Piece[] = [];

    for (const p of this.pieces.values()) {
      if (!p.points) continue;
      held += p.count;
      if (!inUse.has(p)) idle.push(p);
    }
    idle.sort((a, b) => a.lastWanted - b.lastWanted);
    for (const p of idle) {
      if (held <= CACHE_POINTS) break;
      held -= p.count;
      this.drop(p);
      this.pieces.delete(p.key);
    }
  }

  private drop(p: Piece) {
    if (!p.points) return;
    p.points.removeFromParent();
    p.points.geometry.dispose();
    (p.points.material as THREE.Material).dispose();
    p.points = null;
  }

  private report(levels: Record<string, number>) {
    let points = 0;

    for (const plan of this.plans.values())
      for (const p of plan.shown) points += p.count;
    this.onStatus?.({ points, levels });
  }

  /**
   * The part of the z = 0 plane (in the parent's frame) on screen. Null when
   * the camera doesn't look onto the plane.
   */
  private visibleRect(): Rect | null {
    this.root.updateWorldMatrix(true, false);
    const toLocal = this.root.matrixWorld.clone().invert();
    const origin = this.camera.position.clone().applyMatrix4(toLocal);
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
      const p = new THREE.Vector3(nx, ny, 0.5)
        .unproject(this.camera)
        .applyMatrix4(toLocal);
      const t = -origin.z / (p.z - origin.z);

      if (!(t > 0) || !isFinite(t)) return null;
      const x = origin.x + t * (p.x - origin.x);
      const y = origin.y + t * (p.y - origin.y);

      minX = Math.min(minX, x);
      maxX = Math.max(maxX, x);
      minY = Math.min(minY, y);
      maxY = Math.max(maxY, y);
    }

    return { minX, minY, maxX, maxY };
  }
}

interface Rect {
  minX: number;
  minY: number;
  maxX: number;
  maxY: number;
}

/** First index of a sorted array whose value is >= x. */
function lowerBound(a: Uint32Array, x: number): number {
  let lo = 0;
  let hi = a.length;

  while (lo < hi) {
    const mid = (lo + hi) >>> 1;

    if (a[mid] < x) lo = mid + 1;
    else hi = mid;
  }

  return lo;
}
