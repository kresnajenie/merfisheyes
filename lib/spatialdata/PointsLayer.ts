import type { FeatureTiles, NetworkStats, PointsElement } from "./points";

import * as THREE from "three";

import {
  createSmPointCloud,
  updateDotSize,
  updatePointCloudColor,
  updatePointCloudShape,
} from "@/lib/webgl/point-cloud";

/**
 * How much to draw: a point budget per view, shared by the selected genes.
 * "all" is every molecule in view, up to what a browser can hold — about
 * 32 bytes of GPU memory a molecule, so 64 M is some 2 GB.
 */
export type PointsDetail = "auto" | "more" | "all";

export const DETAIL_BUDGET: Record<PointsDetail, number> = {
  auto: 4_000_000,
  more: 16_000_000,
  all: 64_000_000,
};

/** What the layer is drawing, gene by gene, and what it is fetching. */
export interface PointsStatus {
  genes: {
    gene: string;
    /** Level of the current view. */
    level: number;
    levels: number;
    /** Share of the gene's molecules in that level: 1 is every molecule. */
    fraction: number;
    /** Molecules on screen (the current view's, and any still standing in). */
    drawn: number;
    /** The gene's molecules in level 0. */
    total: number;
    /** Molecules the current view needs, and how many of them have arrived. */
    wanted: number;
    loaded: number;
    /** Every molecule of this gene in the current view, at level 0. */
    allInView: number;
  }[];
  detail: PointsDetail;
  network: NetworkStats;
}

export interface GeneStyle {
  color: string;
  /** Dot size as createSmPointCloud takes it. */
  dotSize: number;
  shape: "circle" | "square";
  /** THREE renderOrder of the gene's points. */
  renderOrder: number;
}

// Points kept in memory for tiles that have scrolled out of view.
const CACHE_POINTS = 24_000_000;
const MAX_CONCURRENT_LOADS = 6;
const REPLAN_MS = 120;
// Rows read at a time within a piece, each drawn as it arrives: one chunk.
const STREAM_ROWS = 65_536;
// Zoomed out, molecules would shrink below a pixel and vanish; the layer
// samples each view anyway, so they're kept at least this big on screen.
const MIN_PX = 1;
// Coarser molecules standing in for a view still loading are drawn larger,
// so the area reads as filled rather than sparse, up to this factor.
const STAND_IN_GROWTH = 3;

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
  /** Rows read so far; the piece is drawn up to here. */
  loaded: number;
  loading: boolean;
  lastWanted: number;
  /** Dot size factor: 1, or larger while standing in. */
  growth: number;
}

interface GenePlan {
  level: number;
  /** The gene's molecules in each level. */
  totals: number[];
  /** The current view's pieces, drawn as they stream in. */
  wanted: Piece[];
  /** Pieces of earlier views kept up until the current one has arrived. */
  standIns: Set<Piece>;
  allInView: number;
}

const done = (p: Piece) => p.loaded >= p.count;

/**
 * Draws the molecules of a multiscale points element (docs/points-zarr-spec.md)
 * for the genes selected, reading only what the camera looks at.
 *
 * Per gene it picks the finest level whose molecules in view fit its share
 * of a point budget — every molecule when zoomed in, a uniform sample when
 * zoomed out — and reads those tiles' rows a chunk at a time, drawing each
 * chunk as it arrives. The previous view stays up, its coarser molecules
 * drawn larger, until the new one has fully arrived; since a coarser level
 * is a subset of a finer one, the two overlap without doubling.
 */
export class PointsLayer {
  visible = true;

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
  private detail: PointsDetail = "auto";
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

        if (plan) {
          for (const p of [...plan.wanted, ...plan.standIns]) this.hide(p);
        }
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
        this.genes.set(gene, style);
        this.restyle(gene);
      }
    }
    this.genes = new Map(genes);
    this.dirty = true;
  }

  /** How much to draw per view (see DETAIL_BUDGET). */
  setDetail(detail: PointsDetail) {
    if (detail === this.detail) return;
    this.detail = detail;
    this.dirty = true;
  }

  dispose() {
    this.disposed = true;
    cancelAnimationFrame(this.rafId);
    this.root.removeFromParent();
    for (const p of this.pieces.values()) this.drop(p);
    this.pieces.clear();
  }

  /** A snapshot of what is drawn and fetched, for a status display. */
  status(): PointsStatus {
    const genes: PointsStatus["genes"] = [];

    for (const gene of this.genes.keys()) {
      const plan = this.plans.get(gene);

      if (!plan) continue;
      let drawn = 0;
      let wanted = 0;
      let loaded = 0;

      for (const p of plan.wanted) {
        wanted += p.count;
        loaded += Math.min(p.loaded, p.count);
      }
      for (const p of [...plan.wanted, ...plan.standIns]) {
        if (p.points?.visible) drawn += p.loaded;
      }
      genes.push({
        gene,
        level: plan.level,
        levels: plan.totals.length,
        fraction: plan.totals[0] ? plan.totals[plan.level] / plan.totals[0] : 1,
        drawn,
        total: plan.totals[0] ?? 0,
        wanted,
        loaded,
        allInView: plan.allInView,
      });
    }

    return { genes, detail: this.detail, network: this.element.network() };
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
    for (const p of this.pieces.values()) this.applySize(p);
  }

  private applySize(p: Piece) {
    const style = this.genes.get(this.element.features[p.feature]);

    if (p.points && style) {
      updateDotSize(p.points, Math.max(style.dotSize, this.minDot) * p.growth);
    }
  }

  /** Choose each gene's level and tiles for this view, and queue what's missing. */
  private async plan(rect: Rect | null) {
    const views: {
      gene: string;
      f: number;
      runs: TileRuns[];
      counts: number[];
    }[] = [];

    for (const gene of this.genes.keys()) {
      const f = this.featureIndex.get(gene);

      if (f === undefined) continue;
      const perLevel = await Promise.all(
        this.element.levels.map((_, l) => this.element.featureTiles(l, f)),
      );

      if (this.disposed) return;
      const runs = perLevel.map((ft, l) => this.runs(f, l, ft));

      views.push({
        gene,
        f,
        runs,
        counts: runs.map((r, l) => this.countInView(r, l, rect)),
      });
    }
    // The budget is shared by the selected genes. Smallest first: each takes
    // the finest level that fits an equal share of what is left, so what a
    // small gene doesn't need goes to the bigger ones.
    views.sort((a, b) => a.counts[0] - b.counts[0]);
    let budget = DETAIL_BUDGET[this.detail];

    for (const [i, { gene, f, runs, counts }] of views.entries()) {
      if (!this.genes.has(gene)) continue;
      const share = budget / (views.length - i);
      let level = runs.length - 1;

      for (let l = 0; l < runs.length; l++) {
        if (counts[l] <= share) {
          level = l;
          break;
        }
      }
      budget -= counts[level];
      const wanted = this.piecesInView(f, level, runs[level], rect);
      const plan = this.plans.get(gene) ?? {
        level,
        totals: [],
        wanted: [],
        standIns: new Set<Piece>(),
        allInView: 0,
      };
      const next = new Set(wanted);

      // What is on screen and not part of the new view stands in for it
      for (const p of [...plan.wanted, ...plan.standIns]) {
        if (next.has(p)) continue;
        if (p.points?.visible) plan.standIns.add(p);
        else plan.standIns.delete(p);
      }
      for (const p of wanted) {
        plan.standIns.delete(p);
        p.lastWanted = this.frame;
        if (!p.points && !p.loading) {
          p.loading = true;
          this.queue.push(p);
        }
      }
      plan.level = level;
      plan.wanted = wanted;
      plan.totals = runs.map((r) => r.cum[r.cum.length - 1]);
      plan.allInView = counts[0];
      this.plans.set(gene, plan);
      this.refresh(gene);
    }
    this.pump();
    this.evict();
  }

  /**
   * Show the current view's pieces as far as they have arrived. Until all
   * have, earlier ones stay up, coarser ones drawn larger; then they go.
   */
  private refresh(gene: string) {
    const plan = this.plans.get(gene);

    if (!plan) return;
    for (const p of plan.wanted) {
      if (p.points) p.points.visible = p.loaded > 0;
      this.setGrowth(p, 1);
    }
    if (plan.wanted.every(done)) {
      for (const p of plan.standIns) this.hide(p);
      plan.standIns.clear();

      return;
    }
    for (const p of plan.standIns) {
      this.setGrowth(
        p,
        p.level > plan.level
          ? Math.min(STAND_IN_GROWTH, 2 ** (p.level - plan.level))
          : 1,
      );
    }
  }

  private setGrowth(p: Piece, growth: number) {
    if (p.growth === growth) return;
    p.growth = growth;
    this.applySize(p);
  }

  private hide(p: Piece) {
    if (p.points) p.points.visible = false;
    this.setGrowth(p, 1);
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
        loaded: 0,
        loading: false,
        lastWanted: 0,
        growth: 1,
      };
      this.pieces.set(key, p);
    }

    return p;
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

  /** Read a piece a chunk at a time, drawing each chunk as it arrives. */
  private async load(p: Piece) {
    const gene = this.element.features[p.feature];

    try {
      const style = this.genes.get(gene) ?? {
        color: "#ffffff",
        dotSize: 1,
        shape: "circle" as const,
        renderOrder: 0,
      };
      const positions = new Float32Array(p.count * 3);
      const pc = createSmPointCloud(
        positions,
        p.count,
        style.color,
        Math.max(style.dotSize, this.minDot),
        style.shape,
      );
      const attr = pc.geometry.getAttribute(
        "position",
      ) as THREE.BufferAttribute;

      pc.visible = false;
      pc.renderOrder = style.renderOrder;
      // Bounds aren't known until every chunk is in
      pc.frustumCulled = false;
      pc.geometry.setDrawRange(0, 0);
      this.root.add(pc);
      p.points = pc;

      const end = p.start + p.count;

      for (let lo = p.start; lo < end; ) {
        const hi = Math.min(
          end,
          (Math.floor(lo / STREAM_ROWS) + 1) * STREAM_ROWS,
        );
        const xyz = await this.element.readXyz(p.level, lo, hi);

        if (this.disposed || !p.points) return;
        positions.set(xyz, (lo - p.start) * 3);
        attr.addUpdateRange((lo - p.start) * 3, xyz.length);
        attr.needsUpdate = true;
        p.loaded = hi - p.start;
        pc.geometry.setDrawRange(0, p.loaded);
        this.refresh(gene);
        lo = hi;
      }
    } catch (e) {
      console.warn(`Failed to load molecules ${p.key}:`, e);
    } finally {
      p.loading = false;
    }
  }

  private restyle(gene: string) {
    const f = this.featureIndex.get(gene);
    const style = this.genes.get(gene)!;

    for (const p of this.pieces.values()) {
      if (p.feature !== f || !p.points) continue;
      updatePointCloudColor(p.points, style.color);
      updatePointCloudShape(p.points, style.shape);
      this.applySize(p);
      p.points.renderOrder = style.renderOrder;
    }
  }

  /** Free the least recently wanted pieces that aren't on screen. */
  private evict() {
    const inUse = new Set<Piece>();

    for (const plan of this.plans.values()) {
      plan.wanted.forEach((p) => inUse.add(p));
      plan.standIns.forEach((p) => inUse.add(p));
    }
    let held = 0;
    const idle: Piece[] = [];

    for (const p of this.pieces.values()) {
      if (!p.points) continue;
      held += p.count;
      if (!inUse.has(p) && !p.loading) idle.push(p);
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
