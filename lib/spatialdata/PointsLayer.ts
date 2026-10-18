import type { FeatureTiles, NetworkStats, PointsElement } from "./points";
import type {
  FrameInfo,
  LayerBudget,
  LayerStatus,
  SceneHandle,
  ViewRect,
} from "@/lib/scene/contract";

import * as THREE from "three";

import { BaseLayer } from "@/lib/scene/BaseLayer";
import { hint } from "@/lib/storage/ElementStore";
import {
  hooksEnabled,
  markLayer,
  markLayerChanged,
  markLayerContent,
  markLayerSettled,
} from "@/lib/utils/test-hooks";
import { createSmPointCloud } from "@/lib/webgl/point-cloud";

/**
 * How much to draw: a point budget and a byte budget per view, shared by
 * the selected genes (docs/specs/BUDGETS.md: a view settles within a second
 * at 12 MB/s when it reads at most 12 MB). "all" is every molecule in view,
 * up to what a browser can hold — about 17 bytes of GPU memory a molecule,
 * so 64 M is over 1 GB — and is not bounded in bytes.
 */
export type PointsDetail = "auto" | "more" | "all";

export const DETAIL_BUDGET: Record<PointsDetail, number> = {
  auto: 8_000_000,
  more: 16_000_000,
  all: 64_000_000,
};

export const DETAIL_BYTES: Record<PointsDetail, number> = {
  auto: 12 * 1024 * 1024,
  more: 36 * 1024 * 1024,
  all: Infinity,
};

/** GPU bytes per resident molecule: xyz float32, colour, size, alpha. */
export const BYTES_PER_POINT = 17;

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
  /**
   * Hidden genes (the legend's eye, or cmd+click to show one alone) still
   * count towards the budget and keep what they have loaded, so showing
   * them again is instant and doesn't change anyone's level.
   */
  visible: boolean;
}

// Points kept in memory for tiles that have scrolled out of view, unless
// the scene's budget (setBudget) says less.
const CACHE_POINTS = 40_000_000;
const MAX_CONCURRENT_LOADS = 6;
const REPLAN_MS = 120;
// Rows read at a time within a piece, each drawn as it arrives: one chunk.
const STREAM_ROWS = 65_536;
// Molecules handed to the GPU per frame. Arrived chunks join the scene at
// this pace, so a burst of them never stalls drawing.
const UPLOAD_PER_FRAME = 1_000_000;
// Size and alpha are 1 for every molecule: one buffer, shared by all chunks.
const ONES = new Uint8Array(STREAM_ROWS).fill(255);
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
  /** One small point cloud per chunk read, sharing `material`. */
  points: THREE.Group | null;
  material: THREE.ShaderMaterial | null;
  /** Rows read so far; the piece is drawn up to here. */
  loaded: number;
  loading: boolean;
  /** Cancels the read in flight; an aborted piece is dropped. */
  controller: AbortController | null;
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
export class PointsLayer extends BaseLayer {
  readonly kind = "points" as const;

  private _visible = true;
  private camera!: THREE.PerspectiveCamera;
  private genes = new Map<string, GeneStyle>();
  private plans = new Map<string, GenePlan>();
  private pieces = new Map<string, Piece>();
  private queue: Piece[] = [];
  private activeLoads = 0;
  private frame = 0;
  private lastRect = "";
  private lastPlan = 0;
  private planning = false;
  private dirty = true;
  private detail: PointsDetail = "auto";
  private featureIndex: Map<string, number>;
  private runCache = new Map<string, TileRuns>();
  /** Dot size that is MIN_PX on screen at the current distance. */
  private minDot = 0;
  private template: THREE.ShaderMaterial | null = null;
  /** Chunks read and waiting for their turn to be uploaded. */
  private uploads: { piece: Piece; points: THREE.Points }[] = [];

  constructor(private element: PointsElement) {
    super();
    this.featureIndex = new Map(element.features.map((f, i) => [f, i]));
  }

  get visible() {
    return this._visible;
  }
  set visible(v: boolean) {
    if (v === this._visible) return;
    this._visible = v;
    this.changed();
  }

  protected onAttach(scene: SceneHandle) {
    this.camera = scene.camera;
    this.lastRect = "";
    this.dirty = true;
  }

  /** The contract's view of the layer; `pointsStatus()` has the detail. */
  status(): LayerStatus {
    const ps = this.pointsStatus();
    let drawn = 0;
    let wanted = 0;
    let loaded = 0;

    for (const g of ps.genes) {
      drawn += g.drawn;
      wanted += g.wanted;
      loaded += g.loaded;
    }
    let resident = 0;

    for (const p of this.pieces.values()) if (p.points) resident += p.loaded;

    return {
      phase:
        ps.genes.length === 0 ? "idle" : loaded < wanted ? "loading" : "ready",
      progress: wanted ? Math.min(1, loaded / wanted) : null,
      drawn,
      wanted,
      resident,
      residentBytes: resident * BYTES_PER_POINT,
      bytes: ps.network.bytes,
      requests: ps.network.requests,
    };
  }

  /**
   * Cancel every read, queued or in flight. Pieces that had fully arrived
   * stay drawn; a piece cut short is dropped, so a later view asks again.
   */
  abort() {
    for (const p of this.queue) p.loading = false;
    this.queue = [];
    for (const p of this.pieces.values()) if (p.loading) this.cancel(p);
  }

  /** A budget from the scene: the cache and, for Auto, the bytes per view. */
  setBudget(budget: LayerBudget) {
    super.setBudget(budget);
    this.dirty = true;
    this.changed();
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
    const shown = (g: string) => this.genes.get(g)?.visible;
    const toggled = Array.from(genes.keys()).filter(
      (g) => this.genes.has(g) && shown(g) !== genes.get(g)!.visible,
    );

    this.genes = new Map(genes);
    // Shown or hidden at once, from what is loaded
    for (const gene of toggled) this.refresh(gene);
    this.dirty = true;
    this.changed();
  }

  /** How much to draw per view (see DETAIL_BUDGET). */
  setDetail(detail: PointsDetail) {
    if (detail === this.detail) return;
    this.detail = detail;
    this.dirty = true;
    this.changed();
  }

  dispose() {
    super.dispose();
    for (const p of this.pieces.values()) this.drop(p);
    this.pieces.clear();
    this.template?.dispose();
  }

  /** A snapshot of what is drawn and fetched, gene by gene, for the panel. */
  pointsStatus(): PointsStatus {
    const genes: PointsStatus["genes"] = [];

    for (const [gene, style] of this.genes) {
      const plan = this.plans.get(gene);

      if (!plan || !style.visible) continue;
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

  tick(frame: FrameInfo) {
    if (this.disposed || !this.scene) return;
    this.frame++;
    this.root.visible = this._visible;
    if (!this._visible || this.genes.size === 0) return;
    this.keepVisible();
    this.uploadSome();
    const rect = frame.view;
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
  }

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

    if (p.material && style) {
      p.material.uniforms.dotSize.value =
        Math.max(style.dotSize, this.minDot) * p.growth;
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

    // Hidden genes (the eye, or another gene's solo) are left as they are:
    // no share of the budget, no new reads, and what they have loaded stays,
    // ready to stand in when they are shown again
    for (const [gene, style] of this.genes) {
      const f = this.featureIndex.get(gene);

      if (f === undefined || !style.visible) continue;
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
    // The budgets are shared by the genes on screen. Smallest first: each
    // takes the finest level that fits an equal share of what is left, in
    // molecules and in bytes to read, so what a small gene doesn't need
    // goes to the bigger ones.
    views.sort((a, b) => a.counts[0] - b.counts[0]);
    // A new plan is a new view for the budget hooks: the clock restarts
    markLayerChanged("points");
    let budget = DETAIL_BUDGET[this.detail];
    let bytes =
      this.detail === "auto"
        ? Math.min(DETAIL_BYTES.auto, this.budget.bytesPerView)
        : DETAIL_BYTES[this.detail];
    const coarsest = this.element.levels.length - 1;

    for (const [i, { gene, f, runs, counts }] of views.entries()) {
      if (!this.genes.has(gene)) continue;
      const share = budget / (views.length - i);
      const byteShare = bytes / (views.length - i);
      let level = coarsest;

      for (let l = 0; l < runs.length; l++) {
        if (
          counts[l] <= share &&
          counts[l] * this.element.bytesPerRow(l) <= byteShare
        ) {
          level = l;
          break;
        }
      }
      budget -= counts[level];
      bytes -= counts[level] * this.element.bytesPerRow(level);
      const wanted = this.piecesInView(f, level, runs[level], rect);
      const plan = this.plans.get(gene) ?? {
        level,
        totals: [],
        wanted: [],
        standIns: new Set<Piece>(),
        allInView: 0,
      };
      const next = new Set(wanted);

      // What was loaded for the gene's last view stands in for the new one
      // until it has arrived — also when the gene was hidden meanwhile, so
      // a gene shown again (e.g. soloed, now with the whole budget) keeps
      // its coarser molecules up while the finer ones stream in
      for (const p of [...plan.wanted, ...plan.standIns]) {
        if (next.has(p)) continue;
        if (p.points && p.loaded > 0) plan.standIns.add(p);
        else plan.standIns.delete(p);
      }
      const queueUp = (p: Piece) => {
        p.lastWanted = this.frame;
        if (!p.points && !p.loading) {
          p.loading = true;
          this.queue.push(p);
        }
      };

      for (const p of wanted) {
        plan.standIns.delete(p);
        queueUp(p);
      }
      // A gene with nothing on screen yet shows its coarsest level first:
      // a few small reads, so something is drawn within the first moments
      // while the chosen level streams in
      if (
        level < coarsest &&
        plan.standIns.size === 0 &&
        !wanted.some((p) => p.loaded > 0)
      ) {
        for (const p of this.piecesInView(f, coarsest, runs[coarsest], rect)) {
          if (next.has(p)) continue;
          plan.standIns.add(p);
          queueUp(p);
        }
      }
      plan.level = level;
      plan.wanted = wanted;
      plan.totals = runs.map((r) => r.cum[r.cum.length - 1]);
      plan.allInView = counts[0];
      this.plans.set(gene, plan);
      this.refresh(gene);
    }
    // Reads for pieces no view wants any more stop here, queued or in flight
    const inUse = this.piecesInUse();

    this.queue = this.queue.filter((p) => inUse.has(p));
    for (const p of this.pieces.values()) {
      if (p.loading && !inUse.has(p)) this.cancel(p);
    }
    this.report(true);
    this.changed();
    this.pump();
    this.evict();
  }

  /** Pieces some gene's current view draws or keeps as stand-ins. */
  private piecesInUse(): Set<Piece> {
    const inUse = new Set<Piece>();

    for (const plan of this.plans.values()) {
      plan.wanted.forEach((p) => inUse.add(p));
      plan.standIns.forEach((p) => inUse.add(p));
    }

    return inUse;
  }

  /** Stop a piece's read; what it had is dropped, so it can be asked again. */
  private cancel(p: Piece) {
    p.loading = false;
    p.controller?.abort();
    p.controller = null;
    this.drop(p);
  }

  /**
   * Show the current view's pieces as far as they have arrived. Until all
   * have, earlier ones stay up, coarser ones drawn larger; then they go.
   */
  private refresh(gene: string) {
    const plan = this.plans.get(gene);

    if (!plan) return;
    if (!this.genes.get(gene)?.visible) {
      for (const p of [...plan.wanted, ...plan.standIns]) this.hide(p);

      return;
    }
    for (const p of plan.wanted) {
      if (p.points) p.points.visible = p.loaded > 0;
      this.setGrowth(p, 1);
    }
    this.report();
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

  /**
   * The E2E budget hooks (docs/specs/BUDGETS.md): what is drawn and wanted
   * across the visible genes, first content, and a settled view. A no-op
   * outside NEXT_PUBLIC_E2E=1.
   */
  private report(afterPlan = false) {
    if (!hooksEnabled()) return;
    // While a plan is being built gene by gene, the other genes still hold
    // their old (complete) views, so "settled" is only trustworthy once the
    // whole plan is in place
    const maySettle = afterPlan || !this.planning;
    let drawn = 0;
    let wanted = 0;
    let settled = this.plans.size > 0;

    for (const [gene, plan] of this.plans) {
      if (!this.genes.get(gene)?.visible) continue;
      for (const p of plan.wanted) wanted += p.count;
      for (const p of [...plan.wanted, ...plan.standIns]) {
        if (p.points?.visible) drawn += p.loaded;
      }
      if (!plan.wanted.every(done)) settled = false;
    }
    let resident = 0;

    for (const p of this.pieces.values()) if (p.points) resident += p.loaded;
    const net = this.element.network();

    markLayer("points", {
      drawn,
      wanted,
      resident,
      requests: net.requests,
      bytes: net.bytes,
    });
    if (drawn > 0) markLayerContent("points");
    if (maySettle && settled && wanted > 0) markLayerSettled("points");
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
        material: null,
        loaded: 0,
        loading: false,
        controller: null,
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

      if (!p.loading) continue; // cancelled while queued
      this.activeLoads++;
      this.load(p).finally(() => {
        this.activeLoads--;
        if (!this.disposed) this.pump();
      });
    }
  }

  /** The molecule shader's material, for a gene's style. */
  private material(style: GeneStyle) {
    if (!this.template) {
      // createSmPointCloud's material, taken from a one-point cloud
      const one = createSmPointCloud(
        new Float32Array(3),
        1,
        "#fff",
        1,
        "circle",
      );

      one.geometry.dispose();
      this.template = one.material as THREE.ShaderMaterial;
    }
    const m = this.template.clone();

    m.uniforms.dotSize.value = Math.max(style.dotSize, this.minDot);
    m.uniforms.uShape.value = style.shape === "square" ? 1 : 0;

    return m;
  }

  /**
   * Read a piece a chunk at a time. Each chunk becomes its own small point
   * cloud, uploaded once and already filled, and joins the scene when
   * uploadSome() gets to it.
   */
  private async load(p: Piece) {
    const gene = this.element.features[p.feature];

    try {
      const style = this.genes.get(gene) ?? {
        color: "#ffffff",
        dotSize: 1,
        shape: "circle" as const,
        renderOrder: 0,
        visible: true,
      };
      const group = new THREE.Group();

      group.visible = false;
      this.root.add(group);
      p.points = group;
      p.material = this.material(style);
      p.growth = 1;
      // Coarser levels first: they are small and stand in for the finer
      const controller = new AbortController();

      p.controller = controller;
      hint(controller.signal, { tag: "points", priority: p.level });

      const end = p.start + p.count;

      for (let lo = p.start; lo < end; ) {
        const hi = Math.min(
          end,
          (Math.floor(lo / STREAM_ROWS) + 1) * STREAM_ROWS,
        );
        const xyz = await this.element.readXyz(p.level, lo, hi, {
          signal: controller.signal,
        });

        if (this.disposed || p.points !== group) return;
        const n = hi - lo;
        const geometry = new THREE.BufferGeometry();

        geometry.setAttribute("position", new THREE.BufferAttribute(xyz, 3));
        // The point shader's per-point colour, size and alpha, as normalized
        // bytes; size and alpha are the same shared buffer of ones
        geometry.setAttribute(
          "color",
          new THREE.BufferAttribute(
            fillColor(new Uint8Array(n * 3), style.color),
            3,
            true,
          ),
        );
        geometry.setAttribute(
          "size",
          new THREE.BufferAttribute(ONES.subarray(0, n), 1, true),
        );
        geometry.setAttribute(
          "alpha",
          new THREE.BufferAttribute(ONES.subarray(0, n), 1, true),
        );
        geometry.computeBoundingSphere();
        const points = new THREE.Points(geometry, p.material);

        points.renderOrder = this.genes.get(gene)?.renderOrder ?? 0;
        this.uploads.push({ piece: p, points });
        lo = hi;
      }
    } catch (e) {
      if ((e as Error)?.name !== "AbortError") {
        console.warn(`Failed to load molecules ${p.key}:`, e);
      }
    } finally {
      p.loading = false;
      p.controller = null;
      this.changed();
    }
  }

  /** Put arrived chunks into the scene, up to UPLOAD_PER_FRAME molecules. */
  private uploadSome() {
    let budget = UPLOAD_PER_FRAME;
    const touched = new Set<string>();

    while (this.uploads.length && budget > 0) {
      const { piece, points } = this.uploads.shift()!;

      if (!piece.points) {
        points.geometry.dispose(); // the piece was let go meanwhile

        continue;
      }
      const n = points.geometry.getAttribute("position").count;

      piece.points.add(points);
      piece.loaded += n;
      budget -= n;
      touched.add(this.element.features[piece.feature]);
    }
    for (const gene of touched) this.refresh(gene);
    if (touched.size) this.changed();
  }

  private restyle(gene: string) {
    const f = this.featureIndex.get(gene);
    const style = this.genes.get(gene)!;

    for (const p of this.pieces.values()) {
      if (p.feature !== f || !p.points) continue;
      for (const child of p.points.children) {
        const pts = child as THREE.Points;
        const color = pts.geometry.getAttribute(
          "color",
        ) as THREE.BufferAttribute;

        fillColor(color.array as Uint8Array, style.color);
        color.needsUpdate = true;
        pts.renderOrder = style.renderOrder;
      }
      if (p.material) {
        p.material.uniforms.uShape.value = style.shape === "square" ? 1 : 0;
      }
      this.applySize(p);
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
    const cap = Math.min(
      CACHE_POINTS,
      this.budget.residentBytes / BYTES_PER_POINT,
    );

    for (const p of idle) {
      if (held <= cap) break;
      held -= p.count;
      this.drop(p);
      this.pieces.delete(p.key);
    }
  }

  private drop(p: Piece) {
    p.controller?.abort();
    p.controller = null;
    if (!p.points) return;
    p.points.removeFromParent();
    for (const child of p.points.children) {
      (child as THREE.Points).geometry.dispose();
    }
    p.material?.dispose();
    p.points = null;
    p.material = null;
    p.loaded = 0;
  }
}

type Rect = Pick<ViewRect, "minX" | "minY" | "maxX" | "maxY">;

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

/** Fill an RGB byte array with one colour. */
function fillColor(bytes: Uint8Array, hex: string): Uint8Array {
  const c = new THREE.Color(hex);
  const r = Math.round(c.r * 255);
  const g = Math.round(c.g * 255);
  const b = Math.round(c.b * 255);

  for (let i = 0; i < bytes.length; i += 3) {
    bytes[i] = r;
    bytes[i + 1] = g;
    bytes[i + 2] = b;
  }

  return bytes;
}
