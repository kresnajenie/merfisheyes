import type {
  GeneVisualization,
  MoleculeShape,
} from "@/lib/stores/createSingleMoleculeVisualizationStore";
import type { SingleMoleculeDataset } from "@/lib/SingleMoleculeDataset";
import type { FrameInfo, LayerStatus } from "@/lib/scene/contract";

import * as THREE from "three";
import { toast } from "react-toastify";

import {
  appendToSmPointCloud,
  createSmPointCloud,
  createStreamingSmPointCloud,
  SM_DOT_SIZE_FACTOR,
  updateDotSize,
  updatePointCloudColor,
  updatePointCloudShape,
} from "./point-cloud";

import { BaseLayer } from "@/lib/scene/BaseLayer";
import { VISUALIZATION_CONFIG } from "@/lib/config/visualization.config";
import {
  geneLoadKey,
  markGeneLoadEnd,
  markGeneLoadProgress,
  markGeneLoadStart,
} from "@/lib/utils/test-hooks";

/**
 * What the layer draws, read from the single-molecule visualization store.
 * Passed whole to `sync()`; the layer works out what changed.
 */
export interface MoleculeLayerSelection {
  selectedGenes: Map<string, GeneVisualization>;
  /** Legend order: a gene later in it is drawn on top. */
  legendOrder: Iterable<string>;
  globalScale: number;
  showAssigned: boolean;
  showUnassigned: boolean;
  /** Master switch for the whole layer. */
  visible: boolean;
}

export interface MoleculeLayerOptions {
  /**
   * Called once, when the layer first has molecules on screen — the first
   * whole-loaded gene or the first streamed batch. A scene with nothing else
   * in it fits its camera here, so a slow download is visible as it fills.
   */
  onFirstData?: (layer: MoleculeLayer) => void;
  /**
   * Called when the user closes a gene's loading toast, which cancels its
   * download. The caller takes the gene out of the selection so nothing
   * restarts it.
   */
  onCancelGene?: (gene: string) => void;
  /** Called when a gene's assigned molecules are all drawn. */
  onGeneRendered?: (gene: string, count: number) => void;
}

type Kind = "assigned" | "unassigned";

/** One in-flight download, streamed or whole. */
interface Load {
  kind: Kind;
  gene: string;
  controller: AbortController;
  /** Set for a streamed load: the cloud being filled. */
  cloud: THREE.Points | null;
  /** The test-hook entry this load reports to. */
  hook: number;
  /**
   * Resolves when the load stops holding up a `sync()`: a stream once it has
   * started filling its cloud, a whole load once it has settled.
   */
  settled: Promise<void>;
  markSettled: () => void;
}

// Blue background for in-progress gene-loading toasts, so they read as
// activity rather than a neutral notification. Closing one cancels the load.
const LOADING_TOAST_STYLE = { background: "#1d4ed8", color: "#ffffff" };
const TOAST_POSITION = "bottom-left" as const;

/**
 * Render order below the single-cell cloud (which uses three's default of 0),
 * assigned above its own unassigned, later legend entries on top.
 */
function renderOrderFor(legendIndex: number, kind: Kind): number {
  return -1000 + legendIndex * 2 + (kind === "assigned" ? 1 : 0);
}

function formatMB(bytes: number): string {
  return (bytes / 1024 ** 2).toFixed(1);
}

function disposeCloud(cloud: THREE.Points): void {
  cloud.removeFromParent();
  cloud.geometry.dispose();
  (cloud.material as THREE.Material).dispose();
}

/**
 * Draws the selected genes of a single-molecule dataset as point clouds
 * inside a Three.js group: one cloud per gene for assigned molecules and one
 * for unassigned, loaded on demand, streamed in while they download when the
 * file is large, cancellable, with one progress toast per gene.
 *
 * This is the one molecule render path. The standalone molecule scene and
 * the cell scene's molecule overlay both drive it the same way: construct it
 * on the group that holds their points, call `sync()` with the current
 * selection whenever it changes, `dispose()` when the scene goes away.
 *
 * Downloads outlive `sync()` calls: selecting another gene or changing a
 * colour reconciles what is on screen without restarting anything in flight.
 * A download stops only when its gene is deselected or hidden, its toast is
 * closed, or the layer is disposed.
 */
export class MoleculeLayer extends BaseLayer {
  readonly kind = "molecules" as const;

  private clouds = new Map<string, THREE.Points>();
  private loads = new Map<string, Load>();
  /** Genes with a loading toast on screen. */
  private toasts = new Set<string>();
  private selection: MoleculeLayerSelection | null = null;
  private hadData = false;
  private failedGenes = new Set<string>();

  constructor(
    readonly dataset: SingleMoleculeDataset,
    private options: MoleculeLayerOptions = {},
  ) {
    super();
  }

  /** Streams push their own batches; nothing is done per frame. */
  tick(_frame: FrameInfo) {}

  status(): LayerStatus {
    let drawn = 0;

    for (const cloud of this.clouds.values()) {
      const n = cloud.geometry.drawRange.count;

      drawn += Number.isFinite(n)
        ? n
        : cloud.geometry.getAttribute("position").count;
    }
    const loading = this.loads.size > 0;

    return {
      phase: this.failedGenes.size
        ? "error"
        : loading
          ? "loading"
          : this.clouds.size
            ? "ready"
            : "idle",
      progress: loading ? null : this.clouds.size ? 1 : null,
      drawn,
      wanted: drawn + this.loads.size,
      resident: drawn,
      bytes: 0,
      requests: 0,
      error: this.failedGenes.size
        ? `Failed to load ${[...this.failedGenes].join(", ")}`
        : undefined,
    };
  }

  /** Cancel every download; clouds already drawn (or partly streamed) stay. */
  abort() {
    for (const [key, load] of this.loads) {
      this.loads.delete(key);
      load.controller.abort("keep");
    }
    for (const gene of this.toasts) toast.dismiss(this.toastId(gene));
    this.toasts.clear();
    this.changed();
  }

  /** The clouds on screen, keyed `${gene}:assigned` / `${gene}:unassigned`. */
  get pointClouds(): ReadonlyMap<string, THREE.Points> {
    return this.clouds;
  }

  /** Bounds of the molecules drawn so far, in the parent's space. */
  get bounds(): THREE.Box3 | null {
    const box = new THREE.Box3();

    for (const cloud of this.clouds.values()) box.expandByObject(cloud);

    return box.isEmpty() ? null : box;
  }

  /**
   * Make the scene match the selection: remove deselected genes, restyle
   * clouds that exist, start loading the rest. Resolves when every whole
   * load that is pending afterwards has settled; streams keep filling in the
   * background.
   */
  async sync(selection: MoleculeLayerSelection): Promise<void> {
    if (this.disposed) return;
    this.selection = selection;

    const { selectedGenes } = selection;
    const legendIndex = new Map<string, number>();
    let i = 0;

    for (const gene of selection.legendOrder) legendIndex.set(gene, i++);

    // Deselected genes: stop and drop both of their clouds.
    for (const key of [...this.clouds.keys(), ...this.loads.keys()]) {
      const { gene } = this.parse(key);

      if (!selectedGenes.has(gene)) this.remove(key);
    }

    const pending: Promise<void>[] = [];

    for (const [gene, viz] of selectedGenes) {
      const wanted: Record<Kind, boolean> = {
        assigned: selection.showAssigned && viz.showAssigned,
        unassigned:
          this.dataset.hasUnassigned &&
          selection.showUnassigned &&
          viz.showUnassigned,
      };

      for (const kind of ["assigned", "unassigned"] as const) {
        const key = this.key(gene, kind);

        if (!wanted[kind]) {
          this.remove(key);
          continue;
        }
        const cloud = this.clouds.get(key);

        if (cloud) {
          this.style(cloud, gene, kind, viz, legendIndex);
          continue;
        }
        const load = this.loads.get(key) ?? this.startLoad(gene, kind);

        pending.push(load.settled);
      }
    }
    // Whatever changed above (styles, removals, loads started) is one change
    this.changed();

    await Promise.all(pending);
  }

  /** Push the distance-filter uniforms shared with the cell cloud. */
  setTargetFilter(
    center: THREE.Vector3,
    radius: number,
    feather: number,
    enabled: boolean,
  ): void {
    for (const cloud of this.clouds.values()) {
      const uniforms = (cloud.material as THREE.ShaderMaterial).uniforms;

      if (!uniforms?.uTargetCenter) continue;
      (uniforms.uTargetCenter.value as THREE.Vector3).copy(center);
      uniforms.uTargetRadius.value = radius;
      uniforms.uTargetFeather.value = feather;
      uniforms.uTargetFilterEnabled.value = enabled ? 1.0 : 0.0;
    }
  }

  /** Stop every download, drop every cloud and toast. */
  dispose(): void {
    if (this.disposed) return;
    for (const key of [...this.clouds.keys(), ...this.loads.keys()]) {
      this.remove(key);
    }
    for (const gene of this.toasts) toast.dismiss(this.toastId(gene));
    this.toasts.clear();
    super.dispose();
  }

  // ── keys ──────────────────────────────────────────────────────────────

  private key(gene: string, kind: Kind): string {
    return `${gene}:${kind}`;
  }

  private parse(key: string): { gene: string; kind: Kind } {
    const at = key.lastIndexOf(":");

    return { gene: key.slice(0, at), kind: key.slice(at + 1) as Kind };
  }

  private toastId(gene: string): string {
    return `loading-gene-${gene}`;
  }

  // ── clouds ────────────────────────────────────────────────────────────

  /** Abort the load for `key`, if any, and drop its cloud. */
  private remove(key: string): void {
    const load = this.loads.get(key);

    if (load) {
      this.loads.delete(key);
      load.controller.abort();
    }
    const cloud = this.clouds.get(key);

    if (cloud) {
      this.clouds.delete(key);
      disposeCloud(cloud);
      this.changed();
    }
  }

  private dotSize(viz: GeneVisualization, kind: Kind): number {
    const localScale =
      kind === "assigned" ? viz.localScale : viz.unassignedLocalScale;
    const globalScale = this.selection?.globalScale ?? 1;

    return (
      localScale *
      globalScale *
      VISUALIZATION_CONFIG.SINGLE_MOLECULE_POINT_BASE_SIZE *
      SM_DOT_SIZE_FACTOR
    );
  }

  private colorOf(viz: GeneVisualization, kind: Kind): string {
    return kind === "assigned" ? viz.color : viz.unassignedColor;
  }

  private shapeOf(viz: GeneVisualization, kind: Kind): MoleculeShape {
    return kind === "assigned" ? viz.assignedShape : viz.unassignedShape;
  }

  /** Apply the selection's appearance to a cloud that already exists. */
  private style(
    cloud: THREE.Points,
    gene: string,
    kind: Kind,
    viz: GeneVisualization,
    legendIndex: Map<string, number>,
  ): void {
    const color = this.colorOf(viz, kind);
    const shape = this.shapeOf(viz, kind);
    const dotSize = this.dotSize(viz, kind);
    const last = cloud.userData as {
      smColor?: string;
      smShape?: MoleculeShape;
      smDotSize?: number;
    };

    // A sync runs on every store change (each slider tick); repainting a
    // colour buffer of tens of millions of molecules each time is not free.
    if (last.smDotSize !== dotSize) {
      updateDotSize(cloud, dotSize);
      last.smDotSize = dotSize;
    }
    if (last.smColor !== color) {
      updatePointCloudColor(cloud, color);
      last.smColor = color;
    }
    if (last.smShape !== shape) {
      updatePointCloudShape(cloud, shape);
      last.smShape = shape;
    }
    cloud.renderOrder = renderOrderFor(legendIndex.get(gene) ?? 0, kind);
    cloud.visible = this.selection?.visible ?? true;
  }

  /**
   * Put a cloud on screen with the appearance the selection has NOW — a
   * download may finish long after the sync that started it.
   */
  private place(key: string, cloud: THREE.Points): void {
    const { gene, kind } = this.parse(key);
    const viz = this.selection?.selectedGenes.get(gene);
    const legendIndex = new Map<string, number>();
    let i = 0;

    for (const g of this.selection?.legendOrder ?? []) legendIndex.set(g, i++);
    if (viz) this.style(cloud, gene, kind, viz, legendIndex);
    this.root.add(cloud);
    this.clouds.set(key, cloud);
    this.changed();
  }

  private noteData(): void {
    if (this.hadData) return;
    this.hadData = true;
    this.options.onFirstData?.(this);
  }

  // ── loading ───────────────────────────────────────────────────────────

  /**
   * Molecules to stream for this gene, or 0 when it should load whole:
   * streaming pays off past SINGLE_MOLECULE_STREAM_THRESHOLD, needs the
   * count up front (to size the buffer) and a dataset that streams (S3
   * backed; a locally parsed dataset already has the data in memory).
   */
  private async streamableMolecules(gene: string, kind: Kind): Promise<number> {
    const ds = this.dataset;
    const stream =
      kind === "assigned"
        ? ds.streamCoordinatesByGene
        : ds.streamUnassignedCoordinatesByGene;

    // The S3 factories install real streaming as own properties; the class's
    // own methods hand over in-memory data in one batch, which is a whole
    // load in all but name (local folders, files parsed in the browser).
    const streams =
      typeof stream === "function" &&
      Object.prototype.hasOwnProperty.call(
        ds,
        kind === "assigned"
          ? "streamCoordinatesByGene"
          : "streamUnassignedCoordinatesByGene",
      );

    if (!streams) return 0;

    const expected = await ds.resolveMoleculeCount(gene, kind);

    return expected !== null &&
      expected >= VISUALIZATION_CONFIG.SINGLE_MOLECULE_STREAM_THRESHOLD
      ? expected
      : 0;
  }

  /**
   * Start downloading one cloud. Registered in `loads` under its key so a
   * later sync neither duplicates nor restarts it; removed when it settles.
   */
  private startLoad(gene: string, kind: Kind): Load {
    const key = this.key(gene, kind);
    const controller = new AbortController();
    let markSettled = () => {};
    const settled = new Promise<void>((resolve) => {
      markSettled = resolve;
    });
    const load: Load = {
      kind,
      gene,
      controller,
      cloud: null,
      hook: 0,
      settled,
      markSettled,
    };

    this.loads.set(key, load);
    this.showLoadingToast(gene);

    (async () => {
      const viz = this.selection?.selectedGenes.get(gene);

      if (!viz) return;
      const expected = await this.streamableMolecules(gene, kind);

      if (controller.signal.aborted) return;
      if (expected > 0) {
        await this.streamLoad(load, expected, viz);
      } else {
        await this.wholeLoad(load, viz);
      }
    })()
      .then(() => {
        if (this.loads.get(key) === load) this.loads.delete(key);
        load.markSettled();
        this.failedGenes.delete(gene);
        this.settleToast(gene, key);
        this.changed();
      })
      .catch((error: unknown) => {
        if (this.loads.get(key) === load) this.loads.delete(key);
        load.markSettled();
        const hookKey = geneLoadKey(gene, kind);

        if (
          (error as Error)?.name === "AbortError" ||
          controller.signal.aborted
        ) {
          // Deselected, hidden, or cancelled from the toast: drop whatever
          // was drawn so far, quietly. `abort()` (the contract) keeps it.
          markGeneLoadEnd(hookKey, load.hook, "aborted");
          if (controller.signal.reason !== "keep") this.remove(key);
          this.dismissToastIfIdle(gene);
          this.changed();

          return;
        }
        console.error(`[MoleculeLayer] Loading ${key} failed:`, error);
        markGeneLoadEnd(hookKey, load.hook, "failed");
        this.failedGenes.add(gene);
        this.remove(key);
        this.failToast(gene, kind);
        this.changed();
      });

    return load;
  }

  /** Fill a preallocated cloud batch by batch as the file downloads. */
  private async streamLoad(
    load: Load,
    expected: number,
    viz: GeneVisualization,
  ): Promise<void> {
    const { gene, kind, controller } = load;
    const key = this.key(gene, kind);
    const hookKey = geneLoadKey(gene, kind);
    const cloud = createStreamingSmPointCloud(
      expected,
      this.colorOf(viz, kind),
      this.dotSize(viz, kind),
      this.shapeOf(viz, kind),
    );

    load.cloud = cloud;
    this.place(key, cloud);
    load.hook = markGeneLoadStart(hookKey, "stream", expected);
    // From here the cloud fills on its own; a sync need not wait for it.
    load.markSettled();

    const label = kind === "assigned" ? gene : `${gene} (unassigned)`;
    const stream =
      kind === "assigned"
        ? this.dataset.streamCoordinatesByGene.bind(this.dataset)
        : this.dataset.streamUnassignedCoordinatesByGene.bind(this.dataset);
    const coords = await stream(gene, {
      signal: controller.signal,
      onBatch: ({ coords, loadedMolecules }) => {
        // Batches already in flight can arrive after cancellation; they must
        // not repaint the scene or revive the dismissed toast.
        if (controller.signal.aborted) return;
        appendToSmPointCloud(cloud, coords);
        markGeneLoadProgress(hookKey, load.hook, loadedMolecules);
        this.noteData();
        this.changed();
        this.updateToast(
          gene,
          `${label}: ${loadedMolecules.toLocaleString()} / ${expected.toLocaleString()} molecules`,
        );
      },
    });

    if (controller.signal.aborted) {
      throw new DOMException("Aborted", "AbortError");
    }
    const total = coords.length / 3;

    markGeneLoadEnd(hookKey, load.hook, "done");
    if (kind === "assigned") this.options.onGeneRendered?.(gene, total);
  }

  /** Download the whole file, then put the cloud on screen. */
  private async wholeLoad(load: Load, viz: GeneVisualization): Promise<void> {
    const { gene, kind, controller } = load;
    const key = this.key(gene, kind);
    const hookKey = geneLoadKey(gene, kind);
    const label = kind === "assigned" ? gene : `${gene} (unassigned)`;

    load.hook = markGeneLoadStart(
      hookKey,
      "whole",
      this.dataset.moleculeCountFor(gene, kind),
    );

    const opts = {
      signal: controller.signal,
      onProgress: (loadedBytes: number, totalBytes: number | null) => {
        if (controller.signal.aborted) return;
        this.updateToast(
          gene,
          totalBytes
            ? `${label}: ${formatMB(loadedBytes)} / ${formatMB(totalBytes)} MB`
            : `${label}: ${formatMB(loadedBytes)} MB`,
        );
      },
    };
    const coords =
      kind === "assigned"
        ? await this.dataset.getCoordinatesByGene(gene, opts)
        : await this.dataset.getUnassignedCoordinatesByGene(gene, opts);

    if (controller.signal.aborted) {
      throw new DOMException("Aborted", "AbortError");
    }
    const count = coords.length / 3;

    if (count > 0) {
      const cloud = createSmPointCloud(
        coords,
        count,
        this.colorOf(viz, kind),
        this.dotSize(viz, kind),
        this.shapeOf(viz, kind),
      );

      this.place(key, cloud);
      this.noteData();
    }
    markGeneLoadProgress(hookKey, load.hook, count);
    markGeneLoadEnd(hookKey, load.hook, "done");
    if (kind === "assigned") this.options.onGeneRendered?.(gene, count);
  }

  // ── toasts ────────────────────────────────────────────────────────────

  private showLoadingToast(gene: string): void {
    if (this.toasts.has(gene)) return;
    this.toasts.add(gene);
    toast.loading(`Loading ${gene}...`, {
      toastId: this.toastId(gene),
      position: TOAST_POSITION,
      autoClose: false,
      // Loading toasts have no close button unless asked; it is the cancel.
      closeButton: true,
      style: LOADING_TOAST_STYLE,
      onClose: (reason) => {
        this.toasts.delete(gene);
        // Only the close button (reason === true) cancels; autoClose and
        // programmatic dismissals pass no reason.
        if (reason !== true || this.disposed) return;
        let cancelled = false;

        for (const kind of ["assigned", "unassigned"] as const) {
          const key = this.key(gene, kind);

          if (this.loads.has(key)) {
            cancelled = true;
            this.remove(key);
          }
        }
        if (cancelled) this.options.onCancelGene?.(gene);
      },
    });
  }

  private updateToast(gene: string, text: string): void {
    if (!this.toasts.has(gene)) return;
    toast.update(this.toastId(gene), {
      render: text,
      isLoading: true,
      autoClose: false,
      style: LOADING_TOAST_STYLE,
    });
  }

  /** A load finished: once nothing else is pending for the gene, say so. */
  private settleToast(gene: string, finishedKey: string): void {
    if (!this.toasts.has(gene) || this.loadsPendingFor(gene)) return;
    const cloud = this.clouds.get(finishedKey);
    const count = cloud ? cloud.geometry.drawRange.count : 0;
    const total = Number.isFinite(count)
      ? count
      : (cloud?.geometry.getAttribute("position")?.count ?? 0);

    toast.update(this.toastId(gene), {
      render: `${gene}: ${total.toLocaleString()} molecules loaded`,
      type: "success",
      isLoading: false,
      autoClose: 3000,
      position: TOAST_POSITION,
      style: undefined,
    });
  }

  private failToast(gene: string, kind: Kind): void {
    const label = kind === "assigned" ? gene : `${gene} (unassigned)`;

    if (this.toasts.has(gene)) {
      toast.update(this.toastId(gene), {
        render: `Failed to load ${label}`,
        type: "error",
        isLoading: false,
        autoClose: 5000,
        position: TOAST_POSITION,
        style: undefined,
      });
    } else {
      toast.error(`Failed to load ${label}`, {
        toastId: this.toastId(gene),
        position: TOAST_POSITION,
        autoClose: 5000,
      });
    }
  }

  private dismissToastIfIdle(gene: string): void {
    if (this.toasts.has(gene) && !this.loadsPendingFor(gene)) {
      toast.dismiss(this.toastId(gene));
    }
  }

  private loadsPendingFor(gene: string): boolean {
    return (
      this.loads.has(this.key(gene, "assigned")) ||
      this.loads.has(this.key(gene, "unassigned"))
    );
  }
}
