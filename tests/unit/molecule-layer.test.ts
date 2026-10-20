import { describe, it, expect, vi, beforeEach } from "vitest";
import * as THREE from "three";

import type { GeneVisualization } from "@/lib/stores/createSingleMoleculeVisualizationStore";
import type {
  GeneFetchOptions,
  GeneStreamOptions,
  SingleMoleculeDataset,
} from "@/lib/SingleMoleculeDataset";

// The layer talks to react-toastify directly; capture what it says and the
// close handler it installs, so a test can "click" the toast's X.
const toastMock = vi.hoisted(() => ({
  loading: vi.fn(),
  update: vi.fn(),
  dismiss: vi.fn(),
  error: vi.fn(),
}));

vi.mock("react-toastify", () => ({ toast: toastMock }));

import type { SceneHandle } from "@/lib/scene/contract";

import {
  MoleculeLayer,
  type MoleculeLayerOptions,
  type MoleculeLayerSelection,
} from "@/lib/webgl/MoleculeLayer";
import { VISUALIZATION_CONFIG } from "@/lib/config/visualization.config";
import { SM_DOT_SIZE_FACTOR } from "@/lib/webgl/point-cloud";

const BIG = VISUALIZATION_CONFIG.SINGLE_MOLECULE_STREAM_THRESHOLD;

type Kind = "assigned" | "unassigned";

/** Coordinates for `n` molecules: molecule i sits at (i, i, i). */
function coordsFor(n: number): Float32Array {
  const out = new Float32Array(n * 3);

  for (let i = 0; i < n; i++) out.set([i, i, i], i * 3);

  return out;
}

function abortError(): DOMException {
  return new DOMException("Aborted", "AbortError");
}

/** A promise a test can settle from the outside. */
function deferred<T>() {
  let resolve!: (v: T) => void;
  let reject!: (e: unknown) => void;
  const promise = new Promise<T>((res, rej) => {
    resolve = res;
    reject = rej;
  });

  return { promise, resolve, reject };
}

/**
 * The dataset as the layer sees it. Whole loads resolve immediately unless
 * `holdWhole` is set; streams are driven by hand through `pushBatch` /
 * `endStream` so a test controls exactly when molecules arrive.
 */
class FakeDataset {
  id = "fake";
  name = "fake";
  type = "fake";
  dimensions = 3 as const;
  uniqueGenes: string[];
  hasUnassigned: boolean;
  metadata: Record<string, unknown> = {};
  fetches: Record<string, number> = {};
  streams: Record<string, number> = {};
  failing = new Set<string>();
  holdWhole = new Map<string, ReturnType<typeof deferred<Float32Array>>>();
  private activeStreams = new Map<
    string,
    {
      opts: GeneStreamOptions;
      done: ReturnType<typeof deferred<Float32Array>>;
      delivered: Float32Array[];
    }
  >();

  constructor(
    private counts: Record<string, { assigned: number; unassigned?: number }>,
    opts: { streaming?: boolean; hasUnassigned?: boolean } = {},
  ) {
    this.uniqueGenes = Object.keys(counts);
    this.hasUnassigned = opts.hasUnassigned ?? false;
    if (opts.streaming) {
      // Own properties, as the S3 factories install them.
      Object.defineProperty(this, "streamCoordinatesByGene", {
        value: (gene: string, o: GeneStreamOptions) =>
          this.stream(gene, "assigned", o),
        enumerable: true,
      });
      Object.defineProperty(this, "streamUnassignedCoordinatesByGene", {
        value: (gene: string, o: GeneStreamOptions) =>
          this.stream(gene, "unassigned", o),
        enumerable: true,
      });
    }
  }

  moleculeCountFor(gene: string, kind: Kind = "assigned"): number | null {
    const c = this.counts[gene];

    if (!c) return null;

    return kind === "assigned" ? c.assigned : (c.unassigned ?? null);
  }

  async resolveMoleculeCount(gene: string, kind: Kind = "assigned") {
    return this.moleculeCountFor(gene, kind);
  }

  private whole(gene: string, kind: Kind, opts?: GeneFetchOptions) {
    const key = `${gene}:${kind}`;

    this.fetches[key] = (this.fetches[key] ?? 0) + 1;
    if (this.failing.has(key)) return Promise.reject(new Error("boom"));
    const held = this.holdWhole.get(key);
    const n = this.moleculeCountFor(gene, kind) ?? 0;

    if (held) {
      opts?.signal?.addEventListener("abort", () => held.reject(abortError()));

      return held.promise;
    }
    opts?.onProgress?.(12, 24);

    return Promise.resolve(coordsFor(n));
  }

  getCoordinatesByGene(gene: string, opts?: GeneFetchOptions) {
    return this.whole(gene, "assigned", opts);
  }

  getUnassignedCoordinatesByGene(gene: string, opts?: GeneFetchOptions) {
    return this.whole(gene, "unassigned", opts);
  }

  // The class's own (non-streaming) behaviour, as SingleMoleculeDataset has:
  // hand over everything in one batch.
  async streamCoordinatesByGene(gene: string, opts: GeneStreamOptions) {
    const coords = await this.getCoordinatesByGene(gene);

    await opts.onBatch({ coords, loadedMolecules: coords.length / 3 });

    return coords;
  }

  async streamUnassignedCoordinatesByGene(gene: string, opts: GeneStreamOptions) {
    const coords = await this.getUnassignedCoordinatesByGene(gene);

    await opts.onBatch({ coords, loadedMolecules: coords.length / 3 });

    return coords;
  }

  private stream(gene: string, kind: Kind, opts: GeneStreamOptions) {
    const key = `${gene}:${kind}`;

    this.streams[key] = (this.streams[key] ?? 0) + 1;
    if (this.failing.has(key)) return Promise.reject(new Error("boom"));
    const done = deferred<Float32Array>();

    this.activeStreams.set(key, { opts, done, delivered: [] });
    opts.signal?.addEventListener("abort", () => done.reject(abortError()));

    return done.promise;
  }

  /** Deliver one batch of `n` molecules to an open stream. */
  async pushBatch(gene: string, n: number, kind: Kind = "assigned") {
    const s = this.activeStreams.get(`${gene}:${kind}`)!;
    const coords = coordsFor(n);

    s.delivered.push(coords);
    await s.opts.onBatch({
      coords,
      loadedMolecules: s.delivered.reduce((a, c) => a + c.length / 3, 0),
    });
  }

  /** Close an open stream as the download finishing would. */
  endStream(gene: string, kind: Kind = "assigned") {
    const s = this.activeStreams.get(`${gene}:${kind}`)!;
    const total = s.delivered.reduce((a, c) => a + c.length, 0);
    const all = new Float32Array(total);
    let off = 0;

    for (const c of s.delivered) {
      all.set(c, off);
      off += c.length;
    }
    s.done.resolve(all);
  }

  streamSignal(gene: string, kind: Kind = "assigned"): AbortSignal | undefined {
    return this.activeStreams.get(`${gene}:${kind}`)?.opts.signal;
  }

  asDataset(): SingleMoleculeDataset {
    return this as unknown as SingleMoleculeDataset;
  }
}

function viz(gene: string, patch: Partial<GeneVisualization> = {}): GeneVisualization {
  return {
    gene,
    color: "#ff0000",
    localScale: 1,
    showAssigned: true,
    showUnassigned: true,
    assignedShape: "circle",
    unassignedShape: "circle",
    unassignedColor: "#00ff00",
    unassignedLocalScale: 1,
    colorSynced: true,
    ...patch,
  };
}

function selection(
  genes: GeneVisualization[],
  patch: Partial<MoleculeLayerSelection> = {},
): MoleculeLayerSelection {
  return {
    selectedGenes: new Map(genes.map((g) => [g.gene, g])),
    legendOrder: genes.map((g) => g.gene),
    globalScale: 1,
    showAssigned: true,
    showUnassigned: true,
    visible: true,
    ...patch,
  };
}

const flush = () => new Promise((r) => setTimeout(r, 0));

/** A layer attached to a scene stub whose group is `parent`. */
function mount(
  dataset: SingleMoleculeDataset,
  parent: THREE.Group,
  options?: MoleculeLayerOptions,
): MoleculeLayer {
  const layer = new MoleculeLayer(dataset, options);
  const scene: SceneHandle = {
    group: parent,
    camera: new THREE.PerspectiveCamera(),
    renderer: {} as THREE.WebGLRenderer,
    view: () => null,
    invalidate: () => {},
    attach: (l) => {
      l.attach(scene);

      return () => l.detach();
    },
  };

  layer.attach(scene);

  return layer;
}

/** The onClose handler the layer gave a gene's loading toast. */
function toastClose(gene: string): (reason?: boolean) => void {
  const call = [...toastMock.loading.mock.calls]
    .reverse()
    .find(([, o]) => o.toastId === `loading-gene-${gene}`);

  return call![1].onClose;
}

beforeEach(() => {
  toastMock.loading.mockReset();
  toastMock.update.mockReset();
  toastMock.dismiss.mockReset();
  toastMock.error.mockReset();
});

describe("MoleculeLayer: whole loads", () => {
  it("draws a small gene from one download and removes it on deselect", async () => {
    const ds = new FakeDataset({ a: { assigned: 10 } }, { streaming: true });
    const parent = new THREE.Group();
    const layer = mount(ds.asDataset(), parent);

    await layer.sync(selection([viz("a")]));

    const cloud = layer.pointClouds.get("a:assigned")!;

    expect(cloud).toBeDefined();
    expect(cloud.parent?.parent).toBe(parent);
    expect(cloud.geometry.getAttribute("position").count).toBe(10);
    expect(ds.streams["a:assigned"]).toBeUndefined();
    expect(ds.fetches["a:assigned"]).toBe(1);
    // Download progress reached the toast; then it reported success.
    expect(toastMock.loading).toHaveBeenCalledTimes(1);
    expect(toastMock.update).toHaveBeenLastCalledWith(
      "loading-gene-a",
      expect.objectContaining({ type: "success" }),
    );

    const dispose = vi.spyOn(cloud.geometry, "dispose");

    await layer.sync(selection([]));
    expect(layer.pointClouds.size).toBe(0);
    expect(cloud.parent).toBeNull();
    expect(dispose).toHaveBeenCalled();
  });

  it("never refetches for an appearance change, but applies it", async () => {
    const ds = new FakeDataset({ a: { assigned: 4 } }, { streaming: true });
    const layer = mount(ds.asDataset(), new THREE.Group());

    await layer.sync(selection([viz("a")]));
    await layer.sync(
      selection([viz("a", { color: "#0000ff", localScale: 2, assignedShape: "square" })], {
        globalScale: 3,
      }),
    );

    const cloud = layer.pointClouds.get("a:assigned")!;
    const material = cloud.material as THREE.ShaderMaterial;

    expect(ds.fetches["a:assigned"]).toBe(1);
    // Colour is a uniform: a restyle uploads nothing
    expect(cloud.geometry.getAttribute("color")).toBeUndefined();
    expect((material.uniforms.uColor.value as THREE.Color).getHex()).toBe(0x0000ff);
    expect(material.uniforms.uShape.value).toBe(1);
    expect(material.uniforms.dotSize.value).toBeCloseTo(
      2 *
        3 *
        VISUALIZATION_CONFIG.SINGLE_MOLECULE_POINT_BASE_SIZE *
        SM_DOT_SIZE_FACTOR,
      5,
    );
    // Only the loading toast and its success update — no second load.
    expect(toastMock.loading).toHaveBeenCalledTimes(1);
  });

  it("adds and removes the unassigned cloud with the global and per-gene toggles", async () => {
    const ds = new FakeDataset(
      { a: { assigned: 3, unassigned: 2 } },
      { streaming: true, hasUnassigned: true },
    );
    const layer = mount(ds.asDataset(), new THREE.Group());

    await layer.sync(selection([viz("a")]));
    expect([...layer.pointClouds.keys()].sort()).toEqual(["a:assigned", "a:unassigned"]);

    await layer.sync(selection([viz("a")], { showUnassigned: false }));
    expect([...layer.pointClouds.keys()]).toEqual(["a:assigned"]);

    await layer.sync(selection([viz("a", { showAssigned: false })]));
    expect([...layer.pointClouds.keys()]).toEqual(["a:unassigned"]);
    // The unassigned cloud came back from a fresh download; the assigned one
    // will too when it is wanted again.
    expect(ds.fetches["a:unassigned"]).toBe(2);
  });

  it("orders clouds by the legend, assigned above unassigned, all below the cells", async () => {
    const ds = new FakeDataset(
      { a: { assigned: 1, unassigned: 1 }, b: { assigned: 1, unassigned: 1 } },
      { streaming: true, hasUnassigned: true },
    );
    const layer = mount(ds.asDataset(), new THREE.Group());

    await layer.sync(selection([viz("a"), viz("b")]));
    const order = (key: string) => layer.pointClouds.get(key)!.renderOrder;

    expect(order("a:unassigned")).toBeLessThan(order("a:assigned"));
    expect(order("a:assigned")).toBeLessThan(order("b:unassigned"));
    expect(order("b:assigned")).toBeLessThan(0);

    // Re-selecting a gene moves it to the end of the legend: now on top.
    await layer.sync(selection([viz("a"), viz("b")], { legendOrder: ["b", "a"] }));
    expect(order("b:assigned")).toBeLessThan(order("a:unassigned"));
  });

  it("hides and shows every cloud with the layer switch", async () => {
    const ds = new FakeDataset({ a: { assigned: 1 } }, { streaming: true });
    const layer = mount(ds.asDataset(), new THREE.Group());

    await layer.sync(selection([viz("a")], { visible: false }));
    expect(layer.pointClouds.get("a:assigned")!.visible).toBe(false);
    await layer.sync(selection([viz("a")], { visible: true }));
    expect(layer.pointClouds.get("a:assigned")!.visible).toBe(true);
  });

  it("reports a failed download and tries again on the next sync", async () => {
    const ds = new FakeDataset({ a: { assigned: 5 } }, { streaming: true });
    const layer = mount(ds.asDataset(), new THREE.Group());

    ds.failing.add("a:assigned");
    await layer.sync(selection([viz("a")]));
    await flush();
    expect(layer.pointClouds.size).toBe(0);
    expect(toastMock.update).toHaveBeenLastCalledWith(
      "loading-gene-a",
      expect.objectContaining({ type: "error" }),
    );

    ds.failing.clear();
    await layer.sync(selection([viz("a")]));
    expect(ds.fetches["a:assigned"]).toBe(2);
    expect(layer.pointClouds.get("a:assigned")).toBeDefined();
  });

  it("loads whole when the dataset cannot stream, however big the gene", async () => {
    const ds = new FakeDataset({ a: { assigned: BIG * 2 } }, { streaming: false });
    const layer = mount(ds.asDataset(), new THREE.Group());

    await layer.sync(selection([viz("a")]));
    expect(ds.fetches["a:assigned"]).toBe(1);
    expect(ds.streams["a:assigned"]).toBeUndefined();
    expect(layer.pointClouds.get("a:assigned")!.geometry.getAttribute("position").count).toBe(
      BIG * 2,
    );
  });

  it("cancels a whole download when its toast is closed, and says so", async () => {
    const ds = new FakeDataset({ a: { assigned: 5 } }, { streaming: true });
    const onCancelGene = vi.fn();
    const layer = mount(ds.asDataset(), new THREE.Group(), { onCancelGene });
    const held = deferred<Float32Array>();

    ds.holdWhole.set("a:assigned", held);
    const synced = layer.sync(selection([viz("a")]));

    await flush();
    toastClose("a")(true);
    await synced;
    await flush();
    expect(onCancelGene).toHaveBeenCalledWith("a");
    expect(layer.pointClouds.size).toBe(0);
    // autoClose / programmatic dismissals carry no reason and cancel nothing.
    ds.holdWhole.clear();
    await layer.sync(selection([viz("a")]));
    toastClose("a")();
    expect(onCancelGene).toHaveBeenCalledTimes(1);
    expect(layer.pointClouds.size).toBe(1);
  });
});

describe("MoleculeLayer: streamed loads", () => {
  it("fills a preallocated cloud batch by batch and does not hold up the sync", async () => {
    const ds = new FakeDataset({ a: { assigned: BIG } }, { streaming: true });
    const onFirstData = vi.fn();
    const onGeneRendered = vi.fn();
    const layer = mount(ds.asDataset(), new THREE.Group(), {
      onFirstData,
      onGeneRendered,
    });

    // Resolves once the stream has started, long before it finishes.
    await layer.sync(selection([viz("a")]));
    const cloud = layer.pointClouds.get("a:assigned")!;

    expect(ds.streams["a:assigned"]).toBe(1);
    expect(ds.fetches["a:assigned"]).toBeUndefined();
    expect(cloud.geometry.getAttribute("position").count).toBe(BIG);
    expect(cloud.geometry.drawRange.count).toBe(0);
    expect(onFirstData).not.toHaveBeenCalled();

    await ds.pushBatch("a", 1000);
    expect(cloud.geometry.drawRange.count).toBe(1000);
    expect(onFirstData).toHaveBeenCalledTimes(1);
    expect(toastMock.update).toHaveBeenLastCalledWith(
      "loading-gene-a",
      expect.objectContaining({ render: expect.stringContaining("1,000 / ") }),
    );

    await ds.pushBatch("a", 2000);
    expect(cloud.geometry.drawRange.count).toBe(3000);
    expect(onFirstData).toHaveBeenCalledTimes(1);

    ds.endStream("a");
    await flush();
    expect(onGeneRendered).toHaveBeenCalledWith("a", 3000);
    expect(toastMock.update).toHaveBeenLastCalledWith(
      "loading-gene-a",
      expect.objectContaining({ type: "success", render: "a: 3,000 molecules loaded" }),
    );
    expect(layer.bounds!.max.x).toBe(1999);
  });

  it("keeps a stream running across syncs and restyles the partial cloud", async () => {
    const ds = new FakeDataset({ a: { assigned: BIG }, b: { assigned: 2 } }, { streaming: true });
    const layer = mount(ds.asDataset(), new THREE.Group());

    await layer.sync(selection([viz("a")]));
    await ds.pushBatch("a", 10);
    // Selecting another gene and recolouring the streaming one must not
    // restart the download.
    await layer.sync(selection([viz("a", { color: "#0000ff" }), viz("b")]));
    expect(ds.streams["a:assigned"]).toBe(1);
    expect(ds.streamSignal("a")!.aborted).toBe(false);

    const material = layer.pointClouds.get("a:assigned")!
      .material as THREE.ShaderMaterial;

    expect((material.uniforms.uColor.value as THREE.Color).getHex()).toBe(0x0000ff);
    expect(layer.pointClouds.get("b:assigned")).toBeDefined();
  });

  it("aborts a stream on deselect, drops the partial cloud and ignores a late batch", async () => {
    const ds = new FakeDataset({ a: { assigned: BIG } }, { streaming: true });
    const parent = new THREE.Group();
    const layer = mount(ds.asDataset(), parent);

    await layer.sync(selection([viz("a")]));
    await ds.pushBatch("a", 10);
    const cloud = layer.pointClouds.get("a:assigned")!;

    await layer.sync(selection([]));
    expect(ds.streamSignal("a")!.aborted).toBe(true);
    expect(layer.pointClouds.size).toBe(0);
    expect(cloud.parent).toBeNull();
    await flush();
    expect(toastMock.dismiss).toHaveBeenCalledWith("loading-gene-a");

    // A batch that was already in flight: nothing to draw into, no throw.
    await ds.pushBatch("a", 10);
    expect((parent.children[0] as THREE.Group).children.length).toBe(0);
    expect(cloud.geometry.drawRange.count).toBe(10);
  });

  it("cancels a stream from its toast and takes the gene out of the selection", async () => {
    const ds = new FakeDataset({ a: { assigned: BIG } }, { streaming: true });
    const onCancelGene = vi.fn();
    const layer = mount(ds.asDataset(), new THREE.Group(), { onCancelGene });

    await layer.sync(selection([viz("a")]));
    await ds.pushBatch("a", 5);
    toastClose("a")(true);
    await flush();
    expect(ds.streamSignal("a")!.aborted).toBe(true);
    expect(onCancelGene).toHaveBeenCalledWith("a");
    expect(layer.pointClouds.size).toBe(0);
  });

  it("streams the unassigned molecules the same way, behind the assigned ones", async () => {
    const ds = new FakeDataset(
      { a: { assigned: 3, unassigned: BIG } },
      { streaming: true, hasUnassigned: true },
    );
    const layer = mount(ds.asDataset(), new THREE.Group());

    await layer.sync(selection([viz("a")]));
    expect(ds.fetches["a:assigned"]).toBe(1);
    expect(ds.streams["a:unassigned"]).toBe(1);
    await ds.pushBatch("a", 7, "unassigned");
    const u = layer.pointClouds.get("a:unassigned")!;

    expect(u.geometry.drawRange.count).toBe(7);
    expect(u.renderOrder).toBeLessThan(layer.pointClouds.get("a:assigned")!.renderOrder);
    expect(toastMock.update).toHaveBeenLastCalledWith(
      "loading-gene-a",
      expect.objectContaining({ render: expect.stringContaining("a (unassigned): 7 / ") }),
    );
  });

  it("stops everything on dispose", async () => {
    const ds = new FakeDataset({ a: { assigned: BIG }, b: { assigned: 3 } }, { streaming: true });
    const parent = new THREE.Group();
    const layer = mount(ds.asDataset(), parent);
    const held = deferred<Float32Array>();

    ds.holdWhole.set("b:assigned", held);
    const synced = layer.sync(selection([viz("a"), viz("b")]));

    await flush();
    layer.dispose();
    await synced;
    expect(ds.streamSignal("a")!.aborted).toBe(true);
    expect(parent.children.length).toBe(0);
    expect(layer.pointClouds.size).toBe(0);
    expect(toastMock.dismiss).toHaveBeenCalledWith("loading-gene-a");
    expect(toastMock.dismiss).toHaveBeenCalledWith("loading-gene-b");

    // A disposed layer draws nothing more.
    await layer.sync(selection([viz("a")]));
    expect(parent.children.length).toBe(0);
  });
});
