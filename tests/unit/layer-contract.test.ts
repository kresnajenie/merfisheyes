import { describe, it, expect, vi } from "vitest";
import * as THREE from "three";

import type { Layer, SceneHandle } from "@/lib/scene/contract";
import type { SpatialDataImage } from "@/lib/spatialdata/image";
import type { PointsElement } from "@/lib/spatialdata/points";
import type { SingleMoleculeDataset } from "@/lib/SingleMoleculeDataset";

import { ImageLayer } from "@/lib/spatialdata/ImageLayer";
import { PointsLayer } from "@/lib/spatialdata/PointsLayer";
import { ShapesLayer } from "@/lib/spatialdata/ShapesLayer";
import { MoleculeLayer } from "@/lib/webgl/MoleculeLayer";
import { LayerHost } from "@/lib/scene/LayerHost";

vi.mock("react-toastify", () => ({
  toast: { loading: vi.fn(), update: vi.fn(), dismiss: vi.fn(), error: vi.fn() },
}));

/**
 * The layer contract (docs/specs/LAYER-CONTRACT.md), checked on every layer
 * the registry knows. A new element kind adds its fake and its "change" to
 * CASES and must pass everything here before it is drawn anywhere.
 */

const flush = () => new Promise((r) => setTimeout(r, 0));

/** A scene the layers can join: a group, a camera, a renderer-shaped stub. */
function scene(): { handle: SceneHandle; host: LayerHost; group: THREE.Group; invalidate: ReturnType<typeof vi.fn> } {
  const group = new THREE.Group();
  const camera = new THREE.PerspectiveCamera(75, 1, 0.1, 1000);

  camera.position.set(0, 0, 100);
  camera.lookAt(0, 0, 0);
  camera.updateMatrixWorld();
  camera.updateProjectionMatrix();
  const renderer = {
    domElement: { width: 800, height: 600 },
  } as unknown as THREE.WebGLRenderer;
  const invalidate = vi.fn();
  const host = new LayerHost(group, camera, renderer, invalidate);

  return { handle: host.handle, host, group, invalidate };
}

// ── fakes ──────────────────────────────────────────────────────────────

function fakeImage(): SpatialDataImage {
  const level = (width: number, height: number, pixel: number) => ({
    array: {
      shape: [1, height, width],
      getChunk: async () => ({
        data: new Uint16Array(64 * 64).fill(100),
      }),
    } as never,
    width,
    height,
    chunkWidth: 64,
    chunkHeight: 64,
    pixelWidth: pixel,
    pixelHeight: pixel,
    originX: -width * pixel / 2,
    originY: -height * pixel / 2,
  });

  return {
    name: "dapi",
    channelLabels: ["DAPI"],
    maxValue: 65535,
    levels: [level(256, 256, 1), level(64, 64, 4)],
  };
}

function fakePoints(): PointsElement {
  // One feature, one level, four tiles of 25 molecules each
  const n = 100;
  const xyz = new Float32Array(n * 3);

  for (let i = 0; i < n; i++) {
    xyz[i * 3] = (i % 10) * 10 - 50;
    xyz[i * 3 + 1] = Math.floor(i / 10) * 10 - 50;
  }

  return {
    name: "transcripts",
    features: ["Gene"],
    featureKey: "gene",
    dims: 3,
    levels: [
      { index: 0, nPoints: n, grid: { origin: [-50, -50], tileSize: 50, nx: 2, ny: 2 } },
    ],
    featureCounts: async () => Float64Array.from([n]),
    network: () => ({ requests: 0, bytes: 0, inFlight: 0, cacheHits: 0 }),
    featureTiles: async () => ({
      tiles: Uint32Array.from([0, 1, 2, 3]),
      counts: Uint32Array.from([25, 25, 25, 25]),
      starts: Float64Array.from([0, 25, 50, 75]),
    }),
    readXyz: async (_level, lo, hi) => xyz.subarray(lo * 3, hi * 3),
  };
}

function fakePolygons() {
  return async () => ({
    coords: Float32Array.from([0, 0, 10, 0, 10, 10, 0, 10]),
    starts: Uint32Array.from([0, 4]),
    cell: Int32Array.from([0]),
  });
}

function fakeMolecules(): SingleMoleculeDataset {
  const coords = Float32Array.from([1, 1, 1, 2, 2, 2, 3, 3, 3]);

  return {
    id: "m",
    name: "m",
    uniqueGenes: ["A"],
    hasUnassigned: false,
    dimensions: 3,
    moleculeCountFor: () => 3,
    resolveMoleculeCount: async () => 3,
    getCoordinatesByGene: async () => coords,
    getUnassignedCoordinatesByGene: async () => new Float32Array(0),
  } as unknown as SingleMoleculeDataset;
}

const viz = {
  gene: "A",
  color: "#ff0000",
  localScale: 1,
  showAssigned: true,
  showUnassigned: false,
  assignedShape: "circle" as const,
  unassignedShape: "circle" as const,
  unassignedColor: "#00ff00",
  unassignedLocalScale: 1,
  colorSynced: true,
};

interface Case {
  kind: string;
  make: () => Layer;
  /** Ask the layer to draw something; after a frame or so it has content. */
  show: (layer: Layer) => void | Promise<void>;
  /** A change that must notify subscribers and ask for a frame. */
  change: (layer: Layer) => void;
}

const CASES: Case[] = [
  {
    kind: "image",
    make: () => new ImageLayer(fakeImage()),
    show: () => {},
    change: (l) => {
      (l as ImageLayer).visible = false;
    },
  },
  {
    kind: "points",
    make: () => new PointsLayer(fakePoints()),
    show: (l) =>
      (l as PointsLayer).setGenes(
        new Map([
          [
            "Gene",
            { color: "#fff", dotSize: 1, shape: "circle", renderOrder: 0, visible: true },
          ],
        ]),
      ),
    change: (l) => (l as PointsLayer).setDetail("more"),
  },
  {
    kind: "shapes",
    make: () => new ShapesLayer(fakePolygons()),
    show: () => {},
    change: (l) => {
      (l as ShapesLayer).mode = "filled";
    },
  },
  {
    kind: "molecules",
    make: () => new MoleculeLayer(fakeMolecules()),
    show: (l) =>
      (l as MoleculeLayer).sync({
        selectedGenes: new Map([["A", viz]]),
        legendOrder: ["A"],
        globalScale: 1,
        showAssigned: true,
        showUnassigned: true,
        visible: true,
      }),
    change: (l) =>
      void (l as MoleculeLayer).sync({
        selectedGenes: new Map([["A", { ...viz, color: "#0000ff" }]]),
        legendOrder: ["A"],
        globalScale: 2,
        showAssigned: true,
        showUnassigned: true,
        visible: true,
      }),
  },
];

for (const c of CASES) {
  describe(`layer contract: ${c.kind}`, () => {
    it("names its kind and reports a well-formed status before attaching", () => {
      const layer = c.make();
      const s = layer.status();

      expect(layer.kind).toBe(c.kind);
      expect(["idle", "loading", "ready", "error"]).toContain(s.phase);
      for (const n of [s.drawn, s.wanted, s.resident, s.bytes, s.requests]) {
        expect(Number.isFinite(n) && n >= 0).toBe(true);
      }
      layer.tick({ frame: 1, now: 0, view: null }); // a tick before attach is harmless
      layer.dispose();
    });

    it("joins a scene as one child of its group, leaves it, and can join again", () => {
      const { host, group } = scene();
      const layer = c.make();

      const detach = host.attach(layer);

      expect(group.children.length).toBe(1);
      detach();
      expect(group.children.length).toBe(0);
      host.attach(layer); // a scene rebuild
      expect(group.children.length).toBe(1);
      host.tick(16);
      host.tick(32);
      layer.dispose();
      expect(group.children.length).toBe(0);
    });

    it("pushes a change to subscribers once and asks the scene for a frame", async () => {
      const { host, invalidate } = scene();
      const layer = c.make();
      const listener = vi.fn();

      host.attach(layer);
      await c.show(layer);
      host.tick(16);
      await flush();
      const unsubscribe = layer.subscribe(listener);

      invalidate.mockClear();
      listener.mockClear();
      c.change(layer);
      c.change(layer);
      await flush();
      expect(listener).toHaveBeenCalledTimes(1);
      expect(invalidate).toHaveBeenCalled();

      unsubscribe();
      listener.mockClear();
      c.change(layer);
      await flush();
      expect(listener).not.toHaveBeenCalled();
      layer.dispose();
    });

    it("accepts a budget, survives abort, and frees everything on dispose", async () => {
      const { host, group } = scene();
      const layer = c.make();
      const listener = vi.fn();

      host.attach(layer);
      layer.setBudget({ bytesPerView: 1e6, residentBytes: 1e8 });
      await c.show(layer);
      for (let i = 0; i < 5; i++) {
        host.tick(16 * i);
        await flush();
      }
      layer.abort();
      expect(["idle", "loading", "ready", "error"]).toContain(layer.status().phase);

      layer.subscribe(listener);
      layer.dispose();
      expect(group.children.length).toBe(0);
      listener.mockClear();
      expect(() => host.attach(layer)).toThrow();
      await flush();
      expect(listener).not.toHaveBeenCalled();
    });
  });
}
