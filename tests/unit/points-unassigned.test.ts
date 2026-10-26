import { describe, it, expect, vi } from "vitest";
import * as THREE from "three";

import type { PointsElement } from "@/lib/spatialdata/points";
import type { GeneStyle } from "@/lib/spatialdata/PointsLayer";

import { PointsLayer } from "@/lib/spatialdata/PointsLayer";
import { LayerHost } from "@/lib/scene/LayerHost";

/**
 * A points element that marks molecules in no cell (docs/points-zarr-spec.md
 * §8.1) is drawn as MoleculeLayer draws a legacy folder's: the two kinds
 * apart, each in its own colour, the unassigned just below the assigned, and
 * either one shown or hidden without reading anything again.
 */

const flush = () => new Promise((r) => setTimeout(r, 0));

function scene() {
  const group = new THREE.Group();
  const camera = new THREE.PerspectiveCamera(75, 1, 0.1, 1000);

  camera.position.set(0, 0, 100);
  camera.lookAt(0, 0, 0);
  camera.updateMatrixWorld();
  camera.updateProjectionMatrix();
  const renderer = {
    domElement: { width: 800, height: 600 },
  } as unknown as THREE.WebGLRenderer;

  return { group, host: new LayerHost(group, camera, renderer, vi.fn()) };
}

// One gene, 100 molecules over four tiles; every third one in no cell
const N = 100;
const isUnassigned = (i: number) => i % 3 === 0;
const UNASSIGNED = Array.from({ length: N }, (_, i) => i).filter(isUnassigned)
  .length;

function element(marks: boolean) {
  const xyz = new Float32Array(N * 3);

  for (let i = 0; i < N; i++) {
    xyz[i * 3] = (i % 10) * 10 - 50;
    xyz[i * 3 + 1] = Math.floor(i / 10) * 10 - 50;
  }
  const reads = { xyz: 0, unassigned: 0 };
  const el: PointsElement = {
    name: "transcripts",
    features: ["Gene"],
    featureKey: "gene",
    dims: 3,
    levels: [
      {
        index: 0,
        nPoints: N,
        grid: { origin: [-50, -50], tileSize: 50, nx: 2, ny: 2 },
      },
    ],
    featureCounts: async () => Float64Array.from([N]),
    network: () => ({ requests: 0, bytes: 0, inFlight: 0, cacheHits: 0 }),
    bytesPerRow: () => 3,
    featureTiles: async () => ({
      tiles: Uint32Array.from([0, 1, 2, 3]),
      counts: Uint32Array.from([25, 25, 25, 25]),
      starts: Float64Array.from([0, 25, 50, 75]),
    }),
    readXyz: async (_level, lo, hi) => {
      reads.xyz++;

      return xyz.slice(lo * 3, hi * 3);
    },
    ...(marks && {
      readUnassigned: async (_level: number, lo: number, hi: number) => {
        reads.unassigned++;

        return Uint8Array.from({ length: hi - lo }, (_, i) =>
          isUnassigned(lo + i) ? 1 : 0,
        );
      },
    }),
  };

  return { el, reads };
}

const style = (over: Partial<GeneStyle> = {}): GeneStyle => ({
  color: "#ff0000",
  dotSize: 1,
  shape: "circle",
  renderOrder: -999,
  visible: true,
  assigned: true,
  unassigned: { color: "#00ff00", dotSize: 2, shape: "square" },
  ...over,
});

async function settle(host: LayerHost) {
  for (let i = 1; i < 10; i++) {
    host.tick(16 * i);
    await flush();
  }
}

/** Molecules on screen per colour, with their render orders. */
function drawn(group: THREE.Group) {
  const out = new Map<string, { n: number; order: Set<number> }>();

  group.traverse((o) => {
    if (!(o instanceof THREE.Points) || !o.parent?.visible) return;
    const m = o.material as THREE.ShaderMaterial;

    if (!m.visible) return;
    const key = `#${(m.uniforms.uColor.value as THREE.Color).getHexString()}`;
    const e = out.get(key) ?? { n: 0, order: new Set<number>() };

    e.n += o.geometry.getAttribute("position").count;
    e.order.add(o.renderOrder);
    out.set(key, e);
  });

  return out;
}

describe("a points element's molecules in no cell", () => {
  it("are drawn apart, in their own colour, just below the assigned", async () => {
    const { group, host } = scene();
    const { el } = element(true);
    const layer = new PointsLayer(el);

    expect(layer.marksUnassigned).toBe(true);
    host.attach(layer);
    layer.setGenes(new Map([["Gene", style()]]));
    await settle(host);
    const d = drawn(group);

    expect(d.get("#ff0000")).toEqual({ n: N - UNASSIGNED, order: new Set([-999]) });
    expect(d.get("#00ff00")).toEqual({ n: UNASSIGNED, order: new Set([-1000]) });
    layer.dispose();
  });

  it("are shown and hidden, and either kind alone, without reading again", async () => {
    const { group, host } = scene();
    const { el, reads } = element(true);
    const layer = new PointsLayer(el);

    host.attach(layer);
    layer.setGenes(new Map([["Gene", style({ unassigned: null })]]));
    await settle(host);
    expect([...drawn(group).keys()]).toEqual(["#ff0000"]);
    const before = { ...reads };

    // Unassigned only, then recoloured
    layer.setGenes(new Map([["Gene", style({ assigned: false })]]));
    await settle(host);
    expect([...drawn(group).keys()]).toEqual(["#00ff00"]);
    layer.setGenes(
      new Map([
        [
          "Gene",
          style({
            assigned: false,
            unassigned: { color: "#0000ff", dotSize: 2, shape: "circle" },
          }),
        ],
      ]),
    );
    await settle(host);
    expect(drawn(group).get("#0000ff")?.n).toBe(UNASSIGNED);
    expect(reads).toEqual(before);
    layer.dispose();
  });

  it("an element without the column draws every molecule as assigned", async () => {
    const { group, host } = scene();
    const { el, reads } = element(false);
    const layer = new PointsLayer(el);

    expect(layer.marksUnassigned).toBe(false);
    host.attach(layer);
    layer.setGenes(new Map([["Gene", style()]]));
    await settle(host);
    expect(drawn(group).get("#ff0000")?.n).toBe(N);
    expect(drawn(group).has("#00ff00")).toBe(false);
    expect(reads.unassigned).toBe(0);
    layer.dispose();
  });
});
