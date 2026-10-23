/**
 * The viewer against the shapes of store SpatialData-DB added after RNA:
 * CODEX / PhenoCycler (a dense X, dozens of image channels), MALDI-MSI (a
 * feature index column that is not `_index`, an image placed by a rotating
 * affine) and CosMx (no `obsm/spatial`: the cells' positions are in their
 * shapes). Fixtures: scripts/testdata/make-tiny-modalities.py, built by the
 * `spatialdata` library; the real roots they copy are in tests/data/sddb/roots.
 */
import type { SpatialDataImage } from "@/lib/spatialdata/image";
import type { SpatialDataStore } from "@/lib/spatialdata/store";

import { readdirSync } from "node:fs";
import { open, readFile } from "node:fs/promises";
import path from "node:path";

import * as THREE from "three";
import { describe, expect, it } from "vitest";

import { LayerHost } from "@/lib/scene/LayerHost";
import { centroidSource } from "@/lib/spatialdata/centroids";
import { ImageLayer, MAX_VISIBLE_CHANNELS } from "@/lib/spatialdata/ImageLayer";
import { imageBounds } from "@/lib/spatialdata/image";
import { SpatialDataTableAdapter } from "@/lib/spatialdata/SpatialDataTableAdapter";
import { centroid, outerRing } from "@/lib/spatialdata/wkb";

const DIR = path.resolve(__dirname, "../data/tiny/zarr");

type Range = { offset: number; length: number } | { suffixLength: number };

async function readRange(file: string, range: Range) {
  const handle = await open(file).catch(() => null);

  if (!handle) return undefined;
  try {
    const { size } = await handle.stat();
    const offset =
      "suffixLength" in range ? size - range.suffixLength : range.offset;
    const length = "suffixLength" in range ? range.suffixLength : range.length;
    const { buffer } = await handle.read(Buffer.alloc(length), 0, length, offset);

    return new Uint8Array(buffer);
  } finally {
    await handle.close();
  }
}

/** A tiny store on disk, read the way a public bucket is. */
function diskStore(name: string) {
  const root = path.join(DIR, name);
  const keys = (readdirSync(root, { recursive: true }) as string[]).map((k) =>
    k.split(path.sep).join("/"),
  );
  const store = {
    get: (key: string) =>
      readFile(path.join(root, key)).then(
        (b) => new Uint8Array(b),
        () => undefined,
      ),
    getRange: (key: string, range: Range) =>
      readRange(path.join(root, key), range),
  };

  return { root, keys, store };
}

async function openAdapter(name: string) {
  const { keys, store } = diskStore(name);
  const adapter = new SpatialDataTableAdapter(store as any, keys, undefined);

  await adapter.initialize();

  return adapter;
}

describe("MALDI-MSI table", () => {
  it("names its features from the var index column the dataframe declares", async () => {
    const adapter = await openAdapter("maldi_v3.zarr");

    expect(await adapter.loadGenes()).toEqual([
      "C6H12O6+H",
      "C6H12O6+Na",
      "C18H34O2-H",
      "C5H9NO4+K",
      "C27H46O+H",
    ]);
    expect(adapter.hasSpatialCoordinates()).toBe(true);
    expect(adapter.getDatasetInfo()).toMatchObject({ numCells: 64, numGenes: 5, xFormat: "csc" });
    const ion = await adapter.fetchGeneExpression("C6H12O6+Na");

    expect(ion).toHaveLength(64);
  });
});

describe("CODEX table", () => {
  it("reads a dense X one marker at a time", async () => {
    const adapter = await openAdapter("codex_v3.zarr");

    expect(adapter.getDatasetInfo()).toMatchObject({ numCells: 36, numGenes: 6, xFormat: "dense" });
    expect(await adapter.loadGenes()).toEqual(["DAPI", "CD3", "CD20", "blank", "PanCK", "blank (2)"]);
    const cd3 = await adapter.fetchGeneExpression("CD3");

    expect(cd3).toHaveLength(36);
    expect(cd3!.every((v) => v > 0)).toBe(true);
    expect(adapter.getClusterColumnInfo().names).toEqual(["area"]);
  });
});

describe("CosMx table", () => {
  it("has no coordinates of its own and takes them from the vitessce points", async () => {
    const adapter = await openAdapter("cosmx_v3.zarr");

    expect(adapter.hasSpatialCoordinates()).toBe(false);
    await expect(adapter.loadSpatialCoordinates()).rejects.toThrow(/No spatial coordinates/);
    expect(await adapter.loadGenes()).toEqual(["PanCK", "CD45", "CD68", "Membrane"]);
    expect(await adapter.readInstanceIds()).toHaveLength(25);

    const { keys } = diskStore("cosmx_v3.zarr");
    const sd = { keys } as unknown as SpatialDataStore;

    expect(centroidSource(sd, "shapes/cells_polygons")).toBe("shapes/cells_polygons_for_vitessce");
    expect(centroidSource(sd, "shapes/other")).toBe("shapes/other");
    expect(centroidSource(sd, "points/transcripts")).toBeNull();
    expect(centroidSource(sd, null)).toBeNull();
  });
});

describe("WKB centroids", () => {
  const le = (writer: (v: DataView) => void, size: number) => {
    const bytes = new Uint8Array(size);

    writer(new DataView(bytes.buffer));

    return bytes;
  };

  it("takes a Point as it is", () => {
    const point = le((v) => {
      v.setUint8(0, 1);
      v.setUint32(1, 1, true);
      v.setFloat64(5, 12.5, true);
      v.setFloat64(13, -3, true);
    }, 21);

    expect(centroid(point)).toEqual([12.5, -3]);
    expect(outerRing(point)).toBeNull();
  });

  it("takes the vertex mean of a polygon, the closing vertex once", () => {
    // A square 0,0 → 4,0 → 4,2 → 0,2, closed
    const square = le((v) => {
      v.setUint8(0, 1);
      v.setUint32(1, 3, true);
      v.setUint32(5, 1, true); // one ring
      v.setUint32(9, 5, true); // five vertices
      const pts = [0, 0, 4, 0, 4, 2, 0, 2, 0, 0];

      pts.forEach((p, i) => v.setFloat64(13 + i * 8, p, true));
    }, 13 + 10 * 8);

    expect(centroid(square)).toEqual([2, 1]);
    expect(Array.from(outerRing(square)!)).toEqual([0, 0, 4, 0, 4, 2, 0, 2, 0, 0]);
  });
});

/** A fake image of `count` channels, one per chunk, at a placement. */
function fakeImage(count: number, placement: SpatialDataImage["placement"]) {
  const side = 64;
  const array = {
    shape: [count, side, side],
    chunks: [1, side, side],
    is: (t: string) => t === "uint16",
    getChunk: async (coords: number[]) => ({
      data: new Uint16Array(side * side).fill(100 * (coords[0] + 1)),
    }),
  };

  return {
    name: "codex",
    channelLabels: Array.from({ length: count }, (_, i) => `marker${i}`),
    maxValue: 65535,
    placement,
    levels: [
      {
        array: array as never,
        width: side,
        height: side,
        chunkWidth: side,
        chunkHeight: side,
        pixelWidth: 1,
        pixelHeight: 1,
        originX: 0,
        originY: 0,
      },
    ],
  } satisfies SpatialDataImage;
}

function host() {
  const camera = new THREE.PerspectiveCamera(75, 1, 0.1, 1000);

  camera.position.set(32, 32, 100);
  camera.lookAt(32, 32, 0);
  camera.updateMatrixWorld();
  camera.updateProjectionMatrix();
  const renderer = { domElement: { width: 800, height: 600 } } as unknown as THREE.WebGLRenderer;

  return new LayerHost(new THREE.Group(), camera, renderer, () => {});
}

const flush = () => new Promise((r) => setTimeout(r, 0));

describe("ImageLayer with many channels", () => {
  it("compiles one shader slot per shown channel, never more than the cap", async () => {
    const layer = new ImageLayer(fakeImage(84, { scale: [1, 1], offset: [0, 0] }));
    const shader = (layer as any).fragmentShader as string;

    expect(layer.channels).toHaveLength(84);
    expect(shader.match(/usampler2D uTex\d+;/g)).toHaveLength(MAX_VISIBLE_CHANNELS);
    expect(layer.channels.filter((c) => c.visible).map((c) => c.label)).toEqual(["marker0"]);

    // Turn nine on: the first gives way, eight stay
    for (let i = 1; i <= MAX_VISIBLE_CHANNELS; i++) layer.setChannel(i, { visible: true });
    const on = layer.channels.map((c, i) => (c.visible ? i : -1)).filter((i) => i >= 0);

    expect(on).toEqual([1, 2, 3, 4, 5, 6, 7, 8]);
    expect((layer as any).slots).toEqual([1, 2, 3, 4, 5, 6, 7, 8]);
    expect((layer as any).sharedUniforms.uOn.value.every(Boolean)).toBe(true);
  });

  it("draws a tile once its shown channels have arrived, through their slots", async () => {
    const image = fakeImage(20, { scale: [1, 1], offset: [0, 0] });
    const layer = new ImageLayer(image);
    const h = host();

    h.attach(layer);
    layer.setChannel(5, { visible: true });
    h.tick({ view: { minX: 0, minY: 0, maxX: 64, maxY: 64, unitsPerPixel: 0.1 } } as never);
    await flush();
    await flush();
    const tile = [...(layer as any).tiles.values()][0];

    expect(tile.textures[0]).toBeTruthy();
    expect(tile.textures[5]).toBeTruthy();
    expect(tile.textures[1]).toBeNull();
    expect(tile.mesh.material.uniforms.uTex1.value).toBe(tile.textures[5]);
    h.tick({ view: { minX: 0, minY: 0, maxX: 64, maxY: 64, unitsPerPixel: 0.1 } } as never);
    expect(tile.mesh.visible).toBe(true);
    layer.dispose();
  });

  it("carries a rotating placement as the group's matrix and maps the view back", async () => {
    // Image units → scene: scale 2, turned 90° (x' = -2y + 100, y' = 2x)
    const placement = { scale: [0, 0] as [number, number], offset: [100, 0] as [number, number], skew: [-2, 2] as [number, number] };
    const image = fakeImage(1, placement);
    const layer = new ImageLayer(image);

    expect(imageBounds(image)).toEqual({ minX: -27, minY: -1, maxX: 101, maxY: 127 });
    const m = layer.root.matrix.elements;

    // Column-major: the first column is where the image x axis goes
    expect([m[0], m[1]]).toEqual([0, 2]);
    expect([m[4], m[5]]).toEqual([-2, 0]);
    expect([m[12], m[13]]).toEqual([100, 0]);
    // A scene view over the placed image asks for the image's one tile
    const view = { minX: -27, minY: -1, maxX: 101, maxY: 127, unitsPerPixel: 0.4 };
    const local = (layer as any).toImageView(view);

    // The image's one tile spans −0.5 … 63.5 on both axes in its own units
    expect(local.minX).toBeCloseTo(-0.5);
    expect(local.maxX).toBeCloseTo(63.5);
    expect(local.minY).toBeCloseTo(-0.5);
    expect(local.maxY).toBeCloseTo(63.5);
    expect(local.unitsPerPixel).toBeCloseTo(0.2);
  });
});
