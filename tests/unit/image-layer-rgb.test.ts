/**
 * An H&E image as spatialdata writes it (uint8, channels labelled r, g, b,
 * all three planes in one chunk: `[3, 1024, 1024]`): shown as a colour
 * photograph, and each tile read once — not as three fluorescence channels,
 * of which only the first could be read (the others asked for chunks 1 and 2
 * along an axis that has one).
 */
import type { SpatialDataImage } from "@/lib/spatialdata/image";

import * as THREE from "three";
import { describe, expect, it } from "vitest";

import { LayerHost } from "@/lib/scene/LayerHost";
import { ImageLayer, isRgb } from "@/lib/spatialdata/ImageLayer";

const flush = () => new Promise((r) => setTimeout(r, 0));

function host() {
  const camera = new THREE.PerspectiveCamera(75, 1, 0.1, 1000);

  camera.position.set(0, 0, 100);
  camera.lookAt(0, 0, 0);
  camera.updateMatrixWorld();
  camera.updateProjectionMatrix();
  const renderer = {
    domElement: { width: 800, height: 600 },
  } as unknown as THREE.WebGLRenderer;

  return new LayerHost(new THREE.Group(), camera, renderer, () => {});
}

/** One level of `side` pixels in chunks of `chunk`, `channels` per chunk. */
function image(
  labels: string[],
  { side = 128, chunk = 64, perChunk = labels.length, maxValue = 255 } = {},
) {
  const reads: number[][] = [];
  const count = labels.length;
  const array = {
    shape: [count, side, side],
    chunks: [perChunk, chunk, chunk],
    is: (t: string) => t === (maxValue === 255 ? "uint8" : "uint16"),
    getChunk: async (coords: number[]) => {
      reads.push(coords);
      const data = new Uint8Array(perChunk * chunk * chunk);

      // Plane p of chunk c holds the value 10 · (channel + 1)
      for (let p = 0; p < perChunk; p++) {
        data.fill(10 * (coords[0] * perChunk + p + 1), p * chunk * chunk, (p + 1) * chunk * chunk);
      }

      return { data };
    },
  };
  const img: SpatialDataImage = {
    name: "he",
    channelLabels: labels,
    maxValue,
    placement: { scale: [1, 1], offset: [0, 0] },
    levels: [
      {
        array: array as never,
        width: side,
        height: side,
        chunkWidth: chunk,
        chunkHeight: chunk,
        pixelWidth: 1,
        pixelHeight: 1,
        originX: -side / 2,
        originY: -side / 2,
      },
    ],
  };

  return { img, reads };
}

describe("ImageLayer and colour photographs", () => {
  it("knows an RGB image by its labels, or three unlabelled 8-bit channels", () => {
    expect(isRgb(image(["r", "g", "b"]).img)).toBe(true);
    expect(isRgb(image(["Red", "Green", "Blue"]).img)).toBe(true);
    expect(isRgb(image(["Channel 0", "Channel 1", "Channel 2"]).img)).toBe(true);
    expect(isRgb(image(["0", "1", "2"]).img)).toBe(true);
    expect(isRgb(image(["DAPI", "CD45", "18S"]).img)).toBe(false);
    expect(isRgb(image(["r", "g", "b"], { maxValue: 65535 }).img)).toBe(false);
  });

  it("shows an H&E image as red, green and blue, all on, over 0–255", () => {
    const layer = new ImageLayer(image(["r", "g", "b"]).img);

    expect(layer.channels.map((c) => [c.color, c.visible, c.min, c.max])).toEqual([
      ["#ff0000", true, 0, 255],
      ["#00ff00", true, 0, 255],
      ["#0000ff", true, 0, 255],
    ]);
    expect([0, 1, 2].map((i) => layer.hasWindow(i))).toEqual([true, true, true]);
  });

  it("reads each tile's chunk once and gives every channel its own plane", async () => {
    const { img, reads } = image(["r", "g", "b"]);
    const layer = new ImageLayer(img);
    const h = host();

    h.attach(layer);
    h.tick(16);
    for (let i = 0; i < 10; i++) await flush();

    // 2 × 2 tiles, one read each, all at channel-chunk 0
    expect(reads.length).toBe(4);
    expect(reads.every(([c]) => c === 0)).toBe(true);
    const tiles = (layer as unknown as { tiles: Map<string, { textures: THREE.DataTexture[] }> }).tiles;

    for (const tile of tiles.values()) {
      expect(tile.textures.map((t) => (t.image.data as Uint8Array)[0])).toEqual([10, 20, 30]);
    }
  });

  it("still reads a channel-per-chunk image one channel at a time", async () => {
    const { img, reads } = image(["DAPI", "CD45"], { perChunk: 1, maxValue: 255 });
    const layer = new ImageLayer(img);
    const h = host();

    layer.setChannel(1, { visible: true });
    h.attach(layer);
    // Four reads at a time; the next start on the next frame
    for (let frame = 1; frame <= 4; frame++) {
      h.tick(16 * frame);
      for (let i = 0; i < 10; i++) await flush();
    }

    expect(reads.filter(([c]) => c === 0).length).toBe(4);
    expect(reads.filter(([c]) => c === 1).length).toBe(4);
  });
});
