import type { ImageLevel, SpatialDataImage } from "./image";

import * as THREE from "three";

export interface ImageChannel {
  label: string;
  color: string;
  visible: boolean;
  /** Display window: pixel values mapped to black and to full colour. */
  min: number;
  max: number;
}

const CHANNEL_COLORS = [
  "#3d7dff",
  "#00e060",
  "#ff40d0",
  "#ffd000",
  "#00e0ff",
  "#ff7030",
];
// Finer levels are skipped while a view would need more tiles than this
// (a tilted 3D camera can see most of the image at once).
const MAX_TILES_PER_VIEW = 12;
const MAX_CONCURRENT_LOADS = 4;
// Texture memory kept for tiles that have scrolled out of view.
const CACHE_BYTES = 512 * 1024 * 1024;

interface Tile {
  level: number;
  mesh: THREE.Mesh<THREE.BufferGeometry, THREE.ShaderMaterial>;
  chunkX: number;
  chunkY: number;
  /** Bounds in scene units. */
  x0: number;
  y0: number;
  x1: number;
  y1: number;
  textures: (THREE.DataTexture | null)[];
  loading: boolean[];
  lastWanted: number;
}

/**
 * Draws a multiscale image under the cells of a Three.js scene, fetching only
 * the pyramid level and chunks the camera currently looks at.
 *
 * One mesh per chunk; its shader blends the visible channels additively, each
 * through its own colour and display window. The coarsest level stays loaded
 * as a backdrop while finer chunks stream in on top of it.
 */
export class ImageLayer {
  readonly channels: ImageChannel[];
  readonly maxValue: number;
  visible = true;
  /** Called when channel settings change on their own (auto contrast). */
  onChange?: () => void;

  private levels: ImageLevel[];
  private tiles = new Map<string, Tile>();
  private root = new THREE.Group();
  private blank: THREE.DataTexture;
  private sharedUniforms: {
    uColor: { value: THREE.Color[] };
    uWindow: { value: THREE.Vector2[] };
    uOn: { value: boolean[] };
  };
  private fragmentShader: string;
  private queue: { tile: Tile; channel: number }[] = [];
  private activeLoads = 0;
  private frame = 0;
  private rafId = 0;
  private disposed = false;
  private windowed: boolean[];

  constructor(
    image: SpatialDataImage,
    parent: THREE.Object3D,
    private camera: THREE.PerspectiveCamera,
    private renderer: THREE.WebGLRenderer,
  ) {
    this.levels = image.levels;
    this.maxValue = image.maxValue;
    const count = image.levels[0].array.shape[0];

    this.channels = Array.from({ length: count }, (_, i) => ({
      label: image.channelLabels[i] ?? `Channel ${i}`,
      color: CHANNEL_COLORS[i % CHANNEL_COLORS.length],
      visible: i === 0,
      min: 0,
      max: image.maxValue,
    }));
    this.windowed = new Array(count).fill(false);

    this.sharedUniforms = {
      uColor: { value: this.channels.map((c) => new THREE.Color(c.color)) },
      uWindow: {
        value: this.channels.map((c) => new THREE.Vector2(c.min, c.max)),
      },
      uOn: { value: this.channels.map((c) => c.visible) },
    };
    // Sampler arrays can't be indexed by a loop variable in GLSL ES 3, so
    // the per-channel code is written out.
    this.fragmentShader = `
      precision highp float;
      precision highp usampler2D;
      ${this.channels.map((_, i) => `uniform usampler2D uTex${i};`).join("\n")}
      uniform vec3 uColor[${count}];
      uniform vec2 uWindow[${count}];
      uniform bool uOn[${count}];
      in vec2 vUv;
      out vec4 fragColor;
      void main() {
        vec3 rgb = vec3(0.0);
        ${this.channels
          .map(
            (_, i) => `if (uOn[${i}]) {
          float v = float(texture(uTex${i}, vUv).r);
          rgb += uColor[${i}] * clamp((v - uWindow[${i}].x) / max(uWindow[${i}].y - uWindow[${i}].x, 1.0), 0.0, 1.0);
        }`,
          )
          .join("\n")}
        fragColor = vec4(rgb, 1.0);
      }`;

    this.blank = this.makeTexture(
      image.maxValue === 255 ? new Uint8Array(1) : new Uint16Array(1),
      1,
      1,
    );

    this.root.renderOrder = -1000;
    parent.add(this.root);
    this.tick();
  }

  setChannel(index: number, patch: Partial<ImageChannel>) {
    const channel = Object.assign(this.channels[index], patch);

    if (patch.min !== undefined || patch.max !== undefined) {
      this.windowed[index] = true;
    }
    this.sharedUniforms.uColor.value[index].set(channel.color);
    this.sharedUniforms.uWindow.value[index].set(channel.min, channel.max);
    this.sharedUniforms.uOn.value[index] = channel.visible;
  }

  dispose() {
    this.disposed = true;
    cancelAnimationFrame(this.rafId);
    this.root.removeFromParent();
    for (const tile of this.tiles.values()) this.disposeTile(tile);
    this.tiles.clear();
    this.blank.dispose();
  }

  private tick = () => {
    if (this.disposed) return;
    this.rafId = requestAnimationFrame(this.tick);
    this.update();
  };

  /** Choose the level and chunks for the current view, and show what's loaded. */
  private update() {
    this.frame++;
    this.root.visible = this.visible && this.channels.some((c) => c.visible);
    if (!this.root.visible) return;

    const view = this.visibleRect();
    const base = this.levels.length - 1;
    let level = base;

    if (view) {
      // Finest level that is no finer than the screen, within the tile cap
      for (let l = 0; l < base; l++) {
        if (
          this.levels[l].pixelWidth >= view.unitsPerPixel * 0.75 &&
          this.chunksInView(l, view).length <= MAX_TILES_PER_VIEW
        ) {
          level = l;
          break;
        }
      }
    }

    const wanted = new Set<Tile>();

    for (const l of level === base ? [base] : [base, level]) {
      for (const [cx, cy] of this.chunksInView(l, l === base ? null : view)) {
        const tile = this.getTile(l, cx, cy);

        tile.lastWanted = this.frame;
        wanted.add(tile);
        this.channels.forEach((channel, c) => {
          if (channel.visible && !tile.textures[c] && !tile.loading[c]) {
            tile.loading[c] = true;
            this.queue.push({ tile, channel: c });
          }
        });
      }
    }

    // A tile is drawn once every visible channel has arrived. Coarser tiles
    // stay up underneath, so the view refines instead of blanking.
    for (const tile of this.tiles.values()) {
      tile.mesh.visible =
        tile.level >= level &&
        this.channels.every((ch, c) => !ch.visible || tile.textures[c]);
    }

    this.pump();
    this.evict(wanted);
  }

  /**
   * The part of the image plane (z = 0 in the parent's frame) on screen, and
   * how many scene units one screen pixel covers there. Null when the camera
   * doesn't look onto the plane.
   */
  private visibleRect() {
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

    const canvas = this.renderer.domElement;

    return {
      minX,
      minY,
      maxX,
      maxY,
      unitsPerPixel: Math.min(
        (maxX - minX) / canvas.width,
        (maxY - minY) / canvas.height,
      ),
    };
  }

  /** Chunk indices [x, y] of a level overlapping the view (all, if null). */
  private chunksInView(
    l: number,
    view: { minX: number; minY: number; maxX: number; maxY: number } | null,
  ): [number, number][] {
    const lv = this.levels[l];
    const nx = Math.ceil(lv.width / lv.chunkWidth);
    const ny = Math.ceil(lv.height / lv.chunkHeight);
    let x0 = 0;
    let y0 = 0;
    let x1 = nx - 1;
    let y1 = ny - 1;

    if (view) {
      // Scene → pixel index of this level (pixel i spans [i - 0.5, i + 0.5])
      const px = (x: number) => (x - lv.originX) / lv.pixelWidth + 0.5;
      const py = (y: number) => (y - lv.originY) / lv.pixelHeight + 0.5;

      x0 = Math.max(0, Math.floor(px(view.minX) / lv.chunkWidth));
      x1 = Math.min(nx - 1, Math.floor(px(view.maxX) / lv.chunkWidth));
      y0 = Math.max(0, Math.floor(py(view.minY) / lv.chunkHeight));
      y1 = Math.min(ny - 1, Math.floor(py(view.maxY) / lv.chunkHeight));
    }

    const out: [number, number][] = [];

    for (let cy = y0; cy <= y1; cy++) {
      for (let cx = x0; cx <= x1; cx++) out.push([cx, cy]);
    }

    return out;
  }

  private getTile(l: number, chunkX: number, chunkY: number): Tile {
    const key = `${l}/${chunkY}/${chunkX}`;
    const existing = this.tiles.get(key);

    if (existing) return existing;

    const lv = this.levels[l];
    const px0 = chunkX * lv.chunkWidth;
    const py0 = chunkY * lv.chunkHeight;
    const w = Math.min(lv.chunkWidth, lv.width - px0);
    const h = Math.min(lv.chunkHeight, lv.height - py0);
    const x0 = lv.originX + (px0 - 0.5) * lv.pixelWidth;
    const y0 = lv.originY + (py0 - 0.5) * lv.pixelHeight;
    const x1 = x0 + w * lv.pixelWidth;
    const y1 = y0 + h * lv.pixelHeight;
    // Edge chunks are stored full-size; only the valid part is mapped.
    const u = w / lv.chunkWidth;
    const v = h / lv.chunkHeight;

    const geometry = new THREE.BufferGeometry();

    geometry.setAttribute(
      "position",
      new THREE.Float32BufferAttribute(
        [x0, y0, 0, x1, y0, 0, x1, y1, 0, x0, y1, 0],
        3,
      ),
    );
    geometry.setAttribute(
      "uv",
      new THREE.Float32BufferAttribute([0, 0, u, 0, u, v, 0, v], 2),
    );
    geometry.setIndex([0, 1, 2, 0, 2, 3]);

    const uniforms: Record<string, THREE.IUniform> = { ...this.sharedUniforms };

    this.channels.forEach((_, c) => {
      uniforms[`uTex${c}`] = { value: this.blank };
    });

    const mesh = new THREE.Mesh(
      geometry,
      new THREE.ShaderMaterial({
        glslVersion: THREE.GLSL3,
        uniforms,
        vertexShader: `
          out vec2 vUv;
          void main() {
            vUv = uv;
            gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
          }`,
        fragmentShader: this.fragmentShader,
        side: THREE.DoubleSide,
        depthTest: false,
        depthWrite: false,
      }),
    );

    // Behind the points, finer levels over coarser ones
    mesh.renderOrder = -1000 - l;
    mesh.visible = false;
    this.root.add(mesh);

    const tile: Tile = {
      level: l,
      mesh,
      chunkX,
      chunkY,
      x0,
      y0,
      x1,
      y1,
      textures: this.channels.map(() => null),
      loading: this.channels.map(() => false),
      lastWanted: this.frame,
    };

    this.tiles.set(key, tile);

    return tile;
  }

  /** Start queued chunk loads, skipping any whose tile has left the view. */
  private pump() {
    while (this.activeLoads < MAX_CONCURRENT_LOADS && this.queue.length > 0) {
      const { tile, channel } = this.queue.shift()!;

      // Scrolled away (or the channel was switched off) before its turn
      if (tile.lastWanted < this.frame - 1 || !this.channels[channel].visible) {
        tile.loading[channel] = false;
        continue;
      }
      this.activeLoads++;
      this.loadChunk(tile, channel)
        .catch((e) =>
          console.error("[ImageLayer] Failed to load image chunk:", e),
        )
        .finally(() => {
          tile.loading[channel] = false;
          this.activeLoads--;
        });
    }
  }

  private async loadChunk(tile: Tile, channel: number) {
    const lv = this.levels[tile.level];
    const chunk = await lv.array.getChunk([channel, tile.chunkY, tile.chunkX]);

    if (this.disposed || !this.tiles.has(this.keyOf(tile))) return;

    const data = chunk.data as Uint8Array | Uint16Array;

    if (tile.level === this.levels.length - 1 && !this.windowed[channel]) {
      this.windowed[channel] = true;
      this.setChannel(channel, { min: 0, max: percentile(data, 0.999) });
      this.onChange?.();
    }

    const texture = this.makeTexture(data, lv.chunkWidth, lv.chunkHeight);

    tile.textures[channel] = texture;
    tile.mesh.material.uniforms[`uTex${channel}`].value = texture;
  }

  private makeTexture(
    data: Uint8Array | Uint16Array,
    width: number,
    height: number,
  ) {
    const is8 = data instanceof Uint8Array;
    const texture = new THREE.DataTexture(
      data,
      width,
      height,
      THREE.RedIntegerFormat,
      is8 ? THREE.UnsignedByteType : THREE.UnsignedShortType,
    );

    texture.internalFormat = is8 ? "R8UI" : "R16UI";
    // Integer textures can't be filtered
    texture.minFilter = THREE.NearestFilter;
    texture.magFilter = THREE.NearestFilter;
    texture.unpackAlignment = 1;
    texture.needsUpdate = true;

    return texture;
  }

  /** Drop the textures of out-of-view tiles, oldest first, over the budget. */
  private evict(wanted: Set<Tile>) {
    const bytesOf = (tile: Tile) =>
      tile.textures.reduce(
        (sum, t) => sum + (t ? (t.image.data as Uint16Array).byteLength : 0),
        0,
      );
    let total = 0;

    for (const tile of this.tiles.values()) total += bytesOf(tile);
    if (total <= CACHE_BYTES) return;

    const idle = Array.from(this.tiles.values())
      .filter((t) => !wanted.has(t) && !t.loading.some(Boolean))
      .sort((a, b) => a.lastWanted - b.lastWanted);

    for (const tile of idle) {
      if (total <= CACHE_BYTES) break;
      total -= bytesOf(tile);
      this.disposeTile(tile);
      this.tiles.delete(this.keyOf(tile));
    }
  }

  private keyOf(tile: Tile) {
    return `${tile.level}/${tile.chunkY}/${tile.chunkX}`;
  }

  private disposeTile(tile: Tile) {
    this.root.remove(tile.mesh);
    tile.mesh.geometry.dispose();
    tile.mesh.material.dispose();
    for (const texture of tile.textures) texture?.dispose();
  }
}

/** Value below which the given fraction of the non-zero pixels fall. */
function percentile(data: Uint8Array | Uint16Array, fraction: number): number {
  const histogram = new Uint32Array(65536);
  let count = 0;

  for (let i = 0; i < data.length; i++) {
    if (data[i] > 0) {
      histogram[data[i]]++;
      count++;
    }
  }

  let seen = 0;

  for (let v = 1; v < histogram.length; v++) {
    seen += histogram[v];
    if (seen >= count * fraction) return v;
  }

  return 1;
}
