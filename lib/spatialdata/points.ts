/**
 * Reader for a points element stored as Zarr (docs/points-zarr-spec.md):
 * levels of uniform samples, each in a tile-major and a feature-major order,
 * with an index of row numbers. Only what a view needs is read: the index
 * slices of the selected features, then their rows in the tiles on screen.
 */
import type { SpatialDataStore } from "./store";
import type { ReadOptions } from "@/lib/storage/ElementStore";

// Zarr v3 string arrays need a newer zarrita than the h5ad-zarr path's.
import * as zarr from "zarrita-v3";

import { hint, hinted } from "@/lib/storage/ElementStore";

type Arr = zarr.Array<zarr.DataType, zarr.Readable>;

export interface PointsGrid {
  origin: [number, number];
  tileSize: number;
  nx: number;
  ny: number;
}

export interface PointsLevel {
  index: number;
  nPoints: number;
  grid: PointsGrid;
}

/** One feature's tiles at one level: tile ids, counts, and the rows they start at. */
export interface FeatureTiles {
  tiles: Uint32Array;
  counts: Uint32Array;
  /** feature_major row where each tile's run starts. */
  starts: Float64Array;
}

/** What the element has read: from the server, and from its own cache. */
export interface NetworkStats {
  /** Reads sent to the store (the server, or the dropped files). */
  requests: number;
  bytes: number;
  inFlight: number;
  /** Reads answered from the cache. */
  cacheHits: number;
}

export interface PointsElement {
  name: string;
  features: string[];
  featureKey: string;
  dims: 2 | 3;
  levels: PointsLevel[];
  /** Points per feature in a level (level 0 by default). */
  featureCounts(level?: number): Promise<Float64Array>;
  /** Totals of what has been read so far. */
  network(): NetworkStats;
  featureTiles(level: number, feature: number): Promise<FeatureTiles>;
  /**
   * x, y(, z) of feature_major rows [lo, hi) at a level, interleaved xyz.
   * The signal cancels the read and carries its priority (`hint`).
   */
  readXyz(
    level: number,
    lo: number,
    hi: number,
    opts?: ReadOptions,
  ): Promise<Float32Array>;
  /** Compressed bytes read per row at a level so far (a prior before any). */
  bytesPerRow(level: number): number;
}

const decode = (b: Uint8Array) => JSON.parse(new TextDecoder().decode(b));

/** Index reads: small, needed before anything else can be planned. */
const INDEX = () => hinted({ tag: "points", priority: Infinity });

// Compressed bytes a row costs before any has been read: x and y as
// float32 with blosc zstd come to 2–3 bytes on the real stores.
const BYTES_PER_ROW_PRIOR = 3;

/** Names of points elements in this format (by `points:version`). */
export async function findPointsElements(
  sd: SpatialDataStore,
): Promise<string[]> {
  const names = new Set<string>();

  for (const key of sd.keys) {
    const m = key.match(/^points\/([^/]+)\/zarr\.json$/);

    if (m) names.add(m[1]);
  }
  const out: string[] = [];

  for (const name of names) {
    const bytes = await sd.store.get(`/points/${name}/zarr.json`);

    if (bytes && decode(bytes).attributes?.["points:version"]) out.push(name);
  }

  return out.sort();
}

export async function openPointsElement(
  sd: SpatialDataStore,
  name: string,
): Promise<PointsElement> {
  const store = sd.store;
  const root = zarr.root(store as zarr.Readable).resolve(`points/${name}`);
  const meta = await store.get(`/points/${name}/zarr.json`);

  if (!meta) throw new Error(`No points element "${name}"`);
  const attrs = decode(meta).attributes;
  const featureKey: string = attrs.spatialdata_attrs.feature_key;
  const dims = (attrs.axes as string[]).length as 2 | 3;
  const levels: PointsLevel[] = attrs.multiscales.layout.map(
    (e: any, i: number) => ({
      index: i,
      nPoints: e["points:n_points"],
      grid: {
        origin: e["points:grid"].origin,
        tileSize: e["points:grid"].tile_size,
        ny: e["points:grid"].shape[0],
        nx: e["points:grid"].shape[1],
      },
    }),
  );

  const arrays = new Map<string, Promise<Arr>>();
  const open = (path: string) => {
    let a = arrays.get(path);

    if (!a) {
      a = zarr.open.v3(root.resolve(path), { kind: "array" }) as Promise<Arr>;
      arrays.set(path, a);
    }

    return a;
  };
  const read = async (
    path: string,
    lo?: number,
    hi?: number,
    signal: AbortSignal = INDEX(),
  ) => {
    const a = await open(path);
    const sel = lo === undefined ? [null] : [zarr.slice(lo, hi!)];

    return (await zarr.get(a, sel, { signal })).data as ArrayLike<
      number | bigint
    >;
  };
  // Bytes per row per level, from what the store has read for each
  const rowsRead: number[] = [];
  const levelTag = (level: number) => `points:${level}`;
  const bytesPerRow = (level: number) => {
    const b = store.stats().byTag[levelTag(level)]?.bytes ?? 0;
    const r = rowsRead[level] ?? 0;

    // Weighted with the prior until a few chunks have been read
    return (b + BYTES_PER_ROW_PRIOR * 200_000) / (r + 200_000);
  };

  // Column names → member names, per level and order
  const columnCache = new Map<string, Promise<Map<string, string>>>();
  const columns = (level: number, order: string) => {
    const key = `${level}/${order}`;
    let p = columnCache.get(key);

    if (!p) {
      p = (async () => {
        const names = (await read(`${key}/columns`)) as unknown as string[];

        return new Map(
          Array.from(names, (n, i) => [String(n), `${key}/column_${i}`]),
        );
      })();
      columnCache.set(key, p);
    }

    return p;
  };

  const fmCols = await columns(0, "feature_major");
  const featurePath = fmCols.get(featureKey);

  if (!featurePath) throw new Error(`No "${featureKey}" column`);
  const features = Array.from(
    (await read(`${featurePath}/enum_values`)) as unknown as string[],
    String,
  );

  // Small per-level index arrays, read once
  const rowsCache = new Map<number, Promise<[Float64Array, Float64Array]>>();
  const featureIndex = (level: number) => {
    let p = rowsCache.get(level);

    if (!p) {
      p = Promise.all([
        read(`${level}/index/feature_rows`),
        read(`${level}/index/feature_tile_ptr`),
      ]).then(([rows, ptr]) => [
        Float64Array.from(rows as ArrayLike<bigint>, Number),
        Float64Array.from(ptr as ArrayLike<bigint>, Number),
      ]);
      rowsCache.set(level, p);
    }

    return p;
  };

  const tilesCache = new Map<string, Promise<FeatureTiles>>();
  const interleave = async (
    level: number,
    order: string,
    lo: number,
    hi: number,
    signal?: AbortSignal,
  ) => {
    const cols = await columns(level, order);
    const axes = ["x", "y", "z"].slice(0, dims);
    const sig = signal ?? new AbortController().signal;

    hint(sig, { tag: levelTag(level) });
    const parts = await Promise.all(
      axes.map((a) => read(cols.get(a)!, lo, hi, sig) as Promise<Float32Array>),
    );

    rowsRead[level] = (rowsRead[level] ?? 0) + (hi - lo);
    const n = hi - lo;
    const out = new Float32Array(n * 3);

    for (let i = 0; i < n; i++) {
      out[i * 3] = parts[0][i];
      out[i * 3 + 1] = parts[1][i];
      out[i * 3 + 2] = dims === 3 ? parts[2][i] : 0;
    }

    return out;
  };

  return {
    name,
    features,
    featureKey,
    dims,
    levels,
    network: () => {
      const st = store.stats();
      let requests = 0;
      let bytes = 0;

      for (const [tag, t] of Object.entries(st.byTag)) {
        if (tag === "points" || tag.startsWith("points:")) {
          requests += t.requests;
          bytes += t.bytes;
        }
      }

      return {
        requests,
        bytes,
        inFlight: st.inFlight,
        cacheHits: st.cacheHits,
      };
    },
    bytesPerRow,
    async featureCounts(level = 0) {
      const [rows] = await featureIndex(level);

      return rows.subarray(1).map((v, i) => v - rows[i]);
    },
    featureTiles(level, feature) {
      const key = `${level}/${feature}`;
      let p = tilesCache.get(key);

      if (!p) {
        p = (async () => {
          const [rows, ptr] = await featureIndex(level);
          const lo = ptr[feature];
          const hi = ptr[feature + 1];

          if (hi === lo) {
            return {
              tiles: new Uint32Array(0),
              counts: new Uint32Array(0),
              starts: new Float64Array(0),
            };
          }
          const [tiles, counts] = await Promise.all([
            read(`${level}/index/feature_tile`, lo, hi),
            read(`${level}/index/feature_tile_count`, lo, hi),
          ]);
          const c = Uint32Array.from(counts as ArrayLike<number>);
          const starts = new Float64Array(c.length);
          let s = rows[feature];

          for (let i = 0; i < c.length; i++) {
            starts[i] = s;
            s += c[i];
          }

          return {
            tiles: Uint32Array.from(tiles as ArrayLike<number>),
            counts: c,
            starts,
          };
        })();
        tilesCache.set(key, p);
      }

      return p;
    },
    readXyz: (level, lo, hi, opts) =>
      interleave(level, "feature_major", lo, hi, opts?.signal),
  };
}
