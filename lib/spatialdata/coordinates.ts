/**
 * The coordinate systems of a SpatialData store, and where each element
 * lands in the scene.
 *
 * In SpatialData an element's `coordinateTransformations` is a list of
 * **alternatives**, one per output coordinate system (each entry names its
 * `output`); it is not a chain. A store need not have a `global` system:
 * a Visium store names one system per image, and the spots map into each
 * with a scale. So the viewer chooses one target system, and each element
 * takes the one entry that leads to it (docs/plans/SDDB.md §2):
 *
 * - `global`, when every element reaches it;
 * - else the system the cells' region (the shapes element the table
 *   annotates) is identity into, which is the system the table's
 *   coordinates are in;
 * - else, with no region, the one system every element shares.
 *
 * The scene is drawn in the region's own coordinates (the table's
 * `obsm/spatial`, or the region's own centroids), so an element is placed
 * by its entry into the target followed by the inverse of the region's. A
 * transform is a scale, translation, identity, affine (the 2-D part, which
 * may rotate: a MALDI optical image sits at an angle to its pixels) or a
 * sequence of them; an element with no entry into the target is not placed
 * and says why. Points and shapes are drawn as written, so their placement
 * must still be the identity.
 *
 * A list whose entries name no output (plain OME-NGFF) is a chain, applied
 * in order: that is how a standalone OME-Zarr image is read.
 */
import type { SpatialDataStore } from "./store";

import { defaultTableName, elementNames } from "./store";
import { imageMultiscale } from "./ome";

export interface Transform {
  type: string;
  scale?: number[];
  translation?: number[];
  /** Of an "affine": one row per output axis, one column per input axis plus the translation. */
  affine?: number[][];
  /** Of a "sequence": its steps, applied in order. */
  transformations?: Transform[];
  input?: { name?: string; axes?: { name: string }[] } | string;
  output?: { name?: string; axes?: { name: string }[] } | string;
}

/**
 * A 2-D affine map:
 *   x' = scale[0] · x + skew[0] · y + offset[0]
 *   y' = skew[1] · x + scale[1] · y + offset[1]
 * `skew` is absent (or zero) for the axis-aligned maps every scale and
 * translation produce; an affine with rotation fills it.
 */
export interface Map2 {
  scale: [number, number];
  offset: [number, number];
  skew?: [number, number];
}

export const IDENTITY: Map2 = { scale: [1, 1], offset: [0, 0] };

const skewOf = (m: Map2): [number, number] => m.skew ?? [0, 0];

/** The map applied to a point. */
export function applyMap(m: Map2, x: number, y: number): [number, number] {
  const [b, c] = skewOf(m);

  return [
    m.scale[0] * x + b * y + m.offset[0],
    c * x + m.scale[1] * y + m.offset[1],
  ];
}

/** True when the map keeps x and y axes as they are (no rotation or shear). */
export function isAxisAligned(m: Map2): boolean {
  const [b, c] = skewOf(m);

  return Math.abs(b) < 1e-12 && Math.abs(c) < 1e-12;
}

/** How much the map scales lengths: √|det| (one number, even when rotated). */
export function linearScale(m: Map2): number {
  const [b, c] = skewOf(m);

  return Math.sqrt(Math.abs(m.scale[0] * m.scale[1] - b * c));
}

export type Placement = { ok: true; map: Map2 } | { ok: false; reason: string };

export type TargetRule =
  | "global"
  | "region identity"
  | "only common system"
  | "none";

/** An element's transforms and the names of its intrinsic axes. */
export interface ElementTransforms {
  transforms: Transform[];
  axes: string[];
}

export interface CoordinateSystems {
  /** Every system an element declares, sorted. */
  systems: string[];
  target: string | null;
  rule: TargetRule;
  /** The element whose coordinates the scene is drawn in (`shapes/<name>`). */
  region: string | null;
  /** Where an element (`images/<name>`, …) lands in the scene. */
  place(path: string): Placement;
}

export function outputName(t: Transform): string | null {
  const o = t.output;

  if (typeof o === "string") return o;

  return o?.name ?? null;
}

function inputAxes(t: Transform): string[] | null {
  const i = t.input;

  return typeof i === "object" && i?.axes ? i.axes.map((a) => a.name) : null;
}

/** A map from its 2 × 2 linear part and its offset; `skew` only when needed. */
function map2(
  a: number,
  b: number,
  c: number,
  d: number,
  tx: number,
  ty: number,
): Map2 {
  const m: Map2 = { scale: [a, d], offset: [tx, ty] };

  if (Math.abs(b) >= 1e-12 || Math.abs(c) >= 1e-12) m.skew = [b, c];

  return m;
}

/** First step, then the second. */
export function compose(a: Map2, b: Map2): Map2 {
  const [ab, ac] = skewOf(a);
  const [bb, bc] = skewOf(b);
  const [ox, oy] = applyMap(b, a.offset[0], a.offset[1]);

  // Linear part: B · A
  return map2(
    b.scale[0] * a.scale[0] + bb * ac,
    b.scale[0] * ab + bb * a.scale[1],
    bc * a.scale[0] + b.scale[1] * ac,
    bc * ab + b.scale[1] * a.scale[1],
    ox,
    oy,
  );
}

export function invert(m: Map2): Map2 {
  const [b, c] = skewOf(m);
  const det = m.scale[0] * m.scale[1] - b * c;

  if (!det) throw new Error("the transform is not invertible");
  const ia = m.scale[1] / det;
  const ib = -b / det;
  const ic = -c / det;
  const id = m.scale[0] / det;

  return map2(
    ia,
    ib,
    ic,
    id,
    -(ia * m.offset[0] + ib * m.offset[1]),
    -(ic * m.offset[0] + id * m.offset[1]),
  );
}

function outputAxes(t: Transform): string[] | null {
  const o = t.output;

  return typeof o === "object" && o?.axes ? o.axes.map((a) => a.name) : null;
}

/**
 * The x, y part of an "affine" step. SpatialData writes one row per output
 * axis and one column per input axis plus a last column of translations
 * (a `c` axis row and column carry the identity); the axes are named by the
 * step's input and output, else by `axes`.
 */
function affineStep(t: Transform, axes: string[]): Map2 {
  const rows = t.affine ?? [];
  const inputs = inputAxes(t) ?? axes;
  const outputs = outputAxes(t) ?? inputs;
  const xi = inputs.indexOf("x");
  const yi = inputs.indexOf("y");
  const xo = outputs.indexOf("x");
  const yo = outputs.indexOf("y");
  const at = (r: number, c: number, fill: number) =>
    r >= 0 && c >= 0 && rows[r] && c < rows[r].length ? rows[r][c] : fill;
  const last = (r: number) => (rows[r] ? rows[r][rows[r].length - 1] : 0);

  return map2(
    at(xo, xi, 1),
    at(xo, yi, 0),
    at(yo, xi, 0),
    at(yo, yi, 1),
    xo >= 0 ? last(xo) : 0,
    yo >= 0 ? last(yo) : 0,
  );
}

/**
 * A chain of scale / translation / identity / affine steps (a sequence
 * counts as its steps), on the x and y axes. `axes` names the vector
 * entries when a step does not name its own input axes. Throws on any other
 * kind of step.
 */
export function chain(transforms: Transform[], axes: string[]): Map2 {
  let m: Map2 = IDENTITY;

  for (const t of transforms) {
    const names = inputAxes(t) ?? axes;
    const pick = (v: number[], fill: number): [number, number] => {
      const at = (axis: string) => {
        const i = names.indexOf(axis);

        return i >= 0 && i < v.length ? v[i] : fill;
      };

      return [at("x"), at("y")];
    };

    if (t.type === "identity") continue;
    if (t.type === "sequence") {
      m = compose(m, chain(t.transformations ?? [], names));
    } else if (t.type === "scale" && t.scale) {
      m = compose(m, { scale: pick(t.scale, 1), offset: [0, 0] });
    } else if (t.type === "translation" && t.translation) {
      m = compose(m, { scale: [1, 1], offset: pick(t.translation, 0) });
    } else if (t.type === "affine" && t.affine) {
      m = compose(m, affineStep(t, names));
    } else {
      throw new Error(`unsupported transform "${t.type}"`);
    }
  }

  return m;
}

function isIdentity(m: Map2): boolean {
  return (
    isAxisAligned(m) &&
    Math.abs(m.scale[0] - 1) < 1e-9 &&
    Math.abs(m.scale[1] - 1) < 1e-9 &&
    Math.abs(m.offset[0]) < 1e-9 &&
    Math.abs(m.offset[1]) < 1e-9
  );
}

/** The element's map into `system`, or why there is none. */
export function toSystem(
  el: ElementTransforms,
  system: string | null,
): Placement {
  const named = el.transforms.some((t) => outputName(t) !== null);

  if (!named) {
    // Plain OME-NGFF: one chain, into whatever system the caller means
    try {
      return { ok: true, map: chain(el.transforms, el.axes) };
    } catch (e) {
      return { ok: false, reason: (e as Error).message };
    }
  }
  const entry = el.transforms.find((t) => outputName(t) === system);

  if (!entry) {
    return { ok: false, reason: `no transform into "${system}"` };
  }
  try {
    return { ok: true, map: chain([entry], el.axes) };
  } catch (e) {
    return { ok: false, reason: `${(e as Error).message} into "${system}"` };
  }
}

function outputs(el: ElementTransforms): Set<string> {
  return new Set(
    el.transforms.map(outputName).filter((s): s is string => s !== null),
  );
}

/**
 * Choose the target and place every element, from their transforms alone.
 * `region` is the element the scene's coordinates are those of.
 */
export function coordinateSystems(
  elements: Map<string, ElementTransforms>,
  region: string | null,
): CoordinateSystems {
  const named = Array.from(elements.values()).map(outputs);
  const systems = Array.from(new Set(named.flatMap((s) => [...s]))).sort();
  const every = named.length
    ? named.reduce((a, b) => new Set([...a].filter((s) => b.has(s))))
    : new Set<string>();
  const regionEl = region ? elements.get(region) : undefined;
  const regionIdentity = regionEl
    ? Array.from(outputs(regionEl))
        .filter((s) => {
          const p = toSystem(regionEl, s);

          return p.ok && isIdentity(p.map);
        })
        .sort()
    : [];

  let target: string | null = null;
  let rule: TargetRule = "none";

  if (every.has("global")) {
    [target, rule] = ["global", "global"];
  } else if (regionIdentity.length) {
    [target, rule] = [regionIdentity[0], "region identity"];
  } else if (!regionEl && every.size === 1) {
    [target, rule] = [[...every][0], "only common system"];
  } else if (!systems.length) {
    // Nothing names a system: plain transforms, all in one unnamed space
    rule = "only common system";
  }

  const regionToTarget: Placement = regionEl
    ? toSystem(regionEl, target)
    : { ok: true, map: IDENTITY };

  return {
    systems,
    target,
    rule,
    region: regionEl ? region : null,
    place(path) {
      const el = elements.get(path);

      if (!el) return { ok: false, reason: `${path} is not in the store` };
      if (systems.length && target === null) {
        return { ok: false, reason: "the store's elements share no system" };
      }
      if (!regionToTarget.ok) {
        return { ok: false, reason: `the cells: ${regionToTarget.reason}` };
      }
      const p = toSystem(el, target);

      return p.ok
        ? { ok: true, map: compose(p.map, invert(regionToTarget.map)) }
        : p;
    },
  };
}

/** True when a placement leaves x and y where they are. */
export function placesAsWritten(p: Placement): boolean {
  return p.ok && isIdentity(p.map);
}

const KINDS = ["images", "labels", "shapes", "points"] as const;

async function attrsOf(
  sd: SpatialDataStore,
  path: string,
): Promise<Record<string, any>> {
  const nodes = sd.root?.consolidated_metadata?.metadata;
  const node = nodes?.[path] ?? (await readNode(sd, path));

  return node?.attributes ?? {};
}

async function readNode(sd: SpatialDataStore, path: string) {
  const bytes = await sd.store.get(`/${path}/zarr.json`);

  return bytes ? JSON.parse(new TextDecoder().decode(bytes)) : null;
}

/** An element's transforms as its metadata declares them. */
export function elementTransforms(
  kind: (typeof KINDS)[number],
  attrs: Record<string, any>,
): ElementTransforms {
  if (kind === "images" || kind === "labels") {
    const ms = imageMultiscale(attrs) ?? {};

    return {
      transforms: ms.coordinateTransformations ?? [],
      axes: (ms.axes ?? []).map((a: { name: string } | string) =>
        typeof a === "string" ? a : a.name,
      ),
    };
  }

  return {
    transforms: attrs.coordinateTransformations ?? [],
    axes: attrs.axes ?? (kind === "points" ? ["x", "y", "z"] : ["x", "y"]),
  };
}

/**
 * The element the scene is drawn in: the shapes (else points) element the
 * default table annotates; with no table, the first shapes element, else
 * the first points element (as the scene placed images before).
 */
async function regionOf(sd: SpatialDataStore): Promise<string | null> {
  const has = (kind: string, name: string) =>
    elementNames(sd.keys, kind).includes(name);
  const table = defaultTableName(
    sd.root?.attributes ?? {},
    elementNames(sd.keys, "tables"),
  );

  if (table) {
    const region = (await attrsOf(sd, `tables/${table}`)).region;
    const names: string[] = Array.isArray(region)
      ? region
      : region
        ? [region]
        : [];

    for (const kind of ["shapes", "points"]) {
      const name = names.find((n) => has(kind, n));

      if (name) return `${kind}/${name}`;
    }
  }
  const [shapes] = elementNames(sd.keys, "shapes");

  if (shapes) return `shapes/${shapes}`;
  const [points] = elementNames(sd.keys, "points");

  return points ? `points/${points}` : null;
}

/** Read every element's transforms (from the consolidated metadata). */
export async function readCoordinateSystems(
  sd: SpatialDataStore,
): Promise<CoordinateSystems> {
  const elements = new Map<string, ElementTransforms>();

  for (const kind of KINDS) {
    for (const name of elementNames(sd.keys, kind)) {
      const path = `${kind}/${name}`;

      elements.set(path, elementTransforms(kind, await attrsOf(sd, path)));
    }
  }

  return coordinateSystems(elements, await regionOf(sd));
}

/**
 * The image to draw under the cells: the finest of those that reach the
 * target (a Visium full image; its hires image where there is no full one),
 * ties by name. Null with the reason when none can be placed.
 */
export function chooseImage(
  cs: CoordinateSystems,
  names: string[],
): { name: string } | { name: null; reason: string } {
  const placed = names
    .map((name) => ({ name, p: cs.place(`images/${name}`) }))
    .filter((x): x is { name: string; p: { ok: true; map: Map2 } } => x.p.ok);

  if (!placed.length) {
    const reasons = names.map((n) => {
      const p = cs.place(`images/${n}`);

      return `${n}: ${p.ok ? "" : p.reason}`;
    });

    return { name: null, reason: reasons.join("; ") || "no image" };
  }
  placed.sort(
    (a, b) =>
      linearScale(a.p.map) - linearScale(b.p.map) ||
      a.name.localeCompare(b.name),
  );

  return { name: placed[0].name };
}
