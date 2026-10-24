#!/usr/bin/env python3
"""Check a SpatialData store on disk against the MERFISHEYES profiles
(docs/specs/store-layout.md, table-profile.md, shapes-profile.md,
image-profile.md, points-zarr-spec.md, mapping.md).

  python scripts/validate_store.py <store.zarr> [--mapping mapping.json]
                                   [--sample N] [--json]
  python scripts/validate_store.py --root-json <zarr.json> [--json]

With --root-json the store is its root `zarr.json` alone (a store read
through a URL that cannot be listed): every check the consolidated
metadata can answer runs; the parquet files and the points invariants are
reported as not checked unless their footers are given (`RootStore`).

Prints one line per check — OK, WARN or FAIL — and exits 1 when any MUST
fails. Needs numpy, zarr >= 3, jsonschema; the points checks use
scripts/points_zarr.py.
"""
from __future__ import annotations

import argparse
import gzip
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCHEMAS = HERE.parent / "docs" / "specs" / "schemas"
# The names particular producers use, shared with the viewer
# (lib/spatialdata/conventions.ts); docs/specs/conventions.md
CONVENTIONS = json.loads((HERE.parent / "docs" / "specs" / "conventions.json").read_text())
sys.path.insert(0, str(HERE))


class Report:
    def __init__(self):
        self.lines: list[tuple[str, str, str]] = []

    def ok(self, where, msg):
        self.lines.append(("OK", where, msg))

    def warn(self, where, msg):
        self.lines.append(("WARN", where, msg))

    def fail(self, where, msg):
        self.lines.append(("FAIL", where, msg))

    @property
    def failed(self):
        return any(level == "FAIL" for level, _, _ in self.lines)


def read_json(path: Path):
    if path.suffix == ".gz":
        with gzip.open(path, "rt") as f:
            return json.load(f)
    return json.loads(path.read_text())


class RootStore:
    """A store known only by its root `zarr.json` (consolidated metadata).

    `parquet` maps an element path (`shapes/<name>`, `points/<name>`) to what
    is known of its parquet: `("file", FileMetaData)`, `("dir", FileMetaData
    of one part or None)`, or nothing when it was not read.
    """

    def __init__(self, root: dict, parquet: dict | None = None):
        self.root = root
        self.nodes = (root.get("consolidated_metadata") or {}).get("metadata") or {}
        self.parquet = parquet or {}


def node(store, path: str) -> dict | None:
    if isinstance(store, RootStore):
        return store.root if not path else store.nodes.get(path)
    f = store / path / "zarr.json" if path else store / "zarr.json"
    return read_json(f) if f.is_file() else None


def children(store, path: str) -> list[str]:
    if isinstance(store, RootStore):
        prefix = f"{path}/"
        return sorted({k[len(prefix):].split("/")[0] for k in store.nodes if k.startswith(prefix)})
    d = store / path
    return sorted(p.name for p in d.iterdir() if p.is_dir() and (p / "zarr.json").is_file()) if d.is_dir() else []


def parquet_meta(store, where: str, filename: str):
    """("file" | "dir" | None, FileMetaData | None) for `<where>/<filename>`."""
    if isinstance(store, RootStore):
        return store.parquet.get(where, (None, None))
    f = store / where / filename
    if f.is_dir():
        return "dir", None
    if not f.is_file():
        return None, None
    import pyarrow.parquet as pq

    return "file", pq.read_metadata(f)


def chunk_shape(meta: dict) -> list[int]:
    """The inner chunk shape: the sharding codec's, when sharded."""
    for codec in meta.get("codecs", []):
        if codec.get("name") == "sharding_indexed":
            return list(codec["configuration"]["chunk_shape"])
    return list(meta["chunk_grid"]["configuration"]["chunk_shape"])


def schema_errors(instance, schema_name: str) -> list[str]:
    import jsonschema

    schema = json.loads((SCHEMAS / schema_name).read_text())
    v = jsonschema.Draft202012Validator(schema)
    return [f"{'/'.join(str(p) for p in e.absolute_path) or '/'}: {e.message}" for e in v.iter_errors(instance)]


# ── the store ───────────────────────────────────────────────────────────

def check_root(store, r: Report) -> dict | None:
    root = node(store, "")
    if not root:
        r.fail("root", "no zarr.json at the store root")
        return None
    if root.get("zarr_format") != 3:
        r.fail("root", f"zarr_format is {root.get('zarr_format')}, not 3")
    attrs = root.get("attributes", {})
    if "spatialdata_attrs" not in attrs:
        r.fail("root", "attributes lack spatialdata_attrs")
    else:
        r.ok("root", f"spatialdata_attrs {attrs['spatialdata_attrs']}")
    nodes = root.get("consolidated_metadata", {}).get("metadata")
    if not isinstance(nodes, dict):
        r.fail("root", "no consolidated metadata: a store opened by URL cannot be listed")
        return root
    if isinstance(store, RootStore):
        r.ok("root", f"consolidated metadata names {len(nodes)} nodes (root only: not compared with the files)")
        return root
    on_disk = sorted(
        str(p.parent.relative_to(store)).replace("\\", "/")
        for p in store.rglob("zarr.json")
        if p.parent != store
    )
    missing = [p for p in on_disk if p not in nodes]
    if missing:
        r.fail("root", f"consolidated metadata misses {len(missing)} node(s), e.g. {missing[:3]}: consolidate after the last write")
    else:
        r.ok("root", f"consolidated metadata names all {len(on_disk)} nodes")
    return root


def check_table(store, name: str, r: Report):
    where = f"tables/{name}"
    t = node(store, where)
    attrs = t.get("attributes", {}) if t else {}
    spatial = node(store, f"{where}/obsm/spatial")
    if not spatial:
        # The viewer takes one point per cell from the region shapes element
        # (a Point `<region>_for_vitessce` sibling first, else the polygons'
        # centres), matched by instance id (lib/spatialdata/centroids.ts)
        reg = attrs.get("region")
        names = reg if isinstance(reg, list) else [reg] if reg else []
        source = next(
            (
                f"shapes/{n}{suffix}"
                for n in names
                for suffix in (CONVENTIONS["shapes"]["pointSiblingSuffix"], "")
                if node(store, f"shapes/{n}{suffix}")
            ),
            None,
        )
        if source and attrs.get("instance_key"):
            r.ok(where, f"no obsm/spatial: the cells' positions come from {source}")
        else:
            r.fail(where, "no obsm/spatial, and no region shapes element with an instance_key to take positions from")
    elif len(spatial.get("shape", [])) != 2 or spatial["shape"][1] not in (2, 3):
        r.fail(where, f"obsm/spatial shape {spatial.get('shape')} is not n_obs × 2 or 3")
    else:
        r.ok(where, f"obsm/spatial {spatial['shape']}")
    obs = node(store, f"{where}/obs")
    index = (obs or {}).get("attributes", {}).get("_index", "_index")
    if not node(store, f"{where}/obs/{index}"):
        r.fail(where, f"no obs index column ({index})")
    var = node(store, f"{where}/var")
    vindex = (var or {}).get("attributes", {}).get("_index", "_index")
    if not node(store, f"{where}/var/{vindex}"):
        r.warn(where, f"no var index column ({vindex}): no gene names")
    x = node(store, f"{where}/X")
    layers = children(store, f"{where}/layers")
    if not x:
        r.warn(where, "no X: no gene expression")
    else:
        enc = x.get("attributes", {}).get("encoding-type", "array" if x.get("node_type") == "array" else "?")
        if enc.startswith("csc"):
            cs = chunk_shape(node(store, f"{where}/X/data") or {"chunk_grid": {"configuration": {"chunk_shape": [0]}}})
            r.ok(where, f"X is CSC (data chunks of {cs[0]})")
        elif enc.startswith("csr"):
            # The reader takes any CSC layer of X's shape (SpatialDataTableAdapter)
            shape = x.get("attributes", {}).get("shape")
            csc = [
                l for l in layers
                if str((node(store, f"{where}/layers/{l}") or {}).get("attributes", {}).get("encoding-type", "")).startswith("csc")
                and (node(store, f"{where}/layers/{l}") or {}).get("attributes", {}).get("shape") == shape
            ]
            if csc:
                r.ok(where, f"X is CSR with a CSC copy in layers/{csc[0]}")
            else:
                r.fail(where, "X is CSR with no CSC layer: every gene would scan the whole matrix")
        else:
            r.ok(where, f"X is {enc}")
    for layer in attrs.get("layers", []) if isinstance(attrs.get("layers"), list) else []:
        if layer not in layers:
            r.fail(where, f"attributes name a layer {layer!r} that does not exist")
    for key in ("instance_key", "region_key"):
        v = attrs.get(key)
        if v and not node(store, f"{where}/obs/{v}"):
            r.fail(where, f"{key} {v!r} names no obs column")
    if attrs.get("instance_key"):
        r.ok(where, f"instance_key {attrs['instance_key']!r}")
    else:
        r.warn(where, "no instance_key: polygons cannot be matched to cells")


def check_points(store, name: str, r: Report, sample: int):
    where = f"points/{name}"
    p = node(store, where)
    attrs = (p or {}).get("attributes", {})
    if "points:version" not in attrs:
        if parquet_meta(store, where, "points.parquet")[0]:
            r.warn(where, "parquet points (spatialdata's own layout): not read by the viewer; convert with scripts/points_zarr.py")
        else:
            r.warn(where, "not a multiscale points element")
        return
    errors = schema_errors(attrs, "points.schema.json")
    if errors:
        for e in errors[:5]:
            r.fail(where, f"metadata: {e}")
        return
    r.ok(where, f"metadata valid (points:version {attrs['points:version']}, {len(attrs['multiscales']['layout'])} levels)")
    if isinstance(store, RootStore):
        r.warn(where, "root only: level invariants not checked")
        return
    try:
        import points_zarr as pz

        report = pz.validate(store / where, sample=sample)
        r.ok(where, f"levels {', '.join(f'{k}: {v['n']}' for k, v in report.items())}")
    except Exception as e:  # noqa: BLE001
        r.fail(where, f"invariants: {e}")


def check_shapes(store, name: str, r: Report):
    where = f"shapes/{name}"
    s = node(store, where)
    attrs = (s or {}).get("attributes", {})
    if not attrs.get("coordinateTransformations"):
        r.fail(where, "no coordinateTransformations")
    try:
        kind, meta = parquet_meta(store, where, "shapes.parquet")
    except Exception as e:  # noqa: BLE001
        r.fail(where, f"unreadable: {e}")
        return
    if kind == "dir":
        r.fail(where, "shapes.parquet is a directory of parts: the viewer reads one file")
        return
    if kind is None:
        if isinstance(store, RootStore):
            r.warn(where, "root only: shapes.parquet not read")
        else:
            r.fail(where, "no shapes.parquet")
        return
    try:
        kv = {k.decode(): v.decode() for k, v in (meta.metadata or {}).items()}
        geo = json.loads(kv["geo"]) if "geo" in kv else None
        if not geo:
            r.fail(where, "no GeoParquet 'geo' metadata")
        else:
            r.ok(where, f"GeoParquet, geometry column {geo.get('primary_column')!r}, {meta.num_rows} rows, {meta.num_row_groups} row group(s)")
        pandas = json.loads(kv["pandas"]) if "pandas" in kv else None
        if not (pandas and pandas.get("index_columns")):
            r.warn(where, "no pandas index column: the first non-geometry column is taken as the instance id")
        if meta.num_row_groups == 1 and meta.num_rows > 100_000:
            r.warn(where, "one row group: a region cannot be read by range")
    except Exception as e:  # noqa: BLE001
        r.fail(where, f"unreadable: {e}")
    outlines = CONVENTIONS["shapes"]["cellOutlineElements"]
    sibling = CONVENTIONS["shapes"]["pointSiblingSuffix"]
    if name not in outlines and not name.endswith(sibling):
        r.warn(where, f"not a name the viewer draws as outlines ({', '.join(outlines)})")


def check_image(store, name: str, r: Report):
    where = f"images/{name}"
    g = node(store, where)
    ome = (g or {}).get("attributes", {}).get("ome", {})
    ms = (ome.get("multiscales") or [None])[0]
    if not ms:
        r.fail(where, "no ome.multiscales")
        return
    axes = [a.get("name") for a in ms.get("axes", [])]
    if axes != ["c", "y", "x"]:
        r.fail(where, f"axes {axes}, not c, y, x")
    for i, ds in enumerate(ms.get("datasets", [])):
        a = node(store, f"{where}/{ds.get('path')}")
        if not a:
            r.fail(where, f"level {ds.get('path')} missing")
            continue
        if a.get("data_type") not in ("uint8", "uint16"):
            r.fail(where, f"level {ds['path']} is {a.get('data_type')}, not uint8 or uint16")
        cs = chunk_shape(a)
        if max(cs[-2:]) > 1024:
            r.warn(where, f"level {ds['path']} chunks {cs}: over 1024² means {cs[-1] * cs[-2] * (2 if a.get('data_type') == 'uint16' else 1) / 1e6:.0f} MB a tile")
        if not ds.get("coordinateTransformations"):
            r.fail(where, f"level {ds['path']} has no coordinateTransformations")
        if i == 0:
            r.ok(where, f"{a.get('data_type')} {a.get('shape')}, chunks {cs}, {len(ms.get('datasets', []))} levels")


# ── coordinate systems (store-layout.md §8; lib/spatialdata/coordinates.ts) ─

ELEMENT_KINDS = ("images", "labels", "shapes", "points")


def output_name(t: dict) -> str | None:
    o = t.get("output")
    return o if isinstance(o, str) else (o or {}).get("name")


def element_transforms(kind: str, attrs: dict) -> tuple[list[dict], list[str]]:
    """An element's transforms and the names of its intrinsic axes."""
    if kind in ("images", "labels"):
        ms = ((attrs.get("ome") or attrs).get("multiscales") or [{}])[0]
        axes = [a if isinstance(a, str) else a.get("name") for a in ms.get("axes", [])]
        return ms.get("coordinateTransformations") or [], axes
    default = ["x", "y", "z"] if kind == "points" else ["x", "y"]
    return attrs.get("coordinateTransformations") or [], attrs.get("axes") or default


IDENTITY = (1.0, 0.0, 0.0, 1.0, 0.0, 0.0)


def compose(m, n):
    """First m, then n; maps are (a, b, c, d, tx, ty): x' = a x + b y + tx, y' = c x + d y + ty."""
    a, b, c, d, tx, ty = m
    e, f, g, h, ux, uy = n
    return (
        e * a + f * c, e * b + f * d,
        g * a + h * c, g * b + h * d,
        e * tx + f * ty + ux, g * tx + h * ty + uy,
    )


def invert(m):
    a, b, c, d, tx, ty = m
    det = a * d - b * c
    if not det:
        raise ValueError("the transform is not invertible")
    ia, ib, ic, id_ = d / det, -b / det, -c / det, a / det
    return ia, ib, ic, id_, -(ia * tx + ib * ty), -(ic * tx + id_ * ty)


def chain(transforms: list[dict], axes: list[str]) -> tuple[float, float, float, float, float, float]:
    """The 2-D affine (a, b, c, d, tx, ty) of scale / translation / identity /
    affine steps applied in order (a sequence counts as its steps), as
    lib/spatialdata/coordinates.ts computes it; ValueError on any other step."""
    m = IDENTITY
    for t in transforms:
        names = [a.get("name") for a in (t.get("input") or {}).get("axes", [])] if isinstance(t.get("input"), dict) else []
        names = names or axes

        def pick(v, fill):
            return [v[names.index(a)] if a in names and names.index(a) < len(v) else fill for a in ("x", "y")]

        kind = t.get("type")
        if kind == "identity":
            continue
        if kind == "sequence":
            m = compose(m, chain(t.get("transformations") or [], names))
        elif kind == "scale" and t.get("scale"):
            a, b = pick(t["scale"], 1.0)
            m = compose(m, (a, 0.0, 0.0, b, 0.0, 0.0))
        elif kind == "translation" and t.get("translation"):
            a, b = pick(t["translation"], 0.0)
            m = compose(m, (1.0, 0.0, 0.0, 1.0, a, b))
        elif kind == "affine" and t.get("affine") is not None:
            # One row per output axis, one column per input axis plus the
            # translation; a `c` axis carries the identity
            rows = t["affine"]
            outs = [a.get("name") for a in (t.get("output") or {}).get("axes", [])] if isinstance(t.get("output"), dict) else []
            outs = outs or names
            xi, yi = names.index("x") if "x" in names else -1, names.index("y") if "y" in names else -1
            xo, yo = outs.index("x") if "x" in outs else -1, outs.index("y") if "y" in outs else -1

            def at(r, c, fill):
                return rows[r][c] if r >= 0 and c >= 0 and r < len(rows) and c < len(rows[r]) else fill

            def last(r):
                return rows[r][-1] if 0 <= r < len(rows) and rows[r] else 0.0

            m = compose(
                m,
                (
                    at(xo, xi, 1.0), at(xo, yi, 0.0),
                    at(yo, xi, 0.0), at(yo, yi, 1.0),
                    last(xo), last(yo),
                ),
            )
        else:
            raise ValueError(f'unsupported transform "{kind}"')
    return m


def to_system(el: tuple[list[dict], list[str]], system: str | None):
    """The element's map into `system`, or the reason there is none."""
    transforms, axes = el
    if not any(output_name(t) for t in transforms):
        try:
            return chain(transforms, axes), None
        except ValueError as e:
            return None, str(e)
    entry = next((t for t in transforms if output_name(t) == system), None)
    if entry is None:
        return None, f'no transform into "{system}"'
    try:
        return chain([entry], axes), None
    except ValueError as e:
        return None, f'{e} into "{system}"'


def is_identity(m) -> bool:
    return all(abs(v - e) < 1e-9 for v, e in zip(m, IDENTITY))


def describe(m) -> str:
    """A map in words: its scale and offset, and the rotation if it has one."""
    import math

    a, b, c, d, tx, ty = m
    text = f"scale {math.hypot(a, c):g}, {math.hypot(b, d):g}; offset {tx:g}, {ty:g}"
    if abs(b) > 1e-9 or abs(c) > 1e-9:
        text += f"; rotated {math.degrees(math.atan2(c, a)):.1f}°"
    return text


def default_table(root_attrs: dict, tables: list[str]) -> str | None:
    """The curated default, else `table`, else the first by name
    (where each is named: docs/specs/conventions.json)."""
    named = CONVENTIONS["tables"]["curatedDefaultTable"]
    fallback = CONVENTIONS["tables"]["fallbackTableName"]
    curated = (root_attrs.get(named["attribute"]) or {}).get(named["key"])
    if curated in tables:
        return curated
    if fallback in tables:
        return fallback
    return sorted(tables)[0] if tables else None


def coordinate_systems(store) -> dict:
    """The target system the viewer chooses and each element's placement in
    the cells' coordinates, as lib/spatialdata/coordinates.ts computes them:
    {systems, target, rule, region, elements, place(path) -> (map, reason)}."""
    elements = {}
    for kind in ELEMENT_KINDS:
        for name in children(store, kind):
            node_ = node(store, f"{kind}/{name}") or {}
            elements[f"{kind}/{name}"] = element_transforms(kind, node_.get("attributes", {}))

    root_attrs = (node(store, "") or {}).get("attributes", {})
    table = default_table(root_attrs, children(store, "tables"))
    region = None
    if table:
        reg = (node(store, f"tables/{table}") or {}).get("attributes", {}).get("region")
        names = reg if isinstance(reg, list) else [reg] if reg else []
        for kind in ("shapes", "points"):
            hit = next((n for n in names if f"{kind}/{n}" in elements), None)
            if hit:
                region = f"{kind}/{hit}"
                break
    if region is None:
        first = next((p for p in elements if p.startswith("shapes/")), None) or next(
            (p for p in elements if p.startswith("points/")), None
        )
        region = first

    outs = {p: {output_name(t) for t in el[0]} - {None} for p, el in elements.items()}
    systems = sorted(set().union(*outs.values())) if outs else []
    every = set.intersection(*outs.values()) if outs else set()
    region_identity = (
        sorted(s for s in outs[region] if (m := to_system(elements[region], s)[0]) and is_identity(m)) if region else []
    )
    target, rule = None, "none"
    if "global" in every:
        target, rule = "global", "global"
    elif region_identity:
        target, rule = region_identity[0], "region identity"
    elif region is None and len(every) == 1:
        target, rule = next(iter(every)), "only common system"
    elif not systems:
        rule = "only common system"
    region_map, region_reason = to_system(elements[region], target) if region else (IDENTITY, None)

    def place(path: str):
        if path not in elements:
            return None, f"{path} is not in the store"
        if systems and target is None:
            return None, "the store's elements share no system"
        if region_map is None:
            return None, f"the cells: {region_reason}"
        m, reason = to_system(elements[path], target)
        if m is None:
            return None, reason
        return compose(m, invert(region_map)), None

    return {"systems": systems, "target": target, "rule": rule, "region": region, "elements": elements, "place": place}


def check_coordinates(store, r: Report):
    """store-layout.md §8: one target; images reach it (by a scale, a
    translation or an affine, which may rotate); points and shapes are in
    the cells' frame (they are drawn as written)."""
    cs = coordinate_systems(store)
    if not cs["elements"]:
        return
    if cs["systems"] and cs["target"] is None:
        r.fail("coordinates", f"no target system: the elements share none of {cs['systems']}")
        return
    r.ok("coordinates", f"target {cs['target']!r} ({cs['rule']}), the cells' coordinates are {cs['region']}'s")
    images = [p for p in cs["elements"] if p.startswith("images/")]
    placed = []
    for path in cs["elements"]:
        m, reason = cs["place"](path)
        kind = path.split("/")[0]
        if kind in ("points", "shapes"):
            if m is None:
                r.fail(path, f"not placed with the cells: {reason}")
            elif not is_identity(m):
                r.fail(path, f"declares a frame other than the cells' ({describe(m)}): it is drawn as written, so its transforms must match theirs")
        elif kind == "images":
            if m is None:
                r.warn(path, f"not drawn: {reason}")
            else:
                placed.append(path)
                if abs(m[1]) > 1e-9 or abs(m[2]) > 1e-9:
                    r.ok(path, f"placed by an affine ({describe(m)})")
    if images and not placed:
        r.fail("coordinates", "no image can be placed with the cells")


def check_mapping(path: Path, r: Report):
    try:
        mapping = read_json(path)
    except Exception as e:  # noqa: BLE001
        r.fail("mapping.json", f"unreadable: {e}")
        return
    errors = schema_errors(mapping, "mapping.schema.json")
    for e in errors:
        r.fail("mapping.json", e)
    if not errors:
        what = []
        if mapping.get("links"):
            what.append(f"links by {mapping.get('linkColumn')!r} ({mapping.get('source', 's3')})")
        if mapping.get("image"):
            what.append(f"image ({mapping.get('imageSource', 's3')})")
        r.ok("mapping.json", ", ".join(what) or "names the store only")


def validate(store, mapping: Path | None = None, sample: int = 50_000) -> Report:
    r = Report()
    if isinstance(store, RootStore):
        return _validate(store, r, None, sample)
    store = Path(store)
    if not store.is_dir():
        r.fail("store", f"{store} is not a directory")
        return r
    mapping = mapping or (store.parent / "mapping.json" if (store.parent / "mapping.json").is_file() else None)
    return _validate(store, r, mapping, sample)


def _validate(store, r: Report, mapping: Path | None, sample: int) -> Report:
    root = check_root(store, r)
    if root is None:
        return r
    for name in children(store, "tables"):
        check_table(store, name, r)
    for name in children(store, "points"):
        check_points(store, name, r, sample)
    for name in children(store, "shapes"):
        check_shapes(store, name, r)
    for name in children(store, "images"):
        check_image(store, name, r)
    check_coordinates(store, r)
    if not (children(store, "tables") or children(store, "points") or children(store, "images")):
        r.fail("store", "no table, points element or image to show")
    if mapping:
        check_mapping(mapping, r)
    return r


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("store", nargs="?")
    ap.add_argument("--root-json", type=Path, help="check a store from its root zarr.json alone")
    ap.add_argument("--mapping", type=Path)
    ap.add_argument("--sample", type=int, default=50_000, help="rows checked per level of a points element")
    ap.add_argument("--json", action="store_true", help="the report as JSON")
    args = ap.parse_args(argv)
    if args.root_json:
        r = validate(RootStore(read_json(args.root_json)))
    elif args.store:
        r = validate(Path(args.store), args.mapping, args.sample)
    else:
        ap.error("give a store or --root-json")
    if args.json:
        print(json.dumps([{"level": lv, "where": w, "message": m} for lv, w, m in r.lines], indent=1))
    else:
        for level, where, msg in r.lines:
            print(f"{level:4} {where}: {msg}")
    return 1 if r.failed else 0


if __name__ == "__main__":
    sys.exit(main())
