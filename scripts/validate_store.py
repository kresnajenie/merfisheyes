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
        r.fail(where, "no obsm/spatial")
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
    if name not in ("cell_boundaries", "cell_polygons", "cells"):
        r.warn(where, "not a name the viewer looks for (cell_boundaries, cell_polygons, cells)")


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
