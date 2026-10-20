#!/usr/bin/env python3
"""Census of the SpatialData-DB stores (docs/plans/SDDB.md, step 1).

  python scripts/sddb_census.py [--out perf/sddb] [--workers 6] [--limit N]

Lists the latest SpatialData artifacts of the public Lamin instance
`scverse/spatialdata-db` through the hub REST, reads each store's root
`zarr.json` and its parquet footers through Lamin's storage proxy, and
writes, per store, what the viewer needs to know: elements, image layout,
shapes geometry and row counts, points layout, CSC layers, coordinate
systems and the target the viewer would choose, licence, and the verdict of
`validate_store.py --root-json`.

Writes `<out>/census.json` (one record per store), `<out>/CENSUS.md` and
`<out>/NO-CSC.md` (the stores whose default table has no CSC layer).
Everything fetched is cached under `<out>/.cache/`, so a re-run fetches
nothing; delete it to re-read the instance.

The proxy answers 403 to Python's default user agent, so every request
sends a browser one. Requests are throttled (a few in flight, backing off
on 429 and 5xx): the proxy's rate limit is not documented.

Needs pyarrow (the parquet footers) and the standard library.
"""
from __future__ import annotations

import argparse
import json
import re
import struct
import sys
import threading
import time
import urllib.error
import urllib.request
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import validate_store as vs  # noqa: E402

INSTANCE = "f95ff493-5f8b-4c1a-ad80-b781031e8faf"
API = f"https://aws.eu-central-1.lamin.ai/api/instances/{INSTANCE}/modules/core"
PROXY = "https://lamin.ai/storage/s3/scverse-spatial-eu-central-1//"
UA = "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/141.0 Safari/537.36"
LICENSES = {
    "CC BY-NC-SA 4.0": "noncommercial",
    "CC BY 4.0": "commercial",
    "Creative Commons Attribution 4.0 International": "commercial",
    "NCBI GEO policy": "geo",
}
SUPPORTED = {"identity", "scale", "translation", "sequence"}
TAIL = 65536

_slots = threading.Semaphore(6)


# ── fetching ────────────────────────────────────────────────────────────

def request(url: str, *, method="GET", body=None, headers=None, tries=5):
    """(status, headers, bytes); 404 is an answer, not an error."""
    h = {"User-Agent": UA, **(headers or {})}
    data = json.dumps(body).encode() if body is not None else None
    if data:
        h["Content-Type"] = "application/json"
    for attempt in range(tries):
        try:
            with _slots:
                with urllib.request.urlopen(urllib.request.Request(url, data=data, method=method, headers=h), timeout=120) as r:
                    return r.status, r.headers, r.read()
        except urllib.error.HTTPError as e:
            if e.code in (403, 404, 416):
                return e.code, e.headers, b""
            if e.code not in (429, 500, 502, 503, 504) or attempt == tries - 1:
                raise
        except (urllib.error.URLError, TimeoutError, ConnectionError):
            if attempt == tries - 1:
                raise
        time.sleep(2 ** attempt)
    raise RuntimeError("unreachable")


def list_artifacts() -> list[dict]:
    out, offset = [], 0
    body = {
        "filter": {"and": [{"suffix": {"eq": ".zarr"}}, {"is_latest": {"eq": True}}]},
        "select": ["uid", "key", "size", "description", "n_files", "created_at", "_key_is_virtual", "_overwrite_versions", "ulabels"],
    }
    while True:
        _, _, b = request(f"{API}/artifact?limit=200&limit_to_many=50&offset={offset}", method="POST", body=body)
        page = json.loads(b)
        for a in page:
            a["ulabels"] = [label["name"] for label in a.get("ulabels") or []]
        out += page
        offset += 200
        if len(page) < 200:
            return out


def store_path(a: dict) -> str:
    """The artifact's folder in the bucket (lamindb's auto storage key)."""
    if not a.get("_key_is_virtual"):
        return a["key"]
    uid = a["uid"][:16] if a.get("_overwrite_versions", True) else a["uid"]
    return f".lamindb/{uid}.zarr"


def cached(path: Path, fetch):
    if path.is_file():
        return json.loads(path.read_text())
    value = fetch()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value))
    return value


def footer(url: str) -> dict | None:
    """A parquet file's footer, from a suffix range: rows, groups, columns, geo."""
    import pyarrow as pa
    import pyarrow.parquet as pq

    status, _, b = request(url, headers={"Range": f"bytes=-{TAIL}"})
    if status == 404 or not b:
        return None
    if b[-4:] != b"PAR1":
        return {"error": "not parquet"}
    n = struct.unpack("<I", b[-8:-4])[0]
    if n + 8 > len(b):
        _, _, b = request(url, headers={"Range": f"bytes=-{n + 8}"})
    meta = pq.read_metadata(pa.BufferReader(b"PAR1" + b[-8 - n:]))
    kv = {k.decode(): v.decode() for k, v in (meta.metadata or {}).items()}
    geo = json.loads(kv["geo"]) if "geo" in kv else None
    pandas = json.loads(kv["pandas"]) if "pandas" in kv else None
    return {
        "rows": meta.num_rows,
        "row_groups": meta.num_row_groups,
        "columns": meta.schema.names,
        "geometry_types": sorted({t for c in (geo or {}).get("columns", {}).values() for t in c.get("geometry_types") or []}),
        "geo": geo is not None,
        "index_columns": [c for c in (pandas or {}).get("index_columns", []) if isinstance(c, str)],
    }


def count_parts(base: str) -> int:
    """Parts `part.0 … part.N-1` of a parquet directory, found by probing."""
    def exists(i):
        return request(f"{base}/part.{i}.parquet", method="HEAD")[0] == 200

    if not exists(0):
        return 0
    hi = 1
    while exists(hi):
        hi *= 2
    lo = hi // 2  # exists(lo), not exists(hi)
    while hi - lo > 1:
        mid = (lo + hi) // 2
        lo, hi = (mid, hi) if exists(mid) else (lo, mid)
    return hi


def fetch_store(a: dict, cache: Path) -> dict:
    base = PROXY + store_path(a)
    uid = a["uid"]

    def root():
        status, _, b = request(f"{base}/zarr.json")
        return {"status": status, "root": json.loads(b) if status == 200 else None}

    got = cached(cache / "roots" / f"{uid}.json", root)
    if not got["root"]:
        return {"root": None, "status": got["status"], "parquet": {}}
    nodes = (got["root"].get("consolidated_metadata") or {}).get("metadata") or {}

    def parquet():
        out = {}
        for name in vs.children(vs.RootStore(got["root"]), "shapes"):
            out[f"shapes/{name}"] = {"file": footer(f"{base}/shapes/{name}/shapes.parquet")}
        for name in vs.children(vs.RootStore(got["root"]), "points"):
            if "points:version" in nodes.get(f"points/{name}", {}).get("attributes", {}):
                continue
            d = f"{base}/points/{name}/points.parquet"
            part0 = footer(f"{d}/part.0.parquet")
            if part0:
                out[f"points/{name}"] = {"dir": part0, "parts": count_parts(d)}
            else:
                out[f"points/{name}"] = {"file": footer(d)}
        return out

    return {"root": got["root"], "status": 200, "parquet": cached(cache / "parquet" / f"{uid}.json", parquet)}


# ── reading a root ──────────────────────────────────────────────────────

def transforms_of(attrs: dict) -> list[dict]:
    """An element's coordinateTransformations: images keep them in the multiscale."""
    ms = (attrs.get("ome") or attrs).get("multiscales")
    if ms:
        return ms[0].get("coordinateTransformations") or []
    return attrs.get("coordinateTransformations") or []


def types_in(ts: list[dict]) -> set[str]:
    out = set()
    for t in ts:
        out.add(t.get("type"))
        out |= types_in(t.get("transformations") or [])
    return out


def output(t: dict) -> str | None:
    o = t.get("output")
    return o.get("name") if isinstance(o, dict) else o


FOLDERS = {"visium": "Visium", "visiumhd": "Visium HD", "xenium": "Xenium", "atera": "Atera", "spatialtranscriptomics": "spatial transcriptomics"}


def technology(key: str | None, assay: str | None) -> str:
    """From the artifact's folder (consistent), else the curated Assay."""
    folder = (key or "").split("/")[0] if "/" in (key or "") else ""
    return FOLDERS.get(folder) or (f"{folder} (folder)" if folder else f"{assay or '?'} (top level)")


def chunk_of(meta: dict) -> list[int]:
    return vs.chunk_shape(meta) if meta.get("node_type") == "array" else []


def survey(a: dict, got: dict) -> dict:
    labels = a["ulabels"]
    desc = a.get("description") or ""
    tier = re.search(r"tier (\w+)", desc)
    rec = {
        "uid": a["uid"],
        "key": a["key"],
        "size": a["size"],
        "n_files": a["n_files"],
        "tier": tier.group(1) if tier else None,
        "license": next((label for label in labels if label in LICENSES), None),
    }
    rec["license_class"] = LICENSES.get(rec["license"], "unknown")
    root = got["root"]
    if not root:
        rec["status"] = f"root {got['status']}"
        return rec
    attrs = root.get("attributes", {})
    nodes = (root.get("consolidated_metadata") or {}).get("metadata") or {}
    general = attrs.get("metadata_general") or {}
    rs = vs.RootStore(root)
    rec.update(
        status="ok" if "spatialdata_attrs" in attrs else "not spatialdata",
        sddb_id=(attrs.get("sddb_provenance") or {}).get("sddb_id"),
        technology=technology(a["key"], general.get("Assay")),
        spatialdata=(attrs.get("spatialdata_attrs") or {}).get("version"),
        reader=attrs.get("spatialdata_io_reader"),
        default_table=general.get("Default Table"),
        organism=general.get("Organism"),
        tissue=general.get("Tissue"),
    )

    # elements and their transforms
    elements: dict[str, list[dict]] = {}
    kinds = ("images", "labels", "shapes", "points")
    for kind in kinds:
        for name in vs.children(rs, kind):
            elements[f"{kind}/{name}"] = transforms_of(nodes.get(f"{kind}/{name}", {}).get("attributes", {}))

    rec["images"] = []
    for name in vs.children(rs, "images"):
        g = nodes.get(f"images/{name}", {}).get("attributes", {})
        ms = ((g.get("ome") or g).get("multiscales") or [{}])[0]
        datasets = ms.get("datasets") or []
        level0 = nodes.get(f"images/{name}/{datasets[0]['path']}", {}) if datasets else {}
        rec["images"].append({
            "name": name,
            "axes": [x.get("name") for x in ms.get("axes", [])],
            "dtype": level0.get("data_type"),
            "shape": level0.get("shape"),
            "chunk": chunk_of(level0),
            "levels": len(datasets),
        })
    rec["labels"] = []
    for name in vs.children(rs, "labels"):
        g = nodes.get(f"labels/{name}", {}).get("attributes", {})
        ms = ((g.get("ome") or g).get("multiscales") or [{}])[0]
        datasets = ms.get("datasets") or []
        level0 = nodes.get(f"labels/{name}/{datasets[0]['path']}", {}) if datasets else {}
        rec["labels"].append({"name": name, "dtype": level0.get("data_type"), "shape": level0.get("shape"), "chunk": chunk_of(level0), "levels": len(datasets)})

    rec["shapes"] = []
    for name in vs.children(rs, "shapes"):
        f = (got["parquet"].get(f"shapes/{name}") or {}).get("file")
        geom = (f or {}).get("geometry_types", [])
        cols = (f or {}).get("columns", [])
        if not f:
            kind = "unread"
        elif "radius" in cols and geom in (["Point"], []):
            kind = "circle"
        elif geom and set(geom) <= {"Polygon", "MultiPolygon"}:
            kind = "square" if re.search(r"square_\d+um|_bins?_", name) else "polygon"
        else:
            kind = "/".join(geom) or "untyped"
        rec["shapes"].append({"name": name, "kind": kind, "rows": (f or {}).get("rows"), "row_groups": (f or {}).get("row_groups"), "index": (f or {}).get("index_columns")})

    rec["points"] = []
    for name in vs.children(rs, "points"):
        pattrs = nodes.get(f"points/{name}", {}).get("attributes", {})
        if "points:version" in pattrs:
            rec["points"].append({"name": name, "layout": "points-zarr", "rows": pattrs.get("points:n_points")})
            continue
        p = got["parquet"].get(f"points/{name}") or {}
        if p.get("dir"):
            rec["points"].append({"name": name, "layout": "parquet-dir", "parts": p["parts"], "part0_rows": p["dir"]["rows"], "columns": p["dir"]["columns"]})
        elif p.get("file"):
            rec["points"].append({"name": name, "layout": "parquet-file", "rows": p["file"]["rows"], "columns": p["file"]["columns"]})
        else:
            rec["points"].append({"name": name, "layout": "unread"})

    rec["tables"] = []
    for name in vs.children(rs, "tables"):
        t = f"tables/{name}"
        ta = nodes.get(t, {}).get("attributes", {})
        x = nodes.get(f"{t}/X", {})
        xa = x.get("attributes", {})
        xenc = xa.get("encoding-type") or ("dense" if x.get("node_type") == "array" else None)
        xshape = xa.get("shape") or x.get("shape")
        csc = [
            layer for layer in vs.children(rs, f"{t}/layers")
            if str(nodes.get(f"{t}/layers/{layer}", {}).get("attributes", {}).get("encoding-type", "")).startswith("csc")
            and nodes.get(f"{t}/layers/{layer}", {}).get("attributes", {}).get("shape") == xshape
        ]
        region = ta.get("region")
        spatial = nodes.get(f"{t}/obsm/spatial", {})
        rec["tables"].append({
            "name": name,
            "x": xenc,
            "shape": xshape,
            "csc": csc[0] if csc else ("X" if str(xenc).startswith("csc") else None),
            "region": region if isinstance(region, list) else [region] if region else [],
            "spatial": spatial.get("data_type"),
        })

    # coordinate systems: the validator's rule, which the viewer shares
    # (validate_store.coordinate_systems ↔ lib/spatialdata/coordinates.ts)
    cs = vs.coordinate_systems(rs)
    target = cs["target"]
    types = set().union(*(types_in(ts) for ts in elements.values())) if elements else set()
    placed = {p for p in elements if cs["place"](p)[0] is not None}
    images = [p for p in placed if p.startswith("images/")]
    # A full image is its target's own pixels (identity); others are scaled
    full = [p for p in images if (m := vs.to_system(cs["elements"][p], target)[0]) and vs.is_identity(m)]
    tables = {t["name"]: t for t in rec["tables"]}
    default = tables.get(vs.default_table(attrs, list(tables)) or "")
    regions = [f"{k}/{r}" for r in (default or {}).get("region", []) for k in ("shapes", "points") if f"{k}/{r}" in elements]
    region_systems = [
        frozenset(s for s in {output(t) for t in elements[r]} - {None}
                  if (m := vs.to_system(cs["elements"][r], s)[0]) and vs.is_identity(m))
        for r in regions
    ]
    rec["coords"] = {
        "systems": cs["systems"],
        "has_global": "global" in cs["systems"],
        "target": target,
        "rule": cs["rule"],
        "region": cs["region"],
        "unreachable": sorted(set(elements) - placed),
        "full_images": sorted(full),
        "downscaled_images": sorted(set(images) - set(full)),
        "samples": len(set(region_systems)),
        "regions": len(regions),
        "transform_types": sorted(t for t in types if t),
        "unsupported": sorted(t for t in types if t and t not in SUPPORTED),
        "multi_entry": any(len(ts) > 1 for ts in elements.values()),
    }

    # the validator, from the root and the footers
    parquet = {}
    for where, p in got["parquet"].items():
        if p.get("file"):
            parquet[where] = ("file", _FakeMeta(p["file"]))
        elif p.get("dir"):
            parquet[where] = ("dir", None)
    report = vs.validate(vs.RootStore(root, parquet))
    rec["validator"] = {
        "ok": sum(1 for lv, _, _ in report.lines if lv == "OK"),
        "warn": sum(1 for lv, _, _ in report.lines if lv == "WARN"),
        "fail": [f"{w}: {m}" for lv, w, m in report.lines if lv == "FAIL"],
    }
    rec["profile"] = profile(rec)
    return rec


class _FakeMeta:
    """Enough of pyarrow's FileMetaData for validate_store.check_shapes."""

    def __init__(self, f: dict):
        self.num_rows = f["rows"]
        self.num_row_groups = f["row_groups"]
        kv = {}
        if f.get("geo"):
            kv[b"geo"] = json.dumps({"primary_column": "geometry"}).encode()
        if f.get("index_columns"):
            kv[b"pandas"] = json.dumps({"index_columns": f["index_columns"]}).encode()
        self.metadata = kv


def profile(rec: dict) -> str:
    """The class a store falls in for the viewer; never "unknown" when read."""
    c = rec["coords"]
    if c["target"] is None:
        place = f"no target ({c['rule']})"
    elif c["target"] == "global":
        place = "global"
    else:
        place = "named system"
    if c["samples"] > 1:
        place += f", {c['samples']} samples"
    image = "full image" if c["full_images"] else "downscaled image only" if c["downscaled_images"] else "no image"
    kinds = sorted({s["kind"] for s in rec["shapes"]})
    return f"{rec['technology']} · {place} · {image} · shapes: {'/'.join(kinds) or 'none'}"


# ── output ──────────────────────────────────────────────────────────────

def summary(records: list[dict], website_total: int) -> str:
    read = [r for r in records if r.get("status") == "ok"]
    sddb = [r for r in read if r.get("sddb_id")]
    lines = [
        "# SpatialData-DB census",
        "",
        f"Generated by `scripts/sddb_census.py` on {time.strftime('%Y-%m-%d')} from the hub REST and the storage",
        "proxy (see docs/plans/SDDB.md, step 1). Per-store records: `census.json`.",
        "",
        f"- Latest SpatialData artifacts in the instance: **{len(records)}**; root read: **{len(read)}**;",
        f"  with `sddb_provenance` (converted by SDDB): **{len(sddb)}** (the website lists {website_total}).",
        f"- Not read: {dict(Counter(r.get('status') for r in records if r.get('status') != 'ok'))}",
        "",
    ]

    def table(title, counter, unit="Stores"):
        lines.extend([f"## {title}", "", f"| | {unit} |", "|---|---|"])
        for k, n in counter.most_common():
            lines.append(f"| {k} | {n} |")
        lines.append("")

    table("Technology", Counter(r["technology"] for r in read))
    table("Tier", Counter(r["tier"] or "(none in description)" for r in read))
    table("Licence", Counter(f"{r['license'] or '(none)'} → {r['license_class']}" for r in read))
    table("Profile (technology · placement · image · shapes kinds)", Counter(r["profile"] for r in read))
    table("Target coordinate system rule", Counter(r["coords"]["rule"] for r in read))
    table("Samples (distinct systems of the default table's regions)", Counter(r["coords"]["samples"] for r in read))
    table("Transform types used", Counter(t for r in read for t in r["coords"]["transform_types"]))
    table("Stores using an unsupported transform", Counter("/".join(r["coords"]["unsupported"]) for r in read if r["coords"]["unsupported"]))
    table("Stores with elements declaring several transforms", Counter(r["coords"]["multi_entry"] for r in read))
    table("Images: dtype, max chunk side", Counter(f"{i['dtype']}, {max(i['chunk'][-2:]) if i['chunk'] else '?'}" for r in read for i in r["images"]), "Images")
    table("Shapes kinds", Counter(s["kind"] for r in read for s in r["shapes"]), "Elements")

    def bucket(n):
        if n is None:
            return "unread"
        for edge in (1_000, 10_000, 60_000, 250_000, 1_000_000):
            if n <= edge:
                return f"≤ {edge:,}"
        return "> 1,000,000"

    table("Circle / square rows per element (against the 60,000 per-view cap)", Counter(f"{s['kind']} {bucket(s['rows'])}" for r in read for s in r["shapes"] if s["kind"] in ("circle", "square")), "Elements")
    table("Points layout", Counter(p["layout"] for r in read for p in r["points"]), "Elements")
    default_csc = Counter()
    for r in read:
        t = next((t for t in r["tables"] if t["name"] == r["default_table"]), None)
        default_csc["no default table named" if t is None else "default has CSC" if t["csc"] else "default has no CSC"] += 1
    table("CSC on the default table", default_csc)
    table("CSC per table", Counter(("CSC" if t["csc"] else f"no CSC ({t['x']})") for r in read for t in r["tables"]), "Tables")
    table("The first table by name is not the default (the table profile reads the first)", Counter(
        r["technology"] for r in read if r["tables"] and r["default_table"] and sorted(t["name"] for t in r["tables"])[0] != r["default_table"]))
    table("Validator: stores with a FAIL", Counter(bool(r["validator"]["fail"]) for r in read))
    fails = Counter(re.sub(r"^[^:]+/[^:]+: ", "", f) for r in read for f in r["validator"]["fail"])
    table("Validator FAIL messages", fails, "Failures")
    return "\n".join(lines)


def no_csc(records: list[dict]) -> str:
    """The stores whose default table cannot be read one gene at a time."""
    rows = []
    for r in records:
        if r.get("status") != "ok":
            continue
        t = next((t for t in r["tables"] if t["name"] == r["default_table"]), None)
        if t is None or t["csc"]:
            continue
        with_csc = [x["name"] for x in r["tables"] if x["csc"]]
        rows.append((r["technology"], r.get("organism") or "", r.get("tissue") or "", r["uid"], t["name"], t["shape"], with_csc))
    rows.sort()
    lines = [
        "# Stores whose default table has no CSC layer",
        "",
        "Every gene read from these scans the whole matrix (CSR `X`, no",
        "`layers/csc_counts`). Generated by `scripts/sddb_census.py`; the last",
        "column names the tables in the same store that do have one.",
        "",
        f"**{len(rows)} stores.**",
        "",
        "| Technology | Organism | Tissue | uid | Default table | cells × genes | Tables with CSC |",
        "|---|---|---|---|---|---|---|",
    ]
    for tech, org, tissue, uid, name, shape, with_csc in rows:
        dims = " × ".join(f"{n:,}" for n in shape) if shape else "?"
        lines.append(f"| {tech} | {org} | {tissue} | `{uid}` | `{name}` | {dims} | {', '.join(f'`{c}`' for c in with_csc) or 'none'} |")
    return "\n".join(lines) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", type=Path, default=Path("perf/sddb"))
    ap.add_argument("--workers", type=int, default=6, help="requests in flight")
    ap.add_argument("--limit", type=int, help="only the first N artifacts (a trial run)")
    ap.add_argument("--website-total", type=int, default=1778)
    args = ap.parse_args(argv)
    global _slots
    _slots = threading.Semaphore(args.workers)
    cache = args.out / ".cache"
    t0 = time.time()
    artifacts = cached(cache / "artifacts.json", list_artifacts)
    artifacts = sorted(artifacts, key=lambda a: a["uid"])[: args.limit]
    print(f"{len(artifacts)} artifacts", flush=True)

    done = 0

    def one(a):
        nonlocal done
        try:
            rec = survey(a, fetch_store(a, cache))
        except Exception as e:  # noqa: BLE001
            rec = {"uid": a["uid"], "key": a["key"], "status": f"error: {type(e).__name__}: {e}"}
        done += 1
        if done % 50 == 0:
            print(f"{done}/{len(artifacts)} in {time.time() - t0:.0f}s", flush=True)
        return rec

    with ThreadPoolExecutor(args.workers * 2) as pool:
        records = list(pool.map(one, artifacts))

    args.out.mkdir(parents=True, exist_ok=True)
    with open(args.out / "census.json", "w") as f:
        f.write("[\n" + ",\n".join(json.dumps(r, separators=(",", ":")) for r in records) + "\n]\n")
    (args.out / "CENSUS.md").write_text(summary(records, args.website_total) + "\n")
    (args.out / "NO-CSC.md").write_text(no_csc(records))
    bad = [r for r in records if r.get("status") != "ok"]
    print(f"{len(records)} stores in {time.time() - t0:.0f}s; {len(bad)} not read", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
