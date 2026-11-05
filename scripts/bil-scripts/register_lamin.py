#!/usr/bin/env python3
"""Register datasets already in S3 with the Lamin catalog, in place, and check them.

    python register_lamin.py ace-dip-vat                          # a SpatialData store
    python register_lamin.py --format legacy ace-can-ill ace-dip-tub
    python register_lamin.py --format legacy --catalog bil-completed.json ace-can-ill
    python register_lamin.py --check-only ace-dip-vat             # re-verify, write nothing

Nothing is uploaded, copied or renamed. The bytes were put in S3 by `aws s3 sync`
(s3_sync_store.sbatch, or the legacy s3_sync_*.sbatch); this script records each
dataset at the path it already has, so the URL the viewer opens and the path the
Lamin record names are the same files.

Layout, per format (both under s3://<bucket>/<prefix>/<sample>/):

  spatialdata  <name>.zarr  stats.json  mapping.json        prefix bil-psc-zarr
  legacy       meyes_output/            (cells, chunked)    prefix bil-psc-data2
               sm_output/<sample_id>/   (molecules, one per slice)

Never Artifact.from_spatialdata() and never .load() on these stores: both go
through spatialdata, and spatialdata 0.8 cannot read our points element.

Each artifact is checked after it is saved (and only checked with --check-only):

  record   the Lamin path is the S3 path, the key is not virtual
  listing  a fresh S3 listing has the file count and bytes the record holds
  public   the entry file answers on the public HTTPS URL, with CORS for the app
  lamin    the entry file reads back through Artifact.get(key=...).path

Exit 0 when every check passed, 1 when any failed, 2 on bad arguments.
Needs: lamindb (`lamin login` once), and AWS credentials that can list the
prefix and write its `.lamindb/` marker the first time a prefix is registered.
"""
from __future__ import annotations

import argparse
import gzip
import json
import sys
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field

import lamindb as ln

INSTANCE = "merfisheyes/eyes-catalog"
BUCKET = "merfisheyes-bil"   # defaults; --bucket / --region override (env.sh passes them)
REGION = "us-west-2"
PREFIX = {"spatialdata": "bil-psc-zarr", "legacy": "bil-psc-data2"}
APP = "https://www.merfisheyes.com"
CORS_ORIGIN = APP
PROJECT = "BIL"

# name -> dtype. Kept flat and free-text for now; species/tissue can move to
# bionty ontologies once the catalog settles.
FEATURES = {
    "bil_code": str,
    "data_format": str,      # spatialdata | merfisheyes-chunked
    "dataset_kind": str,     # single_cell | single_molecule | sidecar
    "species": str,
    "tissue": str,
    "platform": str,
    "investigator": str,
    "n_cells": int,
    "n_genes": int,
    "n_molecules": int,
    "pipeline_commit": str,
    "viewer_url": str,
}


@dataclass
class Item:
    """One thing to register: a folder or a file under the storage root."""
    rel: str                 # path below s3://bucket/prefix/
    kind: str                # single_cell | single_molecule | sidecar
    entry: str | None        # file inside a folder to read back; None for a file
    viewer: str | None = None
    values: dict = field(default_factory=dict)
    description: str = ""


@dataclass
class Result:
    sample: str
    key: str
    checks: dict = field(default_factory=dict)
    uid: str = ""
    n_files: int | None = None
    size: int | None = None

    @property
    def ok(self) -> bool:
        return bool(self.checks) and all(v == "ok" for v in self.checks.values())


def https(prefix: str, rel: str) -> str:
    return f"https://{BUCKET}.s3.{REGION}.amazonaws.com/{prefix}/{rel}"


def fetch(url: str, origin: str | None = None) -> tuple[bytes, dict]:
    req = urllib.request.Request(url, headers={"Origin": origin} if origin else {})
    with urllib.request.urlopen(req, timeout=60) as r:
        return r.read(), {k.lower(): v for k, v in r.headers.items()}


def read_json(url: str) -> dict | None:
    try:
        body, _ = fetch(url)
    except urllib.error.HTTPError:
        return None
    return json.loads(gzip.decompress(body) if url.endswith(".gz") else body)


def viewer_url(viewer: str, url: str) -> str:
    return f"{APP}/{viewer}/from-s3?url={urllib.parse.quote(url, safe='')}"


# ── what to register, per format ─────────────────────────────

def plan_spatialdata(sample: str, prefix: str, meta: dict) -> list[Item]:
    stats = read_json(https(prefix, f"{sample}/stats.json"))
    if stats is None:
        raise SystemExit(f"{sample}: no stats.json at {https(prefix, sample + '/stats.json')}")
    store = stats["store"]
    values = {
        "data_format": "spatialdata",
        "n_cells": stats.get("cells"),
        "n_genes": stats.get("genes"),
        "n_molecules": stats.get("molecules"),
        "pipeline_commit": stats.get("pipeline"),
    }
    desc = (f"{meta.get('title') or sample} — SpatialData store, {stats['cells']:,} cells, "
            f"{stats['genes']} genes, {stats.get('molecules') or 0:,} molecules")
    items = [Item(f"{sample}/{store}", "single_cell", "zarr.json", "spatialdata-viewer", values, desc)]
    for name in ("stats.json", "mapping.json"):
        items.append(Item(f"{sample}/{name}", "sidecar", None, None,
                          {"data_format": "spatialdata"}, f"{sample} {name}"))
    return items


def plan_legacy(sample: str, prefix: str, meta: dict) -> list[Item]:
    base = https(prefix, "")
    entries = meta.get("entries")
    if entries:   # the catalog knows exactly which folders the viewer opens
        rels = [(e["datasetType"], e["s3BaseUrl"].removeprefix(base).strip("/"), e.get("label"))
                for e in entries if e["s3BaseUrl"].startswith(base)]
    else:         # no catalog: meyes_output, plus every sm_output/<id>/ the bucket lists
        rels = [("single_cell", f"{sample}/meyes_output", "All Single Cell")]
        root = ln.UPath(f"s3://{BUCKET}/{prefix}/{sample}/sm_output")
        if root.exists():
            rels += [("single_molecule", f"{sample}/sm_output/{p.name}", p.name)
                     for p in root.iterdir() if p.is_dir()]
    items = []
    for kind, rel, label in rels:
        if kind == "single_cell":
            m = read_json(https(prefix, f"{rel}/manifest.json")) or {}
            s = m.get("statistics", {})
            values = {"n_cells": s.get("total_cells"), "n_genes": s.get("total_genes")}
            entry, viewer = "manifest.json", "viewer"
            what = f"cells, {s.get('total_cells') or 0:,} × {s.get('total_genes') or 0} genes"
        else:
            m = read_json(https(prefix, f"{rel}/manifest.json.gz")) or {}
            s = m.get("statistics", {})
            values = {"n_molecules": s.get("total_molecules"), "n_genes": s.get("unique_genes")}
            entry, viewer = "manifest.json.gz", "sm-viewer"
            what = f"molecules, {s.get('total_molecules') or 0:,} across {s.get('unique_genes') or 0} genes"
        values["data_format"] = "merfisheyes-chunked"
        desc = f"{meta.get('title') or sample} — {label}: legacy chunked {what}"
        items.append(Item(rel, kind, entry, viewer, values, desc))
    return items


# ── Lamin side ───────────────────────────────────────────────

def ensure_storage(root: str, dry: bool) -> None:
    if ln.Storage.filter(root=root).exists():
        return
    print(f"  storage  {root} is not registered yet" + (" — would register" if dry else " — registering"))
    if not dry:
        ln.Storage(root=root, description=f"MERFISHEYES datasets in {root}, read by URL by the viewer").save()


def ensure_features() -> None:
    for name, dtype in FEATURES.items():
        if not ln.Feature.filter(name=name).exists():
            ln.Feature(name=name, dtype=dtype).save()


def ensure_project():
    return ln.Project.filter(name=PROJECT).one_or_none() or ln.Project(name=PROJECT).save()


def register(item: Item, root: str, prefix: str, sample: str, meta: dict, project) -> ln.Artifact:
    af = ln.Artifact(f"{root}/{item.rel}", description=item.description).save()
    if af.description != item.description:   # same bytes already registered: refresh the text
        af.description = item.description
        af.save()
    values = {
        "bil_code": sample,
        "dataset_kind": item.kind,
        "species": meta.get("species"),
        "tissue": meta.get("tissue"),
        "platform": meta.get("platform"),
        "investigator": (meta.get("metadata") or {}).get("investigator"),
        **item.values,
    }
    if item.viewer:
        # a legacy cell folder is opened by its URL with a trailing slash, as the catalog writes it
        url = https(prefix, item.rel) + ("/" if item.entry == "manifest.json" else "")
        values["viewer_url"] = viewer_url(item.viewer, url)
    values = {k: v for k, v in values.items() if v not in (None, "")}
    af.features.set_values(values)   # replaces, so a re-run does not duplicate
    af.projects.add(project)
    return af


def check(item: Item, root: str, prefix: str, res: Result) -> None:
    """Fill res.checks; never raises."""
    key = item.rel
    try:
        af = ln.Artifact.filter(key=key, storage__root=root, is_latest=True).one_or_none()
    except Exception as e:
        res.checks["record"] = f"query failed: {e}"
        return
    if af is None:
        res.checks["record"] = "no record"
        return
    res.uid, res.n_files, res.size = af.uid, af.n_files, af.size

    want = f"{root}/{key}"
    res.checks["record"] = ("ok" if af.path.as_posix().rstrip("/") == want and not af._key_is_virtual
                            else f"path {af.path.as_posix()} virtual={af._key_is_virtual}")

    try:
        p = ln.UPath(want)
        if item.entry:   # a folder
            listing = p.fs.find(p.path, detail=True)
            n, size = len(listing), sum(v["size"] for v in listing.values())
        else:
            n, size = None, p.fs.info(p.path)["size"]
        res.checks["listing"] = ("ok" if (n == af.n_files and size == af.size)
                                 else f"S3 {n} files/{size} B, record {af.n_files}/{af.size}")
    except Exception as e:
        res.checks["listing"] = f"failed: {e}"

    entry_rel = f"{key}/{item.entry}" if item.entry else key
    try:
        body, headers = fetch(https(prefix, entry_rel), origin=CORS_ORIGIN)
        cors = headers.get("access-control-allow-origin") in (CORS_ORIGIN, "*")
        consolidated = True
        if item.entry == "zarr.json":
            consolidated = "consolidated_metadata" in json.loads(body)
        res.checks["public"] = ("ok" if cors and consolidated
                                else f"cors={cors} consolidated={consolidated}")
    except Exception as e:
        res.checks["public"] = f"failed: {e}"

    try:
        got = ln.Artifact.get(key=key, storage__root=root, is_latest=True).path
        data = (got / item.entry).read_bytes() if item.entry else got.read_bytes()
        res.checks["lamin"] = "ok" if data else "empty"
    except Exception as e:
        res.checks["lamin"] = f"failed: {e}"


def main() -> int:
    global BUCKET, REGION
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("samples", nargs="+", help="sample names (BIL codes), the folders under the prefix")
    ap.add_argument("--format", choices=sorted(PREFIX), default="spatialdata")
    ap.add_argument("--prefix", help="folder in the bucket (default: bil-psc-zarr / bil-psc-data2 by format)")
    ap.add_argument("--catalog", help="bil-completed.json: titles, species, tissue, and for legacy the entries")
    ap.add_argument("--bucket", default=BUCKET)
    ap.add_argument("--region", default=REGION)
    ap.add_argument("--instance", default=INSTANCE)
    ap.add_argument("--dry-run", action="store_true", help="show what would be registered; write nothing")
    ap.add_argument("--check-only", action="store_true", help="verify existing records; write nothing")
    args = ap.parse_args()
    BUCKET, REGION = args.bucket, args.region

    prefix = (args.prefix or PREFIX[args.format]).strip("/")
    root = f"s3://{BUCKET}/{prefix}"
    catalog = {}
    if args.catalog:
        catalog = {d["bilCode"]: d for d in json.load(open(args.catalog))}
        missing = [s for s in args.samples if s not in catalog]
        if missing:
            print(f"not in {args.catalog}: {', '.join(missing)}", file=sys.stderr)
            return 2

    ln.connect(args.instance)
    ln.settings.verbosity = "error"   # this script prints its own report
    print(f"instance {args.instance} as {ln.setup.settings.user.handle} | root {root} | format {args.format}")
    writing = not (args.dry_run or args.check_only)
    ensure_storage(root, dry=not writing)
    if writing:
        ensure_features()
        project = ensure_project()

    results: list[Result] = []
    for sample in args.samples:
        meta = catalog.get(sample, {})
        plan = (plan_spatialdata if args.format == "spatialdata" else plan_legacy)(sample, prefix, meta)
        print(f"\n{sample}: {len(plan)} artifact(s)")
        for item in plan:
            print(f"  {item.kind:15} {item.rel}")
            res = Result(sample, item.rel)
            if args.dry_run:
                print(f"    would register: {item.description}")
                continue
            if writing:
                try:
                    register(item, root, prefix, sample, meta, project)
                except Exception as e:
                    res.checks["register"] = f"failed: {type(e).__name__}: {e}"
                    results.append(res)
                    continue
            check(item, root, prefix, res)
            results.append(res)
            print(f"    {res.uid}  {res.n_files or 1} file(s)  {res.size or 0:,} B  "
                  + "  ".join(f"{k}={'ok' if v == 'ok' else 'FAIL'}" for k, v in res.checks.items()))
            for k, v in res.checks.items():
                if v != "ok":
                    print(f"      {k}: {v}")

    if args.dry_run:
        return 0
    bad = [r for r in results if not r.ok]
    print(f"\n{len(results) - len(bad)}/{len(results)} artifacts passed every check")
    for r in bad:
        print(f"  FAIL {r.key}: " + "; ".join(f"{k}: {v}" for k, v in r.checks.items() if v != "ok"))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
