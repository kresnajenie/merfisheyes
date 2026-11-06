#!/usr/bin/env python3
"""Write catalog.json for one BIL dataset, beside its store, on every pipeline run.

    python catalog_entry.py <bil_code> <store_dir> [--catalog bil-completed.json]

catalog.json is one item of bil-completed.json — title, description, species,
tissue, platform, externalLink, metadata{investigator, institution, ...} —
for this dataset alone, so it travels with the store (s3_sync_store.sbatch
uploads it) and the Lamin record (register_lamin.py) and the website catalog
(build_bil_catalog.py) are both made from it, never typed in by hand.

Where each field comes from, first that answers:

  1. the BIL API (api.brainimagelibrary.org/retrieve?bildid=<code>)
     — metadata through prisma/scrape-bil-metadata.py, the one parser of it
  2. --catalog: this code's item in an existing bil-completed.json
  3. bil-merfish-datasets.csv beside this script (scrape_bil_merfish.py)

`source` in the output says which one each part came from. When none knows
the code the file is still written, with only bilCode, and the exit is 0 with
a warning: a missing title must not stop a two-day run from publishing.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import importlib.util
import json
import sys
import urllib.error
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
RETRIEVE = "https://api.brainimagelibrary.org/retrieve?bildid={code}"
VIEW = "https://api.brainimagelibrary.org/web/view?bildid={code}"
FIELDS = ("title", "description", "species", "tissue", "platform", "externalLink")


def _metadata_parser():
    """prisma/scrape-bil-metadata.py's fetch_bil_metadata — not a second copy of it."""
    path = REPO / "prisma" / "scrape-bil-metadata.py"
    spec = importlib.util.spec_from_file_location("scrape_bil_metadata", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.fetch_bil_metadata


def from_api(code: str) -> dict:
    try:
        req = urllib.request.Request(RETRIEVE.format(code=code), headers={"Accept": "application/json"})
        with urllib.request.urlopen(req, timeout=30) as r:
            ret = (json.loads(r.read()).get("retjson") or [{}])[0]
    except (urllib.error.URLError, TimeoutError, json.JSONDecodeError, ValueError) as e:
        print(f"  BIL API: {e}", file=sys.stderr)
        return {}
    ds = (ret.get("Dataset") or [{}])[0]
    sp = (ret.get("Specimen") or [{}])[0]
    if not ds.get("title"):
        return {}
    item = {
        "title": ds.get("title"),
        "description": ds.get("abstract") or None,
        "species": (sp.get("species") or "").title() or None,
        "tissue": sp.get("organname") or None,
        "platform": ds.get("technique") or None,
        "externalLink": VIEW.format(code=code),
    }
    try:
        item["metadata"] = _metadata_parser()(code)
    except Exception as e:   # the fields above already came back; metadata is best-effort
        print(f"  BIL API metadata: {e}", file=sys.stderr)
    return item


def from_catalog(code: str, path: str | None) -> dict:
    if not path or not Path(path).exists():
        return {}
    for d in json.load(open(path)):
        if d.get("bilCode") == code:
            return {k: d[k] for k in (*FIELDS, "metadata") if d.get(k)}
    return {}


def from_csv(code: str) -> dict:
    path = HERE / "bil-merfish-datasets.csv"
    if not path.exists():
        return {}
    with path.open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if row.get("bilCode") == code:
                return {k: row[k] for k in FIELDS if row.get(k)}
    return {}


def build(code: str, catalog: str | None, offline: bool) -> dict:
    sources = [("bil-csv", from_csv(code)),
               ("bil-completed.json", from_catalog(code, catalog))]
    if not offline:
        sources.append(("bil-api", from_api(code)))
    item: dict = {"bilCode": code}
    source: dict = {}
    for name, got in sources:   # later sources win: the API over a file over the CSV
        for k, v in got.items():
            if k == "metadata":
                if v:
                    item["metadata"] = {**item.get("metadata", {}), **v}
                    source["metadata"] = name
            elif v:
                item[k] = v
                source[k] = name
    item["source"] = source
    item["written_at"] = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")
    return item


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("bil_code")
    ap.add_argument("store_dir", help="the folder holding <name>.zarr; catalog.json is written there")
    ap.add_argument("--catalog", help="an existing bil-completed.json to fall back on")
    ap.add_argument("--offline", action="store_true", help="do not ask the BIL API")
    args = ap.parse_args()

    out = Path(args.store_dir)
    if not out.is_dir():
        print(f"error: {out} does not exist — run this after prepare", file=sys.stderr)
        return 3
    item = build(args.bil_code, args.catalog, args.offline)
    (out / "catalog.json").write_text(json.dumps(item, indent=2, ensure_ascii=False) + "\n")

    print(f"{out / 'catalog.json'}")
    for k in (*FIELDS, "metadata"):
        v = item.get(k)
        shown = (", ".join(f"{a}={b}" for a, b in v.items()) if isinstance(v, dict) else str(v or "—"))
        print(f"  {k:13} {shown[:90]}  [{item['source'].get(k, 'missing')}]")
    if "title" not in item:
        print(f"WARNING: no source knows {args.bil_code}; catalog.json holds only the code", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
