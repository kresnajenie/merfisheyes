#!/usr/bin/env python3
"""Build bil-completed.json from the catalog.json every pipeline run writes beside its store.

    python build_bil_catalog.py ace-dip-vat ace-dry-dub ... --out bil-completed.json
    python build_bil_catalog.py --samples-file samples.csv --base scripts/new-datasets/bil-completed.json \\
        --out scripts/new-datasets/bil-completed.json

bil-completed.json (what scripts/ingest-bil-projects.ts turns into Explore
projects) used to be edited by hand. Now each dataset's part of it is written
by the pipeline (catalog_store.sbatch -> <sample>/catalog.json) and synced with
its store, so this script only gathers: per sample it reads catalog.json and
stats.json from the public bucket and makes one item whose entry is the
SpatialData store.

--base merges into an existing file, matched by bilCode: the catalog fields
and the entries are replaced, the website's own fields (isPublished,
isFeatured, isInternal, thumbnailUrl, genes, publicationLink) are kept, and
items for other codes pass through untouched. A new item starts unpublished.
--keep-entries keeps a base item's old entries (the legacy folders) after the
store instead of replacing them.

Nothing is written over --base unless --out names the same file.
Exit 0 written; 1 a sample has no stats.json in the bucket; 2 bad arguments.
"""
from __future__ import annotations

import argparse
import json
import sys
import urllib.error
import urllib.request
from pathlib import Path

BUCKET, REGION, PREFIX = "merfisheyes-bil", "us-west-2", "bil-psc-zarr"
CATALOG_FIELDS = ("title", "description", "species", "tissue", "platform", "externalLink", "metadata")
NEW_ITEM = {"isBil": True, "isPublished": False, "isFeatured": False, "isInternal": False}


def url(prefix: str, rel: str) -> str:
    return f"https://{BUCKET}.s3.{REGION}.amazonaws.com/{prefix}/{rel}"


def read_json(u: str) -> dict | None:
    try:
        with urllib.request.urlopen(u, timeout=60) as r:
            return json.loads(r.read())
    except urllib.error.HTTPError as e:
        if e.code in (403, 404):   # a public bucket answers 403 for a missing key
            return None
        raise


def item_for(sample: str, prefix: str) -> dict:
    stats = read_json(url(prefix, f"{sample}/stats.json"))
    if stats is None:
        raise LookupError(f"{sample}: no stats.json at {url(prefix, sample + '/stats.json')}")
    catalog = read_json(url(prefix, f"{sample}/catalog.json")) or {}
    if not catalog.get("title"):
        print(f"  {sample}: no catalog.json title — the item is named by its code", file=sys.stderr)
    item = {k: catalog[k] for k in CATALOG_FIELDS if catalog.get(k)}
    item.setdefault("title", sample)
    item["bilCode"] = sample
    item["entries"] = [{
        "label": "SpatialData",
        "datasetType": "spatialdata",
        "s3BaseUrl": url(prefix, f"{sample}/{stats['store']}"),
        "sortOrder": 0,
    }]
    return item


def merge(base: dict | None, new: dict, keep_entries: bool) -> dict:
    if base is None:
        return {**NEW_ITEM, **new}
    out = {**base, **new}   # catalog fields and entries from the run, the rest from the base
    if keep_entries:
        seen = {e["s3BaseUrl"] for e in new["entries"]}
        old = [e for e in base.get("entries") or [] if isinstance(e, dict) and e.get("s3BaseUrl") not in seen]
        out["entries"] = new["entries"] + [{**e, "sortOrder": i + 1} for i, e in enumerate(old)]
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("samples", nargs="*", help="BIL codes, the folders under the prefix")
    ap.add_argument("--samples-file", help="a sample list (first column), e.g. samples.csv")
    ap.add_argument("--prefix", default=PREFIX)
    ap.add_argument("--base", help="an existing bil-completed.json to merge into")
    ap.add_argument("--keep-entries", action="store_true", help="keep the base item's old entries after the store")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    samples = list(args.samples)
    if args.samples_file:
        for line in Path(args.samples_file).read_text().splitlines():
            code = line.split(",")[0].strip()
            if code and not code.startswith("#") and code != "sample_name":
                samples.append(code)
    if not samples:
        ap.error("no samples: name them or pass --samples-file")

    items = json.loads(Path(args.base).read_text()) if args.base else []
    index = {d.get("bilCode"): i for i, d in enumerate(items)}
    failed = 0
    for s in dict.fromkeys(samples):
        try:
            new = item_for(s, args.prefix.strip("/"))
        except LookupError as e:
            print(f"  {e}", file=sys.stderr)
            failed += 1
            continue
        if s in index:
            items[index[s]] = merge(items[index[s]], new, args.keep_entries)
            how = "updated"
        else:
            index[s] = len(items)
            items.append(merge(None, new, False))
            how = "added"
        print(f"  {how:8} {s:14} {new['title'][:70]}")

    Path(args.out).write_text(json.dumps(items, indent=2, ensure_ascii=False) + "\n")
    print(f"wrote {args.out}: {len(items)} item(s)" + (f", {failed} sample(s) missing" if failed else ""))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
