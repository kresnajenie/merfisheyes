#!/usr/bin/env python3
"""Check the JSON files written beside a store — stats.json, qc.json,
catalog.json, mapping.json — against their schemas (docs/specs/schemas/)
and against each other.

  python scripts/validate_sidecars.py <folder|url> [<folder|url> …] [--json]

A folder is the one holding <name>.zarr; a URL is its public HTTPS base
(e.g. https://merfisheyes-bil.s3.us-west-2.amazonaws.com/bil-psc-zarr/ace-dip-vat/).
A missing file is not a failure: legacy and non-BIL folders have none of
them. Prints one line per check — OK, WARN or FAIL — and exits 1 when any
check fails. The agent API reads these files (docs/specs/MCP-PROTOCOL.md),
so this is the gate for changing what writes them.
"""
from __future__ import annotations

import argparse
import json
import sys
import urllib.error
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from validate_store import Report, schema_errors  # noqa: E402

FILES = {
    "stats.json": "stats.schema.json",
    "qc.json": "qc.schema.json",
    "catalog.json": "catalog.schema.json",
    "mapping.json": "mapping.schema.json",
}


def read(base: str, name: str):
    """The parsed file, None when absent; raises on unreadable JSON."""
    if base.startswith(("http://", "https://")):
        try:
            with urllib.request.urlopen(base.rstrip("/") + "/" + name, timeout=30) as resp:
                return json.loads(resp.read())
        except urllib.error.HTTPError as e:
            if e.code in (403, 404):  # a public bucket answers 403 for a missing key
                return None
            raise
    path = Path(base) / name
    return json.loads(path.read_text()) if path.is_file() else None


def check(base: str) -> Report:
    r = Report()
    docs = {}
    for name, schema in FILES.items():
        try:
            doc = read(base, name)
        except Exception as e:  # noqa: BLE001
            r.fail(name, f"unreadable: {e}")
            continue
        if doc is None:
            r.warn(name, "absent")
            continue
        errors = schema_errors(doc, schema)
        for e in errors:
            r.fail(name, e)
        if not errors:
            r.ok(name, "matches its schema")
            docs[name] = doc

    stats, qc, catalog, mapping = (docs.get(n) for n in FILES)
    codes = {n: c for n, c in (("qc.json", (qc or {}).get("dataset")),
                                ("catalog.json", (catalog or {}).get("bilCode"))) if c}
    if len(set(codes.values())) > 1:
        r.fail("consistency", f"BIL codes disagree: {codes}")
    elif codes:
        r.ok("consistency", f"BIL code {next(iter(codes.values()))}")
    if stats and mapping and mapping.get("cells") != stats.get("store"):
        r.fail("consistency", f"mapping.json cells {mapping.get('cells')!r} is not stats.json store {stats.get('store')!r}")
    if stats and qc and "cells" in stats and stats["cells"] != qc["n_cells"]:
        r.warn("consistency", f"stats.json cells {stats['cells']:,} != qc.json n_cells {qc['n_cells']:,}")
    if qc and qc.get("mmc_levels"):
        missing = [lv for lv in qc["mmc_levels"].split(",") if f"n_{lv}" not in qc]
        if missing:
            r.warn("qc.json", f"mmc_levels names {missing} with no n_<level> column")
    return r


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("bases", nargs="+", metavar="folder|url")
    ap.add_argument("--json", action="store_true", help="one JSON object per base instead of lines")
    args = ap.parse_args(argv)

    failed = False
    for base in args.bases:
        r = check(base)
        failed |= any(level == "FAIL" for level, _, _ in r.lines)
        if args.json:
            print(json.dumps({"base": base, "checks": [dict(zip(("level", "where", "message"), l)) for l in r.lines]}))
        else:
            print(f"== {base}")
            for level, where, msg in r.lines:
                print(f"{level:4} {where}: {msg}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
