#!/usr/bin/env python3
"""Build the public query index: one catalog.parquet of every public dataset MERFISHEYES serves.

    python scripts/metadata/build_index.py OUT [--sddb DIR] [--explore FILE] [--no-fetch]

docs/plans/METADATA.md step 6 (M4, M9); columns: docs/specs/metadata.md §3 and
docs/specs/schemas/metadata.schema.json $defs/indexRow. Runs where lamindb is (meyes-lamin).

Rows, in SpatialData-DB's catalog format (schema 3) plus the meyes_ columns:

  ours   every store registered in Lamin merfisheyes/eyes-catalog whose own metadata passed
         bionty (meyes_provenance.validation.status == "pass") and that a public Explore row
         lists (M9: a row needs a store; M7: nothing private). Values come from the store's
         root zarr.json (meyes_metadata) and its Lamin record (uid, size, card).
  explore  every other dataset live on Explore (no store, or a store not yet checked): what the
         site knows, mapped by build_metadata's rules, meyes_metadata_from "explore", not validated.
  sddb   the rows of SpatialData-DB's release at the pinned tag (conventions.json
         metadata.sddb.release) whose store Explore lists — only what is live on MERFISHEYES —
         unchanged, with meyes_source "sddb" and their MERFISHEYES page.

Writes OUT/catalog.parquet, OUT/ontology_rollup.json (SpatialData-DB's, plus our terms),
OUT/technology.json (theirs, unchanged) and OUT/release_info.json (what went in, what was left
out and why). The result is checked with SpatialData-DB's own validator (their
_catalog_schema.py at the pinned tag) and our rows against $defs/indexRow; exit 1 if either fails.

--sddb DIR      their release files (default: downloaded to ~/.cache/merfisheyes/sddb/<tag>/)
--explore FILE  Explore's public rows as GET /api/agent/v1/catalog/search pages (default:
                fetched from www.merfisheyes.com into ~/.cache/merfisheyes/explore-rows.json, rate-limited;
                delete it to refetch)
"""
from __future__ import annotations

import argparse
import datetime as dt
import importlib.util
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
CONVENTIONS = json.loads((REPO / "docs" / "specs" / "conventions.json").read_text())
META = CONVENTIONS["metadata"]
SCHEMA = json.loads((REPO / "docs" / "specs" / "schemas" / "metadata.schema.json").read_text())
INDEX_VERSION = "1"
SITE = "https://www.merfisheyes.com"
INSTANCE = "merfisheyes/eyes-catalog"
OUR_ROOT = "s3://merfisheyes-bil/bil-psc-zarr"
OUR_HTTPS = "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/bil-psc-zarr"
SDDB_REPO = "theislab/spatialdata-db-client"
SDDB_FILES = ("catalog.parquet", "technology.json", "ontology_rollup.json", "release_info.json")
UA = {"User-Agent": "merfisheyes-index-builder"}


def get(url: str, retries: int = 4) -> bytes:
    for i in range(retries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=120) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code != 429 or i == retries - 1:
                raise
            time.sleep(int(e.headers.get("Retry-After") or 30))
    raise RuntimeError(url)


# ── inputs ────────────────────────────────────────────────────────────

def sddb_release(folder: Path | None) -> Path:
    tag = META["sddb"]["release"]
    folder = folder or Path.home() / ".cache" / "merfisheyes" / "sddb" / tag
    folder.mkdir(parents=True, exist_ok=True)
    for name in (*SDDB_FILES, "_catalog_schema.py"):
        if not (folder / name).exists():
            url = (f"https://raw.githubusercontent.com/{SDDB_REPO}/{tag}/src/sddb/_catalog_schema.py"
                   if name == "_catalog_schema.py" else f"https://github.com/{SDDB_REPO}/releases/download/{tag}/{name}")
            print(f"  fetching {name} ({tag})")
            (folder / name).write_bytes(get(url))
    return folder


def explore_rows(path: Path, fetch: bool) -> list[dict]:
    """Every public Explore row, as the agent API's search pages it (the site's own visibility rule)."""
    if path.exists():
        return json.loads(path.read_text())
    if not fetch:
        raise SystemExit(f"{path} is missing and --no-fetch was given")
    rows, offset = [], 0
    while True:
        page = json.loads(get(f"{SITE}/api/agent/v1/catalog/search?category=all&limit=100&offset={offset}"))
        rows += page["items"]
        print(f"\r  Explore rows {len(rows):,}/{page['total']:,}", end="", flush=True)
        if not page.get("has_more"):
            break
        offset = page["next_offset"]
        time.sleep(1.1)   # anonymous agent calls: 60 a minute
    print()
    path.write_text(json.dumps(rows))
    return rows


def https_of(zarr_url: str) -> str:
    """s3://bucket/key → the bucket's public HTTPS URL (conventions.json catalog.bucketRegions)."""
    m = re.match(r"s3://([^/]+)/(.+)", zarr_url or "")
    if not m:
        return ""
    bucket, key = m.groups()
    region = CONVENTIONS["catalog"]["bucketRegions"].get(bucket, "us-west-2")
    return f"https://{bucket}.s3.{region}.amazonaws.com/{key.rstrip('/')}"


def viewer_url(https: str) -> str:
    return f"{SITE}/spatialdata-viewer/from-s3?url={urllib.parse.quote(https, safe='')}"


# ── our rows ──────────────────────────────────────────────────────────

def slugify(s) -> str:
    return re.sub(r"[^a-z0-9]+", "-", str(s or "").lower()).strip("-")


def group(term_id: str | None, rollup: dict, anchors: list, fallback: str):
    """First anchor in list order that is the term or one of its ancestors (metadata.md §7.6)."""
    if not term_id:
        return None, None
    present = {term_id, *((rollup.get(term_id) or {}).get("ancestors_present") or [])}
    for aid, label in anchors:
        if aid in present:
            return label, aid
    return fallback, None


def derived(m: dict, rollup: dict) -> dict:
    """The columns the index derives from a row's own values: licence flags, the tissue and
    disease groupings (metadata.md §7.6), pathological, low_obs."""
    g = META["groupings"]
    lic = META["licenses"].get(m.get("license") or "")
    tsys, tsys_id = group(m.get("tissue_id"), rollup, g["tissueSystems"], g["fallbacks"]["tissue"])
    torg, torg_id = group(m.get("tissue_id"), rollup, g["tissueOrgans"], None)
    if torg is None:
        torg, torg_id = tsys, tsys_id
    disease = m.get("disease")
    mondo = m.get("disease_id") if str(m.get("disease_id") or "").startswith("MONDO:") else None
    dcls, dcls_id = group(mondo, rollup, g["diseaseClasses"], g["fallbacks"]["disease"])
    dgrp, dgrp_id = group(mondo, rollup, g["diseaseGroups"], g["fallbacks"]["disease"])
    sick = disease not in (None, "unknown", "normal")
    return {
        "license_url": lic["url"] if lic else None,
        "license_noncommercial": lic["noncommercial"] if lic else None,
        "license_redistributable": lic["redistributable"] if lic else None,
        "license_unknown": lic is None,
        "low_obs": None if m.get("n_obs") is None else m["n_obs"] < 100,
        "pathological": None if disease in (None, "unknown") else disease != "normal",
        "tissue_system": tsys, "tissue_system_id": tsys_id,
        "tissue_organ": torg, "tissue_organ_id": torg_id,
        "disease_class": dcls if sick else None, "disease_class_id": dcls_id if sick else None,
        "disease_group": dgrp if sick else None, "disease_group_id": dgrp_id if sick else None,
    }


def our_rows(explore_by_code: dict, rollup: dict) -> tuple[list[dict], list[dict]]:
    import lamindb as ln

    ln.connect(INSTANCE)
    ln.settings.verbosity = "error"
    arts = list(ln.Artifact.filter(storage__root=OUR_ROOT, suffix=".zarr", is_latest=True))
    rows, left_out = [], []
    for a in sorted(arts, key=lambda a: a.key):
        https = f"{OUR_HTTPS}/{a.key}"
        attrs = json.loads(get(f"{https}/zarr.json")).get("attributes", {})
        m = attrs.get("meyes_metadata")
        status = ((attrs.get("meyes_provenance") or {}).get("validation") or {}).get("status")
        code = (m or {}).get("meyes_bil_code") or a.key.split("/")[0]
        explore = explore_by_code.get(code)
        why = ("no meyes_metadata" if not m else f"validation {status}" if status != "pass"
               else "no public Explore row" if not explore else None)
        if why:
            left_out.append({"store": a.key, "why": why})
            continue
        row = {
            **{k: v for k, v in m.items() if k != "version"},
            **derived(m, rollup),
            "uid": a.uid,
            "zarr_url": f"{OUR_ROOT}/{a.key}",
            "card_url": f"https://lamin.ai/{INSTANCE}/artifact/{a.uid}",
            "size_bytes": int(a.size or 0),
            "published_at": a.created_at.isoformat(),
            "validation_status": "pass",
            "tier": None,
            "review_status": None,
            "review_reasons": [],
            "slug": f"{slugify(m.get('assay') or m.get('technology') or 'dataset')}-{slugify(m.get('tissue'))}-{a.uid.lower()}",
            "meyes_id": explore["id"],
            "meyes_source": explore["source"],
            "meyes_page_url": explore["pageUrl"],
            "meyes_viewer_url": viewer_url(https),
            "meyes_format": "spatialdata",
            "meyes_metadata_from": "store",
        }
        rows.append(row)
    return rows, left_out


# ── live datasets with no checked store ───────────────────────────────

def explore_detail(row_id: str, cache: Path) -> dict:
    f = cache / f"{row_id}.json"
    if not f.exists():
        f.write_bytes(get(f"{SITE}/api/agent/v1/catalog/{row_id}"))
        time.sleep(2.1)   # detail calls read S3: 30 a minute
    return json.loads(f.read_text())


def explore_only_rows(rows: list[dict], cache: Path, rollup: dict) -> list[dict]:
    """Every live Explore row the index has no store row for (M9 amended, M13): what the site
    knows, mapped to columns by the pipeline's own rules (build_metadata.from_catalog), marked
    meyes_metadata_from "explore" and validation_status "not_validated"."""
    sys.path.insert(0, str(HERE))
    import build_metadata as bm

    cache.mkdir(parents=True, exist_ok=True)
    out = []
    for r in rows:
        d = explore_detail(r["id"], cache) if r.get("source") != "sddb" else r
        # "source" there is the Explore category, not catalog.json's field → origin map
        cols, _ = bm.from_catalog({**d, "metadata": d.get("metadata") or {}, "source": {}})
        if d.get("disease"):
            cols["disease"] = d["disease"]
        cols = bm.complete(cols)
        entries = d.get("entries") or r.get("entries") or []
        store = next((e for e in entries if e.get("datasetType") == "spatialdata" and e.get("url")), None)
        kinds = {e.get("kind") for e in entries}
        key = (f"{META['keyPrefixes']['bil']}:{d['bilCode']}" if d.get("bilCode")
               else f"{META['keyPrefixes']['meyes']}:{r['id']}")
        m = re.match(r"https://([^.]+)\.s3\.[^/]+\.amazonaws\.com/(.+)", (store or {}).get("url") or "")
        out.append({
            **cols,
            "canonical_key": key,
            "study_id": cols.get("study_id") or key,
            "title": d.get("title"),
            "meyes_description": d.get("description"),
            "meyes_institution": d.get("institute") or cols.get("meyes_institution"),
            "meyes_tags": d.get("tags") or None,
            "n_obs": d.get("cells"),
            "n_features": d.get("genes") if isinstance(d.get("genes"), int) else d.get("geneCount") or None,
            "license": d.get("license") if d.get("license") in META["licenses"] else cols.get("license"),
            "zarr_url": f"s3://{m.group(1)}/{m.group(2).rstrip('/')}" if m else None,
            "validation_status": "not_validated",
            "meyes_id": r["id"],
            "meyes_source": r["source"],
            "meyes_page_url": r["pageUrl"],
            "meyes_kind": "single_cell" if "single_cell" in kinds else ("single_molecule" if kinds else None),
            "meyes_viewer_url": (viewer_url(store["url"]) if store else
                                 next((e.get("viewerUrl") for e in entries if e.get("viewerUrl")), None)),
            "meyes_format": "spatialdata" if store else "merfisheyes-chunked",
            "meyes_metadata_from": "explore",
        })
        out[-1].update(derived(out[-1], rollup))
    return out


def group_rereleases(df: pd.DataFrame) -> int:
    """Rows with the same meyes_matrix_hash hold the same data (M12): SpatialData-DB's
    same_section_group (the newest member's uid, else its canonical_key) and relation
    `rerelease`. Returns the number of groups."""
    if "meyes_matrix_hash" not in df:
        return 0
    n = 0
    for h, members in df[df["meyes_matrix_hash"].notna()].groupby("meyes_matrix_hash"):
        if len(members) < 2:
            continue
        newest = members.sort_values("published_at").iloc[-1]
        gid = newest["uid"] if isinstance(newest["uid"], str) else newest["canonical_key"]
        df.loc[members.index, "same_section_group"] = gid
        df.loc[members.index, "same_section_relation"] = "rerelease"
        n += 1
    return n


# ── assembly and checks ───────────────────────────────────────────────

def load_validator(folder: Path):
    spec = importlib.util.spec_from_file_location("sddb_catalog_schema", folder / "_catalog_schema.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("out", type=Path)
    ap.add_argument("--sddb", type=Path)
    ap.add_argument("--explore", type=Path)
    ap.add_argument("--no-fetch", action="store_true", help="use cached files only")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    sddb_dir = sddb_release(a.sddb)
    sddb_schema = load_validator(sddb_dir)
    theirs = pd.read_parquet(sddb_dir / "catalog.parquet")
    rollup = json.loads((sddb_dir / "ontology_rollup.json").read_text())
    explore = explore_rows(a.explore or Path.home() / ".cache" / "merfisheyes" / "explore-rows.json", not a.no_fetch)
    by_code = {r["bilCode"]: r for r in explore if r.get("bilCode")}
    by_store = {}
    for r in explore:
        for e in r.get("entries") or []:
            if e.get("url"):
                by_store.setdefault(e["url"].rstrip("/"), r)

    ours, left_out = our_rows(by_code, rollup)
    print(f"  ours: {len(ours)} rows; left out {len(left_out)}")

    # Their rows, unchanged, with our columns where Explore lists the same store
    theirs = theirs.copy()
    https = [https_of(u) if isinstance(u, str) else "" for u in theirs["zarr_url"]]
    hit = [by_store.get(u) for u in https]
    theirs["meyes_source"] = "sddb"
    theirs["meyes_kind"] = "single_cell"
    theirs["meyes_id"] = [r["id"] if r else None for r in hit]
    theirs["meyes_page_url"] = [r["pageUrl"] if r else None for r in hit]
    theirs["meyes_viewer_url"] = [viewer_url(u) if u else None for u in https]
    theirs["meyes_format"] = "spatialdata"
    theirs["meyes_metadata_from"] = "sddb"
    # The index is an account of what is live on MERFISHEYES: a SpatialData-DB row is in it only
    # when Explore lists its store (SpatialData-DB itself is the reference, not the inventory)
    not_live = int(theirs["meyes_id"].isna().sum())
    theirs = theirs[theirs["meyes_id"].notna()].reset_index(drop=True)

    # Every other live Explore row: on the site, so in the index, with what the site knows
    covered = {r["meyes_id"] for r in ours} | set(theirs["meyes_id"])
    rest = explore_only_rows([r for r in explore if r["id"] not in covered],
                             Path.home() / ".cache" / "merfisheyes" / "explore-detail", rollup)
    print(f"  live without a checked store: {len(rest)} rows (from Explore)")

    cols = list(SCHEMA["$defs"]["indexRow"]["properties"])
    ours_df = pd.DataFrame(ours + rest)
    df = pd.concat([theirs, ours_df], ignore_index=True, sort=False)
    groups = group_rereleases(df)
    extra = [c for c in df.columns if c not in cols]
    df = df.reindex(columns=cols + extra)
    for c, p in SCHEMA["$defs"]["indexRow"]["properties"].items():
        if p.get("x-sddb-dtype", "").startswith("list") or c in ("meyes_sections", "meyes_de_columns",
                                                                   "meyes_co_investigators", "meyes_tags",
                                                                   "meyes_cell_type_levels", "modalities"):
            continue
        t = p.get("x-sddb-dtype")
        if t in ("Int64", "Float64", "boolean", "string"):
            df[c] = df[c].astype(t)
    # list columns: a list, or None — never NaN (an all-empty column would otherwise be float)
    for c, p in SCHEMA["$defs"]["indexRow"]["properties"].items():
        if p.get("type") == "array" or "array" in (p.get("type") or []) or str(p.get("x-sddb-dtype", "")).startswith("list"):
            df[c] = pd.Series([list(v) if isinstance(v, (list, tuple)) or hasattr(v, "tolist") and not isinstance(v, str)
                               else None for v in df[c]], dtype=object)
    for c in ("meyes_dims", "meyes_n_sections", "meyes_n_molecule_genes"):
        df[c] = df[c].astype("Int64")
    for c in ("timepoint_days",):
        df[c] = df[c].astype("Float64")

    # Checks: their validator on the whole index, our schema on our rows
    problems = []
    try:
        sddb_schema.validate(df, kind="catalog")
    except Exception as e:  # noqa: BLE001
        problems.append(f"SpatialData-DB validator: {e}")
    import jsonschema

    row_schema = {"$defs": SCHEMA["$defs"], **SCHEMA["$defs"]["indexRow"]}
    v = jsonschema.Draft202012Validator(row_schema)
    loose = jsonschema.Draft202012Validator({**row_schema, "required": ["meyes_id", "meyes_source", "meyes_page_url",
                                                                         "title", "meyes_format", "meyes_metadata_from"]})
    for r in ours + rest:
        clean = {k: val for k, val in r.items() if val is not None}
        for e in (v if r.get("meyes_metadata_from") == "store" else loose).iter_errors(clean):
            problems.append(f"{r.get('canonical_key')}: {'/'.join(map(str, e.absolute_path))}: {e.message}")

    # Write
    released = dt.date.today().isoformat()
    info = {
        "catalog_schema_version": 3,
        "meyes_index_version": INDEX_VERSION,
        "release": f"meyes-index-{released}",
        "built_at": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "sddb_release": META["sddb"]["release"],
        "sddb_release_info": json.loads((sddb_dir / "release_info.json").read_text()),
        "rows": int(len(df)),
        "rows_by_source": {k: int(v) for k, v in df["meyes_source"].value_counts().items()},
        "sddb_rows_on_explore": int(len(theirs)),
        "sddb_rows_not_on_explore": not_live,
        "explore_public_rows": len(explore),
        "rows_by_metadata_from": {k: int(v) for k, v in df["meyes_metadata_from"].value_counts().items()},
        "rows_by_format": {k: int(v) for k, v in df["meyes_format"].value_counts().items()},
        "rows_from_explore_only": [{"id": r["meyes_id"], "source": r["meyes_source"], "title": r["title"],
                                    "format": r["meyes_format"]} for r in rest],
        "rerelease_groups": groups,
        "every_live_row_indexed": len(explore) == int(df["meyes_id"].notna().sum()),
        "left_out_stores": left_out,
        "checks": problems or "pass",
    }
    import pyarrow as pa
    import pyarrow.parquet as pq

    table = pa.Table.from_pandas(df, preserve_index=False)
    table = table.replace_schema_metadata({**(table.schema.metadata or {}),
                                           b"meyes_release_info": json.dumps(info).encode()})
    pq.write_table(table, a.out / "catalog.parquet", compression="snappy")   # as SpatialData-DB writes theirs; hyparquet reads it without plugins
    merged = {**rollup}
    for r in ours:
        for col in ("tissue_id", "disease_id"):
            t = r.get(col)
            if t and t not in merged:
                merged[t] = {"ancestors_present": [], "name": r.get(col.replace("_id", ""))}
    (a.out / "ontology_rollup.json").write_text(json.dumps(merged, indent=1, sort_keys=True))
    (a.out / "technology.json").write_bytes((sddb_dir / "technology.json").read_bytes())
    (a.out / "release_info.json").write_text(json.dumps(info, indent=1, ensure_ascii=False))
    print(f"  {a.out}/catalog.parquet: {len(df):,} rows "
          f"({', '.join(f'{k} {v}' for k, v in info['rows_by_source'].items())}); "
          f"{info['sddb_rows_not_on_explore']:,} SpatialData-DB rows left out (not on Explore); "
          f"{len(rest)} from Explore only; every live row indexed: {info['every_live_row_indexed']}; "
          f"re-release groups: {groups}")
    for p in problems[:20]:
        print(f"  FAIL {p}")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
