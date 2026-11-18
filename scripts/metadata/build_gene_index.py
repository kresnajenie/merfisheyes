#!/usr/bin/env python3
"""Build the gene index beside the query index: which live datasets measure which genes.

    python scripts/metadata/build_gene_index.py INDEX_DIR --stores ROOT [--stores ROOT ...]

docs/plans/METADATA.md step 7 (M5). INDEX_DIR holds the catalog.parquet that
build_index.py wrote; this adds, in SpatialData-DB's formats (their _catalog_schema.py):

  gene_index.parquet      one row per (our dataset, gene): symbol, feature_id, uid,
                          total_counts, fraction_obs_detected (GENE_INDEX_COLUMNS), plus
                          meyes_id. For a store (meyes_metadata_from "store") the two numbers
                          come from its table's X; for a dataset the index knows only from
                          Explore, the genes come from the site's gene list and the numbers
                          are null (measured, not quantified). SpatialData-DB's own rows are
                          in their gene_index.parquet, which this does not copy.
  gene_search.parquet     one row per gene entity over EVERY index row — ours and
                          SpatialData-DB's live ones: entity_id, feature_id, symbols,
                          n_datasets, bitset (FEATURE_SEARCH-style; bit i, little-endian,
                          is catalog.parquet's row i). Their release's bitsets, re-ordered to
                          our rows (their off-site datasets dropped), joined with ours.
  gene_search_meta.json   meyes_ids and dataset_uids in row order, order_hash (SHA-256 of the
                          meyes_ids joined by "," — SpatialData-DB hashes their uids the same
                          way; not every row of ours has a uid), n_entities, n_datasets.

A gene of ours joins their entity by symbol when one of their entities for the dataset's
organism (Ensembl prefix) lists it; else its entity is the symbol itself, as theirs are for
symbols they could not resolve.

--stores ROOT  where the local copies of our stores are: ROOT/<code>/**/<name>.zarr
               (e.g. ~/data/bil-data/processed); a store not found there is reported.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
import urllib.request
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
META = json.loads((REPO / "docs" / "specs" / "conventions.json").read_text())["metadata"]
SITE = "https://www.merfisheyes.com"
UA = {"User-Agent": "merfisheyes-index-builder"}
# Ensembl gene-id prefix of each organism we serve (to pick the right entity for a symbol)
PREFIX = {"human": "ENSG0", "mouse": "ENSMUSG", "macaque": "ENSMMUG", "white-tufted-ear marmoset": "ENSCJAG",
          "rat": "ENSRNOG", "zebrafish": "ENSDARG"}
CHUNK = 20_000


def get_json(url: str):
    for i in range(4):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=120) as r:
                return json.loads(r.read())
        except urllib.error.HTTPError as e:
            if e.code != 429 or i == 3:
                raise
            time.sleep(int(e.headers.get("Retry-After") or 30))


def explore_genes(meyes_id: str, cache: Path) -> list[str]:
    f = cache / f"{meyes_id}.genes.json"
    if not f.exists():
        genes, offset = [], 0
        while True:
            p = get_json(f"{SITE}/api/agent/v1/catalog/{meyes_id}/genes?limit=1000&offset={offset}")
            genes += p["items"]
            if not p.get("has_more"):
                break
            offset = p["next_offset"]
        f.write_text(json.dumps(genes))
        time.sleep(1.1)
    return json.loads(f.read_text())


def find_store(zarr_url: str, roots: list[Path]) -> Path | None:
    parts = zarr_url.rstrip("/").split("/")
    code, name = parts[-2], parts[-1]
    for root in roots:
        hits = sorted((root / code).glob(f"**/{name}")) if (root / code).exists() else []
        if hits:
            return hits[0]
    return None


def store_genes(store: Path) -> pd.DataFrame:
    """symbol, total_counts, fraction_obs_detected for the store's default table."""
    import anndata as ad
    import scipy.sparse as sp

    sys.path.insert(0, str(REPO / "scripts"))
    import prepare_folder

    a = ad.read_zarr(str(store / "tables" / prepare_folder.default_table_name(store)))
    x = sp.csc_matrix(a.X)
    return pd.DataFrame({
        "symbol": list(a.var_names),
        "total_counts": np.asarray(x.sum(axis=0)).ravel().astype(float),
        "fraction_obs_detected": np.diff(x.indptr) / max(a.n_obs, 1),
    })


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("index_dir", type=Path)
    ap.add_argument("--stores", type=Path, action="append", default=[])
    ap.add_argument("--sddb", type=Path, help="their release files (default: the cache build_index.py uses)")
    a = ap.parse_args()

    tag = META["sddb"]["release"]
    sddb = a.sddb or Path.home() / ".cache" / "merfisheyes" / "sddb" / tag
    for name in ("gene_search.parquet", "gene_search_meta.json"):
        if not (sddb / name).exists():
            print(f"  fetching {name} ({tag})")
            req = urllib.request.Request(f"https://github.com/theislab/spatialdata-db-client/releases/download/{tag}/{name}", headers=UA)
            (sddb / name).write_bytes(urllib.request.urlopen(req, timeout=600).read())

    cat = pd.read_parquet(a.index_dir / "catalog.parquet",
                          columns=["uid", "meyes_id", "organism", "zarr_url", "meyes_metadata_from", "canonical_key"])
    n = len(cat)
    row_of_uid = {u: i for i, u in enumerate(cat["uid"]) if isinstance(u, str)}

    # 1. Their entities: symbol → entity ids, by organism prefix
    theirs = pq.read_table(sddb / "gene_search.parquet")
    tmeta = json.loads((sddb / "gene_search_meta.json").read_text())
    ent_ids = theirs.column("entity_id").to_pylist()
    by_symbol: dict[str, list[str]] = defaultdict(list)
    for eid, syms in zip(ent_ids, theirs.column("symbols").to_pylist()):
        for s in syms or []:
            by_symbol[s].append(eid)

    def entity_for(symbol: str, organism: str | None) -> str:
        cands = by_symbol.get(symbol.strip().lower(), [])
        pre = PREFIX.get(organism or "")
        for e in cands:
            if pre and e.startswith(pre):
                return e
        return symbol.strip()

    # 2. Our genes, per index row
    cache = Path.home() / ".cache" / "merfisheyes" / "explore-detail"
    cache.mkdir(parents=True, exist_ok=True)
    gi_parts, ours_bits = [], defaultdict(set)
    our_symbols: dict[str, set] = defaultdict(set)
    missing = []
    ours = cat[cat["meyes_metadata_from"] != "sddb"]
    for i, r in ours.iterrows():
        if r["meyes_metadata_from"] == "store":
            store = find_store(r["zarr_url"], a.stores)
            if store is None:
                missing.append(r["canonical_key"])
                continue
            g = store_genes(store)
        else:
            g = pd.DataFrame({"symbol": explore_genes(r["meyes_id"], cache)})
            g["total_counts"] = np.nan
            g["fraction_obs_detected"] = np.nan
        g["feature_id"] = [entity_for(s, r["organism"]) for s in g["symbol"]]
        g["uid"] = r["uid"] if isinstance(r["uid"], str) else None
        g["meyes_id"] = r["meyes_id"]
        gi_parts.append(g)
        for s, e in zip(g["symbol"], g["feature_id"]):
            ours_bits[e].add(i)
            our_symbols[e].add(s.strip().lower())
        print(f"\r  our datasets: {len(gi_parts)}/{len(ours)}", end="", flush=True)
    print()
    gi = pd.concat(gi_parts, ignore_index=True)[["symbol", "feature_id", "uid", "total_counts",
                                                  "fraction_obs_detected", "meyes_id"]]
    gi = gi.astype({"symbol": "string", "feature_id": "string", "uid": "string",
                    "total_counts": "Float64", "fraction_obs_detected": "Float64", "meyes_id": "string"})
    pq.write_table(pa.Table.from_pandas(gi, preserve_index=False), a.index_dir / "gene_index.parquet", compression="zstd")

    # 3. gene_search over every index row: theirs re-ordered to our rows, ours OR-ed in
    src_cols = np.array([row_of_uid.get(u, -1) for u in tmeta["dataset_uids"]])
    keep = src_cols >= 0
    nbytes_in = (len(tmeta["dataset_uids"]) + 7) // 8
    nbytes = (n + 7) // 8
    out_rows = []   # (entity_id, feature_id, symbols, n, bitset)
    seen = set()
    bits_col = theirs.column("bitset")
    feats = theirs.column("feature_id").to_pylist()
    syms_all = theirs.column("symbols").to_pylist()
    for start in range(0, len(ent_ids), CHUNK):
        stop = min(start + CHUNK, len(ent_ids))
        raw = np.frombuffer(b"".join(bits_col.slice(start, stop - start).to_pylist()), dtype=np.uint8)
        old = np.unpackbits(raw.reshape(stop - start, nbytes_in), axis=1, bitorder="little")[:, :len(src_cols)]
        new = np.zeros((stop - start, n), dtype=np.uint8)
        new[:, src_cols[keep]] = old[:, keep]
        for k in range(stop - start):
            eid = ent_ids[start + k]
            if eid in ours_bits:
                new[k, list(ours_bits[eid])] = 1
                seen.add(eid)
        counts = new.sum(axis=1)
        packed = np.packbits(new, axis=1, bitorder="little")
        for k in np.flatnonzero(counts):
            eid = ent_ids[start + k]
            syms = sorted(set(syms_all[start + k] or []) | our_symbols.get(eid, set()))
            out_rows.append((eid, feats[start + k], syms, int(counts[k]), packed[k].tobytes()))
        print(f"\r  their entities: {stop:,}/{len(ent_ids):,}", end="", flush=True)
    print()
    for eid in sorted(set(ours_bits) - seen):
        row = np.zeros(n, dtype=np.uint8)
        row[list(ours_bits[eid])] = 1
        out_rows.append((eid, eid, sorted(our_symbols[eid]), int(row.sum()),
                         np.packbits(row, bitorder="little").tobytes()))
    table = pa.table({
        "entity_id": pa.array([r[0] for r in out_rows], pa.string()),
        "feature_id": pa.array([r[1] for r in out_rows], pa.string()),
        "symbols": pa.array([r[2] for r in out_rows], pa.list_(pa.string())),
        "n_datasets": pa.array([r[3] for r in out_rows], pa.int32()),
        "bitset": pa.array([r[4] for r in out_rows], pa.binary()),
    })
    pq.write_table(table, a.index_dir / "gene_search.parquet", compression="zstd")
    ids = [str(x) for x in cat["meyes_id"]]
    meta = {"meyes_ids": ids, "dataset_uids": [u if isinstance(u, str) else None for u in cat["uid"]],
            "order_hash": hashlib.sha256(",".join(ids).encode()).hexdigest(),
            "n_entities": len(out_rows), "n_datasets": n, "nbytes": nbytes,
            "sddb_release": tag, "bit_order": "little"}
    (a.index_dir / "gene_search_meta.json").write_text(json.dumps(meta))

    # 4. Check with SpatialData-DB's own validator
    import importlib.util

    spec = importlib.util.spec_from_file_location("s", sddb / "_catalog_schema.py")
    s = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(s)
    problems = []
    try:
        s.validate(gi.drop(columns=["meyes_id"]), kind="gene_index")
    except Exception as e:  # noqa: BLE001
        problems.append(f"gene_index: {e}")
    new_ents = len(set(ours_bits) - seen)
    print(f"  gene_index.parquet: {len(gi):,} rows over {gi['meyes_id'].nunique()} of our datasets "
          f"({gi['total_counts'].notna().sum():,} with counts)")
    print(f"  gene_search.parquet: {len(out_rows):,} entities over {n:,} datasets "
          f"({len(seen):,} of ours joined their entities, {new_ents:,} new); missing stores: {missing or 'none'}")
    for p in problems:
        print(f"  FAIL {p}")
    return 1 if problems or missing else 0


if __name__ == "__main__":
    sys.exit(main())
