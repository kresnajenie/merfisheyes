#!/usr/bin/env python3
"""
collect_dataset_stats.py

Collect per-dataset preprocessing statistics into one growing CSV, one row per
dataset. Scrapes whichever outputs exist under ${MEYES_BASE}/<dataset>/:

  meyes_output/manifest.json     → n_cells, n_genes_detected  (primary)
  sm_output/manifest.json.gz     → n_genes_imaged (total molecules)
  combined_output/               → legacy CSV path (MERFISH combine step)
  --h5ad <file>                  → single-sample h5ad fallback
  mmc_output/mapping_output.csv  → MapMyCells hierarchical
  mmc_output_corr/               → MapMyCells correlation → mean_corr

Behavior: any dataset scraped this run REPLACES its existing rows; untouched
datasets are left as-is. Writes are fcntl-locked so concurrent jobs don't
clobber the file.

Usage:
  python collect_dataset_stats.py --datasets ace-pie-bit ace-pie-bin \\
      --species mouse \\
      --out /bil/data/meyes/_dataset_stats/dataset_stats.csv
"""
from __future__ import annotations

import argparse
import contextlib
import datetime as dt
import fcntl
import gzip
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd


DATASET_BASE_ORDER = [
    "dataset", "species", "n_samples", "n_cells", "n_genes", "n_genes_detected",
    "n_genes_imaged",
    "total_counts", "mean_counts_per_cell", "median_counts_per_cell",
    "mask_percentile", "n_artifact", "pct_artifact",
    "mmc_methods", "mmc_levels", "mean_corr",
]
CHUNK = 100_000


# ─────────────────────────────────────────────────────────────
# meyes_output manifest scraper
# ─────────────────────────────────────────────────────────────
def _scrape_meyes_manifest(meyes_output: Path) -> dict | None:
    """Read meyes_output/manifest.json. Returns None if not present."""
    p = meyes_output / "manifest.json"
    if not p.exists():
        return None
    with open(p) as f:
        m = json.load(f)
    stats = m.get("statistics", {})
    return {
        "n_cells":          stats.get("total_cells"),
        "n_genes_detected": stats.get("total_genes"),
        "dataset_type":     m.get("type"),
    }


# ─────────────────────────────────────────────────────────────
# sm_output manifest scraper
# ─────────────────────────────────────────────────────────────
def _scrape_sm_manifest(sm_output: Path) -> dict | None:
    """Read sm_output/manifest.json.gz. Returns None if not present."""
    p = sm_output / "manifest.json.gz"
    if not p.exists():
        return None
    with gzip.open(p, "rt", encoding="utf-8") as f:
        m = json.load(f)
    stats = m.get("statistics", {})
    return {
        "n_genes_imaged": stats.get("total_molecules"),
    }


# ─────────────────────────────────────────────────────────────
# Legacy combine-derived metrics (CSV or h5ad)
# ─────────────────────────────────────────────────────────────
def _samples_and_cells(meta_path: Path):
    header = pd.read_csv(meta_path, nrows=0).columns.tolist()
    if "_sample_id" in header:
        sid = pd.read_csv(meta_path, usecols=["_sample_id"])
        return int(sid["_sample_id"].nunique()), int(len(sid))
    n_cells = sum(len(c) for c in pd.read_csv(meta_path, usecols=[0], chunksize=CHUNK))
    return pd.NA, int(n_cells)


def _gene_stats_csv(cbg_path: Path):
    gene_cols = None
    detected = None
    per_cell_sums = []
    for chunk in pd.read_csv(cbg_path, chunksize=CHUNK):
        cols = [c for c in chunk.columns if c != "cell"]
        if gene_cols is None:
            gene_cols = cols
            detected = np.zeros(len(gene_cols), dtype=bool)
        mat = chunk[gene_cols].to_numpy()
        per_cell_sums.append(mat.sum(axis=1))
        detected |= (mat > 0).any(axis=0)
    if gene_cols is None:
        raise ValueError(f"{cbg_path} has no data rows")
    per_cell = np.concatenate(per_cell_sums)
    return _pack_gene_stats(len(gene_cols), int(detected.sum()), per_cell)


def _gene_stats_h5ad(h5ad_path: Path):
    import anndata
    import scipy.sparse as sp

    adata = anndata.read_h5ad(h5ad_path)
    X = adata.X
    sample = X[:10].toarray() if sp.issparse(X) else np.asarray(X[:10])
    is_raw = np.issubdtype(X.dtype, np.integer) or np.all(sample == np.floor(sample))
    if not is_raw and "X_raw" in adata.obsm:
        X = adata.obsm["X_raw"]
    if sp.issparse(X):
        per_cell = np.asarray(X.sum(axis=1)).ravel()
        n_detected = int(np.asarray((X != 0).sum(axis=0)).ravel().astype(bool).sum())
    else:
        X = np.asarray(X)
        per_cell = X.sum(axis=1)
        n_detected = int((X > 0).any(axis=0).sum())
    return (1, X.shape[0]) + _pack_gene_stats(X.shape[1], n_detected, per_cell)


def _pack_gene_stats(n_genes, n_detected, per_cell):
    return (
        n_genes,
        n_detected,
        float(per_cell.sum()),
        round(float(per_cell.mean()), 3),
        round(float(np.median(per_cell)), 3),
    )


def _artifact_stats(combined: Path, percentile: int, n_cells: int):
    mask_path = combined / f"artifact_mask_p{percentile}.csv"
    if not mask_path.exists():
        return pd.NA, pd.NA
    cols = pd.read_csv(mask_path, nrows=0).columns.tolist()
    flag_col = "is_artifact" if "is_artifact" in cols else cols[-1]
    flags = pd.read_csv(mask_path, usecols=[flag_col])[flag_col]
    is_art = flags.astype(str).str.strip().str.lower().isin(("true", "1"))
    n_art = int(is_art.sum())
    pct = round(n_art / n_cells * 100, 3) if n_cells else pd.NA
    return n_art, pct


# ─────────────────────────────────────────────────────────────
# MapMyCells output parsing
# ─────────────────────────────────────────────────────────────
def _read_mmc(csv_path: Path):
    if not csv_path.exists():
        return None
    method = "hierarchical"
    readable = None
    with open(csv_path) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            low = line.lower()
            if "mapping method:" in low:
                method = line.split("mapping method:", 1)[1].strip().split()[0]
                method = method.rstrip(";").strip("'\"")
            if "readable taxonomy hierarchy" in low:
                with contextlib.suppress(Exception):
                    readable = json.loads(line.split("=", 1)[1].strip())
    df = pd.read_csv(csv_path, comment="#")
    if readable:
        levels = [lvl for lvl in readable if f"{lvl}_name" in df.columns]
    else:
        levels = [c[:-len("_name")] for c in df.columns if c.endswith("_name")]
    if not levels:
        return None
    return {"method": method, "df": df, "levels": levels}


def _mmc_summary(hier, corr):
    out = {}
    methods = []
    if hier:
        methods.append("hierarchical")
    if corr:
        methods.append("correlation")
    out["mmc_methods"] = ",".join(methods) if methods else pd.NA
    canonical = hier or corr
    if canonical:
        out["mmc_levels"] = ",".join(canonical["levels"])
        df = canonical["df"]
        for level in canonical["levels"]:
            out[f"n_{level}"] = int(df[f"{level}_name"].nunique())
    if hier:
        df = hier["df"]
        for level in hier["levels"]:
            col = f"{level}_bootstrapping_probability"
            if col in df.columns:
                vals = pd.to_numeric(df[col], errors="coerce")
                out[f"mean_boot_{level}"] = round(float(vals.mean()), 4)
    if corr and "correlation_coefficient" in corr["df"].columns:
        vals = pd.to_numeric(corr["df"]["correlation_coefficient"], errors="coerce")
        out["mean_corr"] = round(float(vals.mean()), 4)
    return out


# ─────────────────────────────────────────────────────────────
# Per-dataset scrape
# ─────────────────────────────────────────────────────────────
def scrape_dataset(name, meyes_base, percentile, species, h5ad_path, now):
    base = Path(meyes_base) / name
    meyes_output = base / "meyes_output"
    sm_output    = base / "sm_output"
    combined     = base / "combined_output"

    # ── cell/gene stats: meyes manifest → combined CSV → h5ad ──
    meyes_info = _scrape_meyes_manifest(meyes_output)
    n_samples = n_genes = n_genes_detected = pd.NA
    total_counts = mean_pc = median_pc = pd.NA
    n_art = pct_art = pd.NA

    if meyes_info is not None:
        n_cells         = meyes_info["n_cells"]
        n_genes_detected = meyes_info["n_genes_detected"]
        n_genes          = n_genes_detected
        # artifact stats still available if combined_output is present
        if combined.is_dir():
            n_art, pct_art = _artifact_stats(combined, percentile, n_cells)

    elif combined.is_dir() and (combined / "cell_metadata.csv").exists() \
            and (combined / "cell_by_gene.csv").exists():
        n_samples, n_cells  = _samples_and_cells(combined / "cell_metadata.csv")
        n_genes, n_genes_detected, total_counts, mean_pc, median_pc = \
            _gene_stats_csv(combined / "cell_by_gene.csv")
        n_art, pct_art = _artifact_stats(combined, percentile, n_cells)

    elif h5ad_path is not None and Path(h5ad_path).exists():
        n_samples, n_cells, n_genes, n_genes_detected, total_counts, mean_pc, median_pc = \
            _gene_stats_h5ad(Path(h5ad_path))

    else:
        # SM-only dataset (no cell segmentation output yet)
        sm_info = _scrape_sm_manifest(sm_output)
        if sm_info is None:
            raise FileNotFoundError(
                f"No usable output found for {name} "
                f"(checked meyes_output, combined_output, h5ad, sm_output)")
        print(f"[sm   ] {name}: SM-only dataset, no cell stats available")
        n_cells = pd.NA

    # ── molecule count from sm manifest ────────────────────────
    sm_info = _scrape_sm_manifest(sm_output)
    n_genes_imaged = sm_info["n_genes_imaged"] if sm_info else pd.NA

    # ── MMC ────────────────────────────────────────────────────
    hier = _read_mmc(base / "mmc_output"      / "mapping_output.csv")
    corr = _read_mmc(base / "mmc_output_corr" / "mapping_output.csv")

    row = {
        "dataset":               name,
        "species":               species or pd.NA,
        "n_samples":             n_samples,
        "n_cells":               n_cells,
        "n_genes":               n_genes,
        "n_genes_detected":      n_genes_detected,
        "n_genes_imaged":        n_genes_imaged,
        "total_counts":          total_counts,
        "mean_counts_per_cell":  mean_pc,
        "median_counts_per_cell": median_pc,
        "mask_percentile":       percentile,
        "n_artifact":            n_art,
        "pct_artifact":          pct_art,
        "scraped_at":            now,
    }
    row.update(_mmc_summary(hier, corr))
    return row


# ─────────────────────────────────────────────────────────────
# Locked accumulator merge
# ─────────────────────────────────────────────────────────────
@contextlib.contextmanager
def _lock(out_path: Path):
    lock_path = out_path.with_name(out_path.name + ".lock")
    with open(lock_path, "w") as fh:
        fcntl.flock(fh, fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(fh, fcntl.LOCK_UN)


def _order_cols(df, preferred):
    front = [c for c in preferred if c in df.columns]
    rest  = sorted(c for c in df.columns if c not in front)
    return df[front + rest]


def _merge_write(out_path: Path, new_df: pd.DataFrame, preferred, sort_by):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    scraped = set(new_df["dataset"])
    with _lock(out_path):
        if out_path.exists():
            existing = pd.read_csv(out_path)
            existing = existing[~existing["dataset"].isin(scraped)]
        else:
            existing = pd.DataFrame()
        combined = pd.concat([existing, new_df], ignore_index=True)
        combined = _order_cols(combined, preferred)
        combined = combined.sort_values(sort_by).reset_index(drop=True)
        combined.to_csv(out_path, index=False)
    return combined


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--datasets",        nargs="+", required=True)
    p.add_argument("--out",             required=True, help="dataset_stats.csv accumulator")
    p.add_argument("--species",         default="")
    p.add_argument("--mask-percentile", type=int, default=25)
    p.add_argument("--h5ad",            default=None,
                   help="h5ad file for combine-derived metrics when no combined_output exists")
    p.add_argument("--meyes-base",
                   default=os.environ.get("MEYES_BASE", "/bil/data/meyes"))
    return p.parse_args()


def main():
    args   = parse_args()
    meyes_base = Path(args.meyes_base)
    now    = dt.datetime.now().isoformat(timespec="seconds")

    rows = []
    for name in args.datasets:
        try:
            row = scrape_dataset(
                name, meyes_base, args.mask_percentile, args.species, args.h5ad, now)
        except (FileNotFoundError, ValueError) as e:
            print(f"[skip ] {name}: {e}", file=sys.stderr)
            continue
        rows.append(row)
        print(f"[ok   ] {name}: cells={row.get('n_cells')} "
              f"genes_detected={row.get('n_genes_detected')} "
              f"molecules={row.get('n_genes_imaged')} "
              f"mean_corr={row.get('mean_corr')}")

    if not rows:
        sys.exit("Nothing to write — no datasets produced usable rows.")

    ds = _merge_write(Path(args.out).expanduser(),
                      pd.DataFrame(rows), DATASET_BASE_ORDER, ["dataset"])
    print(f"\nWrote {args.out}: {len(ds)} dataset row(s)")


if __name__ == "__main__":
    main()
