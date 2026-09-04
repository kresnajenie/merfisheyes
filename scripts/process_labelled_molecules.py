#!/usr/bin/env python3
"""
Labelled single-molecule data to the chunked binary format.

Points are molecules, each carrying N categorical labels and NO expression
matrix — the third dataset type alongside single_cell and single_molecule.
Output is the same layout the single-cell viewer already reads, minus `expr/`:

    output/
    ├── manifest.json          has_expression: false, normalized: false
    ├── coords/spatial.bin.gz  raw coordinates (µm), float32
    ├── obs/
    │   ├── metadata.json
    │   └── {column}.bin.gz    dictionary-encoded categorical, bin-v1
    └── palettes/{column}.json

Coordinates stay RAW so molecules co-register with cell meshes, which are
already in µm.

Usage:
    python process_labelled_molecules.py molecules.parquet output/ \
        --category gene=feature_name \
        --category domain_anno=rna_domain_anno \
        --category domain_id=rna_domain_id \
        --category cell=cell_id --label cell=cell_name \
        --drop-unassigned cell_id \
        --drop-control-genes feature_name
"""

import argparse
import gzip
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

# Single-source the binary layout and the filter rules — these mirror
# lib/utils/obs-binary.ts and lib/utils/gene-filters.ts, and must not drift.
from process_spatial_data import (  # noqa: E402
    OBS_BINARY_FORMAT,
    build_obs_palette,
    encode_obs_binary,
    fmt_elapsed,
    peak_rss_gb,
    write_coordinate_binary,
)
from process_single_molecule import (  # noqa: E402
    compute_unassigned_mask,
    should_filter_gene,
)

_t0 = None


def log(msg):
    el = f"[{fmt_elapsed(time.perf_counter() - _t0)}]" if _t0 else ""
    print(f"{el} {msg}", flush=True)


def parse_kv(pairs, flag):
    """['gene=feature_name', ...] -> {'gene': 'feature_name'}"""
    out = {}
    for p in pairs or []:
        if "=" not in p:
            raise SystemExit(f"{flag} expects NAME=COLUMN, got {p!r}")
        name, col = p.split("=", 1)
        out[name.strip()] = col.strip()
    return out


def disambiguate(keys, labels):
    """Display labels for a keyed column, suffixed where labels collide.

    Cell names repeat across cells (two polar bodies are both 'pb'), so the
    menu keys on the id and shows 'pb (1)' / 'pb (3)' rather than merging them.
    """
    df = pd.DataFrame({"key": keys, "label": labels}).drop_duplicates()
    dup = df["label"].duplicated(keep=False)
    df["display"] = np.where(dup, df["label"] + " (" + df["key"].astype(str) + ")",
                             df["label"])
    return dict(zip(df["key"], df["display"]))


def main():
    global _t0
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", help="Input parquet file")
    ap.add_argument("output_folder")
    ap.add_argument("--category", action="append", metavar="NAME=COLUMN",
                    help="A categorical column to expose. Repeatable.")
    ap.add_argument("--label", action="append", metavar="NAME=COLUMN",
                    help="Human-readable label column for a --category keyed on "
                         "an id. Repeatable.")
    ap.add_argument("--x", default="x_location")
    ap.add_argument("--y", default="y_location")
    ap.add_argument("--z", default="z_location", help="Omit or pass '' for 2D")
    ap.add_argument("--drop-unassigned", action="append", metavar="COLUMN",
                    help="Drop rows unassigned in COLUMN (numeric <= 0, NaN, or "
                         "'unassigned'/'none'/…). Repeatable.")
    ap.add_argument("--drop-control-genes", action="append", metavar="COLUMN",
                    help="Drop rows whose COLUMN value is a control/blank probe. "
                         "Repeatable.")
    ap.add_argument("--name", help="Dataset name (default: input stem)")
    args = ap.parse_args()

    _t0 = time.perf_counter()
    categories = parse_kv(args.category, "--category")
    labels = parse_kv(args.label, "--label")
    if not categories:
        raise SystemExit("At least one --category is required")

    src = Path(args.input)
    out = Path(args.output_folder)
    name = args.name or src.stem

    coord_cols = [args.x, args.y] + ([args.z] if args.z else [])
    needed = list(dict.fromkeys(coord_cols + list(categories.values())
                                + list(labels.values())
                                + (args.drop_unassigned or [])
                                + (args.drop_control_genes or [])))

    log(f"=== STEP 1: Reading {src.name} ===")
    import pyarrow.parquet as pq
    table = pq.read_table(src, columns=needed)
    n_in = table.num_rows
    log(f"  {n_in:,} rows, {len(needed)} columns")

    # ── Row filter ────────────────────────────────────────────────────────
    log("=== STEP 2: Filtering rows ===")
    keep = np.ones(n_in, dtype=bool)
    for col in args.drop_unassigned or []:
        mask = ~compute_unassigned_mask(table[col].to_numpy(zero_copy_only=False))
        log(f"  --drop-unassigned {col}: dropping {int((~mask).sum()):,}")
        keep &= mask
    for col in args.drop_control_genes or []:
        vals = np.asarray(table[col].to_pylist())
        uniq = pd.unique(vals)
        bad = {u for u in uniq if should_filter_gene(str(u))}
        mask = ~np.isin(vals, list(bad)) if bad else np.ones(n_in, dtype=bool)
        log(f"  --drop-control-genes {col}: dropping {int((~mask).sum()):,} "
            f"on {len(bad)} control probe(s)")
        keep &= mask
    n = int(keep.sum())
    log(f"  {n_in:,} -> {n:,} molecules")
    if n == 0:
        raise SystemExit("Every row was filtered out")

    # ── Coordinates ───────────────────────────────────────────────────────
    log("=== STEP 3: Writing coordinates ===")
    coords = np.stack([table[c].to_numpy(zero_copy_only=False)[keep]
                       for c in coord_cols], axis=1).astype(np.float32)
    finite = np.isfinite(coords).all(axis=1)
    if not finite.all():
        log(f"  dropping {int((~finite).sum()):,} rows with non-finite coordinates")
        coords = coords[finite]
        idx = np.flatnonzero(keep)[finite]
        keep = np.zeros(n_in, dtype=bool)
        keep[idx] = True
        n = len(coords)
    write_coordinate_binary(coords, out / "coords" / "spatial.bin.gz")
    lo, hi = coords.min(axis=0), coords.max(axis=0)
    log(f"  extent {' × '.join(f'{a:.1f}..{b:.1f}' for a, b in zip(lo, hi))} (raw units)")

    # ── Categorical columns ───────────────────────────────────────────────
    log(f"=== STEP 4: Writing {len(categories)} obs column(s) ===")
    obs_dir, pal_dir = out / "obs", out / "palettes"
    obs_meta = {}
    for i, (cname, srccol) in enumerate(categories.items(), 1):
        vals = np.asarray(table[srccol].to_pylist(), dtype=object)[keep]
        # A column keyed on an id renders via its label column, with colliding
        # labels suffixed so distinct entities stay distinct in the menu.
        if cname in labels:
            lab = np.asarray(table[labels[cname]].to_pylist(), dtype=object)[keep]
            display = disambiguate(vals, lab)
            vals = np.asarray([display[k] for k in vals], dtype=object)

        payload, n_unique, cats = encode_obs_binary(vals, categorical=True)
        obs_dir.mkdir(parents=True, exist_ok=True)
        with gzip.open(obs_dir / f"{cname}.bin.gz", "wb") as f:
            f.write(payload)

        palette = build_obs_palette(cats, None, None)
        pal_dir.mkdir(parents=True, exist_ok=True)
        with open(pal_dir / f"{cname}.json", "w") as f:
            json.dump(palette, f, indent=2)

        obs_meta[cname] = {"type": "categorical", "unique_values": n_unique,
                           "format": OBS_BINARY_FORMAT}
        log(f"  [{i}/{len(categories)}] {cname} <- {srccol}: {n_unique} unique")

    with open(obs_dir / "metadata.json", "w") as f:
        json.dump(obs_meta, f, indent=2)

    # ── Manifest ──────────────────────────────────────────────────────────
    log("=== STEP 5: Manifest ===")
    manifest = {
        "version": "1.0",
        "normalized": False,
        # No expression matrix: the adapter must skip expr/index.json.
        "has_expression": False,
        "dataset_id": "local_dataset",
        "name": name,
        "type": "labelled_single_molecule",
        "statistics": {
            "total_cells": n,
            "total_genes": 0,
            "spatial_dimensions": coords.shape[1],
            "available_embeddings": [],
            "cluster_count": len(categories),
        },
        "processing": {
            "spatial_scaling_factor": 1.0,
            "created_by": "process_labelled_molecules.py",
        },
        "files": {
            "coordinates": ["spatial"],
            "expression_chunks": 0,
            "observation_columns": list(obs_meta.keys()),
            "de_stats": [],
        },
    }
    with open(out / "manifest.json", "w") as f:
        json.dump(manifest, f, indent=2)

    size = sum(f.stat().st_size for f in out.rglob("*") if f.is_file())
    log("=" * 60)
    log(f"  Molecules:  {n:,}")
    cols = ", ".join(f"{k} ({v['unique_values']})" for k, v in obs_meta.items())
    log(f"  Columns:    {cols}")
    log(f"  Output:     {out}  ({size / 1e6:.1f} MB)")
    log(f"  Total time: {fmt_elapsed(time.perf_counter() - _t0)}   Peak RSS: {peak_rss_gb():.2f} GB")
    log("=" * 60)


if __name__ == "__main__":
    main()
