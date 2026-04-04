#!/usr/bin/env python3
"""Generate synthetic single-molecule spatial transcriptomics datasets.

Produces parquet and CSV files at configurable scales.
48 combinations by default (12 molecule counts × 2 gene counts × 2 filetypes).

Distribution:
  - 50% of molecules → "Mal" gene
  - Of remaining 50%: 40% assigned to real genes, 60% unassigned (control probes)
  - So overall: 50% Mal, 20% other genes, 30% control probes

Usage:
    python generate_single_molecule.py
    python generate_single_molecule.py --molecules 10000000 30000000 --genes 1000 --filetypes parquet
    python generate_single_molecule.py --dry-run
"""

import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import (
    SINGLE_MOLECULE_SIZES, SINGLE_MOLECULE_GENES, SINGLE_MOLECULE_FILETYPES,
    MAL_FRACTION, ASSIGNED_FRACTION, CONTROL_PROBES, MAX_FILE_SIZE_GB,
    estimate_parquet_size_mb, estimate_csv_size_mb,
)

CHUNK_SIZE = 1_000_000  # molecules per chunk during generation


def make_gene_pool(n_genes: int) -> list[str]:
    """Build the pool of real gene names (including Mal)."""
    genes = ["Mal"]
    for i in range(1, n_genes):
        genes.append(f"gene_{i}")
    return genes


def make_unassigned_pool() -> list[str]:
    return list(CONTROL_PROBES)


def estimate_size(filetype: str, n_molecules: int) -> float:
    if filetype == "parquet":
        return estimate_parquet_size_mb(n_molecules)
    return estimate_csv_size_mb(n_molecules)


def fmt_size(mb: float) -> str:
    if mb < 1024:
        return f"{mb:.1f} MB"
    return f"{mb / 1024:.2f} GB"


def output_path_for(base_dir: str, filetype: str, n_mol: int, n_genes: int) -> str:
    tag = f"sm_{n_mol}m_{n_genes}g"
    return os.path.join(base_dir, filetype, f"{tag}.{filetype}")


# ---------------------------------------------------------------------------
# Chunk generator
# ---------------------------------------------------------------------------
def generate_chunk(
    chunk_size: int,
    gene_pool: list[str],
    unassigned_pool: list[str],
    mal_frac: float,
    assigned_frac: float,
    rng: np.random.Generator,
) -> dict[str, np.ndarray]:
    """Generate one chunk of molecule data.

    Returns dict with keys: feature_name, x_location, y_location, z_location
    """
    n = chunk_size
    # Determine how many molecules go to each category
    n_mal = int(n * mal_frac)
    n_remaining = n - n_mal
    n_assigned_other = int(n_remaining * assigned_frac)
    n_unassigned = n - n_mal - n_assigned_other

    # Gene assignments
    genes = np.empty(n, dtype=object)
    genes[:n_mal] = "Mal"
    # Other assigned genes (excluding Mal)
    other_genes = [g for g in gene_pool if g != "Mal"]
    if other_genes:
        genes[n_mal:n_mal + n_assigned_other] = rng.choice(other_genes, size=n_assigned_other)
    else:
        genes[n_mal:n_mal + n_assigned_other] = "Mal"
    # Unassigned
    genes[n_mal + n_assigned_other:] = rng.choice(unassigned_pool, size=n_unassigned)

    # Shuffle so categories are mixed
    perm = rng.permutation(n)
    genes = genes[perm]

    # Coordinates
    x = rng.uniform(-5000, 5000, size=n).astype(np.float32)
    y = rng.uniform(-5000, 5000, size=n).astype(np.float32)
    z = rng.uniform(-50, 50, size=n).astype(np.float32)

    return {
        "feature_name": genes,
        "x_location": x,
        "y_location": y,
        "z_location": z,
    }


# ---------------------------------------------------------------------------
# Parquet writer
# ---------------------------------------------------------------------------
def generate_parquet(
    n_molecules: int, n_genes: int, path: str,
    mal_frac: float, assigned_frac: float,
):
    import pyarrow as pa
    import pyarrow.parquet as pq

    rng = np.random.default_rng(42)
    gene_pool = make_gene_pool(n_genes)
    unassigned_pool = make_unassigned_pool()

    os.makedirs(os.path.dirname(path), exist_ok=True)

    schema = pa.schema([
        ("feature_name", pa.string()),
        ("x_location", pa.float32()),
        ("y_location", pa.float32()),
        ("z_location", pa.float32()),
    ])

    writer = pq.ParquetWriter(
        path, schema,
        compression="snappy",
        version="1.0",  # hyparquet compatibility
    )

    written = 0
    while written < n_molecules:
        chunk_n = min(CHUNK_SIZE, n_molecules - written)
        chunk = generate_chunk(chunk_n, gene_pool, unassigned_pool, mal_frac, assigned_frac, rng)
        table = pa.table({
            "feature_name": pa.array(chunk["feature_name"].tolist(), type=pa.string()),
            "x_location": pa.array(chunk["x_location"]),
            "y_location": pa.array(chunk["y_location"]),
            "z_location": pa.array(chunk["z_location"]),
        }, schema=schema)
        writer.write_table(table)
        written += chunk_n
        if written % (CHUNK_SIZE * 10) == 0 or written == n_molecules:
            pct = written / n_molecules * 100
            print(f"    {written:>15,} / {n_molecules:,}  ({pct:.0f}%)")

    writer.close()


# ---------------------------------------------------------------------------
# CSV writer
# ---------------------------------------------------------------------------
def generate_csv(
    n_molecules: int, n_genes: int, path: str,
    mal_frac: float, assigned_frac: float,
):
    rng = np.random.default_rng(42)
    gene_pool = make_gene_pool(n_genes)
    unassigned_pool = make_unassigned_pool()

    os.makedirs(os.path.dirname(path), exist_ok=True)

    written = 0
    with open(path, "w") as f:
        f.write("feature_name,x_location,y_location,z_location\n")
        while written < n_molecules:
            chunk_n = min(CHUNK_SIZE, n_molecules - written)
            chunk = generate_chunk(chunk_n, gene_pool, unassigned_pool, mal_frac, assigned_frac, rng)
            lines = []
            for i in range(chunk_n):
                lines.append(
                    f"{chunk['feature_name'][i]},"
                    f"{chunk['x_location'][i]:.4f},"
                    f"{chunk['y_location'][i]:.4f},"
                    f"{chunk['z_location'][i]:.4f}"
                )
            f.write("\n".join(lines) + "\n")
            written += chunk_n
            if written % (CHUNK_SIZE * 10) == 0 or written == n_molecules:
                pct = written / n_molecules * 100
                print(f"    {written:>15,} / {n_molecules:,}  ({pct:.0f}%)")


GENERATORS = {
    "parquet": generate_parquet,
    "csv": generate_csv,
}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Generate synthetic single-molecule data")
    parser.add_argument("--output-dir", default=os.path.join(DEFAULT_SYNTH_DIR, "single_molecule"),
                        help="Output directory")
    parser.add_argument("--molecules", type=int, nargs="+", default=SINGLE_MOLECULE_SIZES,
                        help="Molecule counts to generate")
    parser.add_argument("--genes", type=int, nargs="+", default=SINGLE_MOLECULE_GENES,
                        help="Gene counts to generate")
    parser.add_argument("--filetypes", nargs="+", default=SINGLE_MOLECULE_FILETYPES,
                        choices=SINGLE_MOLECULE_FILETYPES, help="File types to generate")
    parser.add_argument("--mal-fraction", type=float, default=MAL_FRACTION,
                        help="Fraction of molecules assigned to Mal (default 0.5)")
    parser.add_argument("--assigned-fraction", type=float, default=ASSIGNED_FRACTION,
                        help="Fraction of non-Mal molecules assigned to real genes (default 0.4)")
    parser.add_argument("--max-file-size-gb", type=float, default=MAX_FILE_SIZE_GB,
                        help="Skip files estimated larger than this (GB)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Show plan without generating files")
    args = parser.parse_args()

    combinations = [
        (nm, ng, ft)
        for nm in args.molecules
        for ng in args.genes
        for ft in args.filetypes
    ]
    print(f"Benchmark single-molecule data generator")
    print(f"  Combinations: {len(combinations)} "
          f"({len(args.molecules)} molecule sizes × {len(args.genes)} gene counts × {len(args.filetypes)} filetypes)")
    print(f"  Mal fraction: {args.mal_fraction:.0%}")
    print(f"  Assigned fraction (of non-Mal): {args.assigned_fraction:.0%}")
    print(f"  Overall: {args.mal_fraction:.0%} Mal, "
          f"{(1-args.mal_fraction)*args.assigned_fraction:.0%} other genes, "
          f"{(1-args.mal_fraction)*(1-args.assigned_fraction):.0%} unassigned")
    print()

    total_mb = 0.0
    skipped = 0
    plan = []
    for nm, ng, ft in combinations:
        est_mb = estimate_size(ft, nm)
        path = output_path_for(args.output_dir, ft, nm, ng)
        skip = est_mb > args.max_file_size_gb * 1024
        plan.append((nm, ng, ft, est_mb, path, skip))
        if skip:
            skipped += 1
        else:
            total_mb += est_mb

    print(f"{'Molecules':>15}  {'Genes':>8}  {'Type':>10}  {'Est Size':>12}  {'Status'}")
    print("-" * 70)
    for nm, ng, ft, est_mb, path, skip in plan:
        status = "SKIP (too large)" if skip else "generate"
        print(f"{nm:>15,}  {ng:>8,}  {ft:>10}  {fmt_size(est_mb):>12}  {status}")

    print()
    print(f"Total to generate: {fmt_size(total_mb)}  |  Skipped: {skipped}")

    if args.dry_run:
        print("\n(dry run – no files created)")
        return

    generated = 0
    for nm, ng, ft, est_mb, path, skip in plan:
        if skip:
            continue
        if os.path.exists(path):
            print(f"\n[SKIP] Already exists: {path}")
            generated += 1
            continue

        print(f"\n[{generated+1}/{len(plan)-skipped}] Generating {ft} — "
              f"{nm:,} molecules × {ng:,} genes (est {fmt_size(est_mb)})")
        t0 = time.time()
        try:
            GENERATORS[ft](nm, ng, path, args.mal_fraction, args.assigned_fraction)
            elapsed = time.time() - t0
            actual_mb = os.path.getsize(path) / (1024 * 1024)
            print(f"  Done in {elapsed:.1f}s — actual size: {fmt_size(actual_mb)}")
            generated += 1
        except Exception as e:
            print(f"  ERROR: {e}")

    # Write manifest
    manifest = {
        "generated": time.strftime("%Y-%m-%d %H:%M:%S"),
        "mal_fraction": args.mal_fraction,
        "assigned_fraction": args.assigned_fraction,
        "datasets": [],
    }
    for nm, ng, ft, est_mb, path, skip in plan:
        if skip or not os.path.exists(path):
            continue
        actual_mb = os.path.getsize(path) / (1024 * 1024)
        manifest["datasets"].append({
            "n_molecules": nm,
            "n_genes": ng,
            "filetype": ft,
            "path": path,
            "file_size_mb": round(actual_mb, 2),
        })

    manifest_path = os.path.join(args.output_dir, "manifest_single_molecule.json")
    os.makedirs(os.path.dirname(manifest_path), exist_ok=True)
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"\nManifest written to {manifest_path}")
    print(f"Generated {generated} files.")


if __name__ == "__main__":
    DEFAULT_SYNTH_DIR = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
        "benchmark_data", "synthetic"
    )
    main()
