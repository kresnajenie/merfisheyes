#!/usr/bin/env python3
"""Map cells to a reference taxonomy using MapMyCells.

Emits exactly one artifact: ``mapping_output.csv``, shaped the way
``process_spatial_data.py --mmc-csv`` consumes it — one row per cell, in the
input's cell order, with a ``cell_id`` column plus the taxonomy label columns.
Every other column in that CSV becomes an obs column downstream.

Input is an ``.h5ad`` with raw counts (in ``X``, or ``obsm['X_raw']`` when ``X``
has been normalized).

Note: the intermediate query file is unavoidable — MapMyCells'
``FromSpecifiedMarkersRunner`` takes a *path*, not an in-memory object — but it
is written to a temp dir and removed, never into ``output_dir``.
"""
import argparse
import logging
import multiprocessing
import os
import sys
import tempfile
from pathlib import Path

# Prevent thread oversubscription — must be set before importing NumPy/MKL
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import anndata
import numpy as np
import pandas as pd


TAXONOMY_CONFIG = {
    "mouse": {
        "precomputed_stats": "mouse/precomputed_stats_ABC_revision_230821.h5",
        "markers": "mouse/mouse_markers_230821.json",
        "gene_mapping": "mouse/gene.csv",
        "drop_level": "CCN20230722_SUPT",
        "normalization": "raw",
    },
    "human": {
        "precomputed_stats": "human/precomputed_stats.siletti.training.h5",
        "markers": "human/query_markers.n10.20240221800.json",
        "gene_mapping": "human/gene.csv",
        "drop_level": "CCN202210140_SUPC",
        "normalization": "raw",
    },
}

logger = logging.getLogger(__name__)


def resolve_reference_dir(species, cli_arg=None):
    """Locate the reference data and validate ONLY the requested species.

    (Validating every species would make a mouse-only reference set fail with a
    confusing complaint about missing human files.)
    """
    if cli_arg is not None:
        reference_dir = Path(cli_arg).expanduser()
    elif "MERFISHEYES_REFERENCE_DIR" in os.environ:
        reference_dir = Path(os.environ["MERFISHEYES_REFERENCE_DIR"]).expanduser()
    else:
        reference_dir = Path("~/data/merfisheyes-test-data/mapmycells-reference/").expanduser()

    cfg = TAXONOMY_CONFIG[species]
    missing = [
        str(reference_dir / cfg[key])
        for key in ("precomputed_stats", "markers", "gene_mapping")
        if not (reference_dir / cfg[key]).exists()
    ]
    if missing:
        raise FileNotFoundError(
            f"Missing {species} reference file(s) under {reference_dir}:\n  "
            + "\n  ".join(missing)
        )

    return reference_dir


def select_raw_matrix(adata):
    """Return the raw-count matrix, preferring X and falling back to obsm['X_raw']."""
    import scipy.sparse as sp

    X = adata.X
    sample = X[:10].toarray() if sp.issparse(X) else np.asarray(X[:10])
    is_raw = np.issubdtype(X.dtype, np.integer) or np.all(sample == np.floor(sample))

    if is_raw:
        return X

    if "X_raw" not in adata.obsm:
        raise ValueError(
            "adata.X appears normalized and adata.obsm['X_raw'] is missing. "
            "Provide an h5ad with raw counts in X or store them in obsm['X_raw']."
        )

    return adata.obsm["X_raw"]


def build_query_h5ad(input_h5ad, gene_mapping_path, tmp_dir):
    """Write a query h5ad whose var index is Ensembl IDs.

    Stays in AnnData/sparse throughout — the matrix is never densified and never
    round-tripped through a per-gene pandas DataFrame.
    """
    adata = anndata.read_h5ad(input_h5ad)
    X = select_raw_matrix(adata)

    gene_map = pd.read_csv(gene_mapping_path)
    symbol_to_ensembl = dict(zip(gene_map["gene_symbol"], gene_map["gene_identifier"]))

    var_names = adata.var.index.astype(str)
    keep = np.array([s in symbol_to_ensembl for s in var_names], dtype=bool)
    n_total = int(keep.size)
    n_mapped = int(keep.sum())

    if n_mapped == 0:
        raise ValueError(
            "No genes could be translated to Ensembl IDs. "
            f"Checked {n_total} genes against {gene_mapping_path}."
        )

    logger.info(
        "Gene translation: %d/%d mapped to Ensembl (%d dropped).",
        n_mapped, n_total, n_total - n_mapped,
    )

    ensembl_ids = [symbol_to_ensembl[s] for s in var_names[keep]]

    query = anndata.AnnData(
        X=X[:, keep],
        obs=pd.DataFrame(index=adata.obs.index.astype(str)),
        var=pd.DataFrame(index=pd.Index(ensembl_ids, name=None)),
    )

    query_path = Path(tmp_dir) / "query.h5ad"
    query.write_h5ad(query_path)

    return query_path, n_mapped, n_total


def run_mapping(query_path, output_dir, reference_dir, species, n_processors, flatten=False):
    from cell_type_mapper.cli.from_specified_markers import FromSpecifiedMarkersRunner

    cfg = TAXONOMY_CONFIG[species]
    csv_result_path = Path(output_dir) / "mapping_output.csv"

    # extended_result_path is required by the runner even though we discard it.
    json_tmp = tempfile.NamedTemporaryFile(suffix=".json", delete=False)
    json_tmp.close()

    # Use the SLURM allocation when present, else physical RAM.
    slurm_mem = os.environ.get("SLURM_MEM_PER_NODE")
    if slurm_mem:
        total_ram_gb = int(slurm_mem) // 1024  # SLURM reports MB
    else:
        total_ram_gb = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES") // (1024 ** 3)
    max_gb = max(1, total_ram_gb // 2)

    logger.info(
        "System: %d cores, %dGB RAM → n_processors=%d, max_gb=%d",
        multiprocessing.cpu_count(), total_ram_gb, n_processors, max_gb,
    )

    if flatten:
        bootstrap_iteration, bootstrap_factor = 1, 1.0
        logger.info("Flatten mode: skipping hierarchy, bootstrap_iteration=1")
    else:
        bootstrap_iteration, bootstrap_factor = 100, 0.9

    try:
        runner = FromSpecifiedMarkersRunner(
            args=[],
            input_data={
                "query_path": str(query_path),
                "extended_result_path": json_tmp.name,
                "csv_result_path": str(csv_result_path),
                "drop_level": cfg["drop_level"],
                "flatten": flatten,
                "tmp_dir": os.environ.get("TMPDIR", "/tmp"),
                "max_gb": max_gb,
                "precomputed_stats": {"path": str(reference_dir / cfg["precomputed_stats"])},
                "query_markers": {"serialized_lookup": str(reference_dir / cfg["markers"])},
                "type_assignment": {
                    "normalization": cfg["normalization"],
                    "n_processors": n_processors,
                    "chunk_size": 2000,
                    "bootstrap_iteration": bootstrap_iteration,
                    "bootstrap_factor": bootstrap_factor,
                },
            },
        )
        runner.run()
    finally:
        os.unlink(json_tmp.name)

    return csv_result_path


def verify_output(csv_path, expected_cells):
    """process_spatial_data.py requires exactly one CSV row per cell — fail here
    rather than deep inside the chunking stage."""
    df = pd.read_csv(csv_path, comment="#")
    if len(df) != expected_cells:
        raise ValueError(
            f"{csv_path.name} has {len(df):,} rows but the input has "
            f"{expected_cells:,} cells — downstream requires an exact match."
        )
    logger.info("Wrote %s (%d rows, %d columns)", csv_path, len(df), len(df.columns))


def parse_args():
    parser = argparse.ArgumentParser(
        description="Map cells to a reference taxonomy using MapMyCells."
    )
    parser.add_argument("input_h5ad", help="Input .h5ad file with raw counts.")
    parser.add_argument("output_dir", help="Directory to write mapping_output.csv into.")
    parser.add_argument("--reference_dir", default=None,
                        help="Reference data directory (else $MERFISHEYES_REFERENCE_DIR).")
    parser.add_argument("--species", default="mouse", choices=list(TAXONOMY_CONFIG.keys()),
                        help="Species for taxonomy mapping (default: mouse).")
    available = multiprocessing.cpu_count()
    parser.add_argument("--n_processors", type=int, default=max(1, available // 2),
                        help=f"Processors to use (default: half of {available}).")
    parser.add_argument("--flatten", action="store_true",
                        help="Map straight to leaf nodes (also sets bootstrap_iteration=1; much faster).")
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    return parser.parse_args()


def main():
    args = parse_args()
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        stream=sys.stdout,
    )

    input_h5ad = Path(args.input_h5ad).expanduser()
    if not (input_h5ad.is_file() and input_h5ad.suffix == ".h5ad"):
        raise SystemExit(f"Input must be an .h5ad file: {input_h5ad}")

    reference_dir = resolve_reference_dir(args.species, args.reference_dir)
    output_dir = Path(args.output_dir).expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)

    gene_mapping_path = reference_dir / TAXONOMY_CONFIG[args.species]["gene_mapping"]

    # The query h5ad lives only for the duration of the run.
    with tempfile.TemporaryDirectory(prefix="mmc_query_") as tmp_dir:
        query_path, n_mapped, n_total = build_query_h5ad(
            input_h5ad, gene_mapping_path, tmp_dir
        )
        n_cells = anndata.read_h5ad(query_path, backed="r").shape[0]
        csv_path = run_mapping(
            query_path, output_dir, reference_dir, args.species,
            args.n_processors, flatten=args.flatten,
        )

    verify_output(csv_path, n_cells)
    print(f"Gene overlap: {n_mapped}/{n_total} ({n_mapped / n_total * 100:.1f}%)")
    print(f"Cells mapped: {n_cells:,}")
    print(f"Output: {csv_path}")


if __name__ == "__main__":
    main()
