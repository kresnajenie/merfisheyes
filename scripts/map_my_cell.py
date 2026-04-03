import argparse
import json
import logging
import os
import sys
import time
from pathlib import Path

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
    }
}


logger = logging.getLogger(__name__)

# loads combined metadata
def load_inputs(metadata_path, cbg_path):
    metadata_df = pd.read_csv(metadata_path)

    meta_id_col = None
    for candidate in ("EntityID", "id", "cell_id"):
        if candidate in metadata_df.columns:
            meta_id_col = candidate
            break
    if meta_id_col is None:
        raise ValueError(
            "Metadata CSV has no recognized ID column. "
            f"Expected one of 'EntityID', 'id', 'cell_id'. "
            f"Columns present: {list(metadata_df.columns)}"
        )

    metadata_df[meta_id_col] = metadata_df[meta_id_col].astype(str)

    cbg_df = pd.read_csv(cbg_path)
    cbg_df["cell"] = cbg_df["cell"].astype(str)

    metadata_ids = set(metadata_df[meta_id_col])
    cbg_ids = set(cbg_df["cell"])
    
    # failsafe for if you didnt upload a combined file or csv's from different samples
    only_in_metadata = metadata_ids - cbg_ids
    only_in_cbg = cbg_ids - metadata_ids

    if only_in_metadata or only_in_cbg:
        logger.warning(
            "Cell ID mismatch between metadata and cell-by-gene. "
            "%d only in metadata, %d only in cbg. Returning intersection.",
            len(only_in_metadata),
            len(only_in_cbg),
        )

    common_ids = metadata_ids & cbg_ids
    metadata_df = metadata_df[metadata_df[meta_id_col].isin(common_ids)].reset_index(drop=True)
    cbg_df = cbg_df[cbg_df["cell"].isin(common_ids)].reset_index(drop=True)

    return metadata_df, cbg_df


def resolve_reference_dir(cli_arg=None):
    if cli_arg is not None:
        reference_dir = Path(cli_arg).expanduser()
    elif "MERFISHEYES_REFERENCE_DIR" in os.environ:
        reference_dir = Path(os.environ["MERFISHEYES_REFERENCE_DIR"]).expanduser()
    else:
        reference_dir = Path("~/data/merfisheyes-test-data/mapmycells-reference/").expanduser()

    for species_cfg in TAXONOMY_CONFIG.values():
        for key in ("precomputed_stats", "markers", "gene_mapping"):
            expected = reference_dir / species_cfg[key]
            if not expected.exists():
                raise FileNotFoundError(
                    f"Required reference file not found: {expected}"
                )

    return reference_dir


def parse_args():
    parser = argparse.ArgumentParser(
        description="Map cells to a reference taxonomy using MapMyCells."
    )
    parser.add_argument(
        "input_dir",
        help="Path to combined_output folder from combine_slices_v3.py.",
    )
    parser.add_argument(
        "output_dir",
        help="Directory to write output files.",
    )
    parser.add_argument(
        "--reference_dir",
        default=None,
        help="Path to reference data directory (optional).",
    )
    parser.add_argument(
        "--species",
        default="mouse",
        choices=list(TAXONOMY_CONFIG.keys()),
        help="Species for taxonomy mapping (default: mouse).",
    )
    parser.add_argument(
        "--n_processors",
        type=int,
        default=4,
        help="Number of processors to use (default: 4).",
    )
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    return parser.parse_args()

# converts cell id names to ensemblID names for matching
def translate_genes_to_ensembl(cbg_df, gene_mapping_path):
    gene_map = pd.read_csv(gene_mapping_path)
    symbol_to_ensembl = dict(zip(gene_map["gene_symbol"], gene_map["gene_identifier"]))
    
    gene_cols = [c for c in cbg_df.columns if c != "cell"]
    translated = {g: symbol_to_ensembl.get(g) for g in gene_cols}
    
    mapped = {k: v for k, v in translated.items() if v is not None}
    dropped = [k for k, v in translated.items() if v is None]
    
    logger.info("Gene translation: %d mapped, %d dropped (no Ensembl match).", len(mapped), len(dropped))
    if dropped:
        logger.info("Dropped genes: %s", dropped[:10])
    
    cbg_df = cbg_df[["cell"] + list(mapped.keys())].copy()
    cbg_df = cbg_df.rename(columns=mapped)
    return cbg_df

def build_h5ad(cbg_df, metadata_df, output_dir):
    gene_columns = [c for c in cbg_df.columns if c != "cell"]
    X = cbg_df[gene_columns].to_numpy(dtype=np.float32)
    
    # Create globally unique cell index using sample_id + cell_id
    meta_id_col = next(c for c in ("EntityID", "id", "cell_id") if c in metadata_df.columns)
    if "_sample_id" in metadata_df.columns:
        unique_index = metadata_df["_sample_id"].astype(str) + "_" + metadata_df[meta_id_col].astype(str)
    else:
        unique_index = cbg_df["cell"].astype(str)
    
    obs = pd.DataFrame(index=unique_index)
    var = pd.DataFrame(index=gene_columns)
    adata = anndata.AnnData(X=X, obs=obs, var=var)
    query_path = Path(output_dir) / "query.h5ad"
    adata.write_h5ad(query_path)
    return query_path


def run_mapping(validated_h5ad_path, output_dir, reference_dir, species, n_processors):
    from cell_type_mapper.cli.from_specified_markers import FromSpecifiedMarkersRunner
    # expandable for human samples in the future (species = human)
    config = TAXONOMY_CONFIG[species]
    output_dir = Path(output_dir)
    config_dict = {
        "query_path": str(validated_h5ad_path),
        "extended_result_path": str(output_dir / "mapping_output.json"),
        "csv_result_path": str(output_dir / "mapping_output.csv"),
        "drop_level": config["drop_level"],
        "tmp_dir": os.environ.get("TMPDIR", "/tmp"),  # Use SLURM local scratch or /tmp
        "precomputed_stats": {"path": str(reference_dir / config["precomputed_stats"])},
        "query_markers": {"serialized_lookup": str(reference_dir / config["markers"])},
        "type_assignment": {
            "normalization": config["normalization"],
            "n_processors": n_processors,
            "chunk_size": 500,
        },
    }
    runner = FromSpecifiedMarkersRunner(args=[], input_data=config_dict)
    runner.run()
    return output_dir / "mapping_output.csv"


def join_results(mapping_csv_path, metadata_df, output_dir):
    mapping_df = pd.read_csv(mapping_csv_path, comment='#')
    mapping_df["cell_id"] = mapping_df["cell_id"].astype(str)

    meta_id_col = next(c for c in ("EntityID", "id", "cell_id") if c in metadata_df.columns)

    # Build the same compound key used in build_h5ad
    if "_sample_id" in metadata_df.columns:
        metadata_df = metadata_df.copy()
        metadata_df["_join_key"] = metadata_df["_sample_id"].astype(str) + "_" + metadata_df[meta_id_col].astype(str)
    else:
        metadata_df = metadata_df.copy()
        metadata_df["_join_key"] = metadata_df[meta_id_col].astype(str)

    merged = metadata_df.merge(mapping_df, left_on="_join_key", right_on="cell_id", how="left")
    merged = merged.drop(columns=["_join_key"])

    enriched_path = Path(output_dir) / "enriched_metadata.csv"
    merged.to_csv(enriched_path, index=False)
    return enriched_path

def plot_cell_types(enriched_metadata_path, output_dir):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    df = pd.read_csv(enriched_metadata_path)

    if "class_name" not in df.columns or df["class_name"].isna().all():
        logger.warning("class_name column empty or missing, skipping cell type plots.")
        return

    # Pick top 3 most common cell type classes to highlight
    top3 = df["class_name"].value_counts().head(3).index.tolist()

    for cell_type in top3:
        plt.style.use("dark_background")
        fig, ax = plt.subplots(figsize=(14, 14))

        # Plot all cells in grey first
        ax.scatter(df["center_x"], df["center_y"],
                   c="dimgrey", s=0.3, alpha=0.3, rasterized=True)

        # Overlay the highlighted cell type in color
        mask = df["class_name"] == cell_type
        ax.scatter(df.loc[mask, "center_x"], df.loc[mask, "center_y"],
                   c="crimson", s=0.8, alpha=0.8, rasterized=True,
                   label=f"{cell_type} (n={mask.sum():,})")

        ax.set_title(f"{cell_type} — {mask.sum():,} of {len(df):,} cells")
        ax.set_aspect("equal")
        ax.legend(loc='upper right', fontsize=8, framealpha=0.3)

        plt.tight_layout()
        safe_name = cell_type.replace(" ", "_").replace("/", "-")
        out_path = Path(output_dir) / f"check_celltype_{safe_name}.png"
        plt.savefig(out_path, dpi=150)
        plt.close()
        logger.info("Saved %s", out_path)

def print_summary(enriched_metadata_path, query_h5ad_path, original_cbg_gene_count):
    enriched_df = pd.read_csv(enriched_metadata_path)
    validated_adata = anndata.read_h5ad(query_h5ad_path, backed="r")
    n_genes = len(validated_adata.var.index)
    gene_overlap_pct = n_genes / original_cbg_gene_count * 100

    print(f"Gene overlap: {n_genes}/{original_cbg_gene_count} ({gene_overlap_pct:.1f}%)")
    print(f"Total cells mapped: {len(enriched_df)}")

    cluster_col = "cluster_label" if "cluster_label" in enriched_df.columns else enriched_df.columns[-1]
    print(f"\nTop 5 assigned cluster labels ({cluster_col}):")
    print(enriched_df[cluster_col].value_counts().head(5).to_string())


def timed(label):
    """Context manager to time a step and print duration."""
    class Timer:
        def __enter__(self):
            self.start = time.time()
            print(f"\n[STEP] {label}...")
            return self
        def __exit__(self, *args):
            elapsed = time.time() - self.start
            m, s = divmod(elapsed, 60)
            print(f"[DONE] {label} — {int(m)}m {s:.1f}s")
    return Timer()


if __name__ == "__main__":
    args = parse_args()
    total_start = time.time()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        stream=sys.stdout,
    )

    # Log thread settings
    print(f"MKL_NUM_THREADS={os.environ.get('MKL_NUM_THREADS', '(not set)')}")
    print(f"OMP_NUM_THREADS={os.environ.get('OMP_NUM_THREADS', '(not set)')}")
    print(f"n_processors={args.n_processors}")
    print(f"TMPDIR={os.environ.get('TMPDIR', '(not set)')}")

    input_dir = Path(args.input_dir).expanduser()
    metadata_path = input_dir / "cell_metadata.csv"
    cbg_path = input_dir / "cell_by_gene.csv"

    if not metadata_path.exists():
        raise FileNotFoundError(f"cell_metadata.csv not found in {input_dir}")
    if not cbg_path.exists():
        raise FileNotFoundError(f"cell_by_gene.csv not found in {input_dir}")

    reference_dir = resolve_reference_dir(args.reference_dir)
    output_dir = Path(args.output_dir).expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)

    with timed("1/6 Load inputs"):
        metadata_df, cbg_df = load_inputs(metadata_path, cbg_path)
        original_gene_count = len(cbg_df.columns) - 1
        print(f"     Cells: {len(metadata_df):,}  Genes: {original_gene_count:,}")

    with timed("2/6 Translate genes to Ensembl"):
        gene_mapping_path = reference_dir / TAXONOMY_CONFIG[args.species]["gene_mapping"]
        cbg_df = translate_genes_to_ensembl(cbg_df, gene_mapping_path)
        print(f"     Genes after translation: {len(cbg_df.columns) - 1:,}")

    with timed("3/6 Build query H5AD"):
        query_h5ad_path = build_h5ad(cbg_df, metadata_df, output_dir)
        print(f"     Written to: {query_h5ad_path}")

    with timed("4/6 Run MapMyCells mapping"):
        mapping_csv_path = run_mapping(
            query_h5ad_path, output_dir, reference_dir, args.species, args.n_processors
        )

    with timed("5/6 Join results"):
        enriched_metadata_path = join_results(mapping_csv_path, metadata_df, output_dir)

    with timed("6/6 Plot cell types"):
        plot_cell_types(enriched_metadata_path, output_dir)

    print_summary(enriched_metadata_path, query_h5ad_path, original_gene_count)

    total_elapsed = time.time() - total_start
    m, s = divmod(total_elapsed, 60)
    print(f"\n{'='*40}")
    print(f"  Total pipeline time: {int(m)}m {s:.1f}s")
    print(f"{'='*40}")
