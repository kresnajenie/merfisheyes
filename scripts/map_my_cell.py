import argparse
import json
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
        help="Path to combined_output folder from combine_slices_v3.py, or a .h5ad file.",
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
    available_cores = multiprocessing.cpu_count()
    default_cores = max(1, available_cores // 2)
    parser.add_argument(
        "--n_processors",
        type=int,
        default=default_cores,
        help=f"Number of processors to use (default: {default_cores}, half of {available_cores} available).",
    )
    parser.add_argument(
        "--flatten",
        action="store_true",
        help="Map directly to leaf nodes, skipping hierarchical traversal. "
             "Also sets bootstrap_iteration=1 for much faster mapping. "
             "(hierarchical method only)",
    )
    parser.add_argument(
        "--method",
        default="hierarchical",
        choices=["hierarchical", "correlation"],
        help="Mapping method. 'hierarchical' uses Allen's bootstrapped marker "
             "taxonomy traversal (FromSpecifiedMarkersRunner). 'correlation' "
             "Pearson-correlates each cell against reference cluster mean "
             "expression profiles from the precomputed_stats file.",
    )
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    return parser.parse_args()


def load_h5ad(h5ad_path):
    """Load an h5ad file and return (metadata_df, cbg_df) matching load_inputs() shape."""
    import scipy.sparse as sp

    adata = anndata.read_h5ad(h5ad_path)

    # Decide which matrix holds raw counts
    X = adata.X
    sample = X[:10].toarray() if sp.issparse(X) else X[:10]
    is_raw = np.issubdtype(X.dtype, np.integer) or np.all(sample == np.floor(sample))

    if not is_raw:
        if "X_raw" not in adata.obsm:
            raise ValueError(
                "adata.X appears normalized and adata.obsm['X_raw'] is missing. "
                "Provide an h5ad with raw counts in X or store them in obsm['X_raw']."
            )
        X = adata.obsm["X_raw"]

    dense = X.toarray() if sp.issparse(X) else np.asarray(X)

    cell_ids = adata.obs.index.astype(str)
    gene_names = adata.var.index.tolist()

    cbg_df = pd.DataFrame(dense, columns=gene_names)
    cbg_df.insert(0, "cell", cell_ids.values)

    metadata_df = pd.DataFrame({"cell_id": cell_ids.values})

    if "X_spatial" in adata.obsm:
        spatial = adata.obsm["X_spatial"]
        if spatial.shape[1] >= 2:
            # Pick the two columns with the largest range as x/y for plotting.
            # Some h5ad files store a narrow-range slice index in column 0,
            # which produces a flat line if used as a spatial axis.
            ranges = np.ptp(spatial, axis=0)
            top2 = np.argsort(ranges)[-2:][::-1]  # two widest-range columns
            metadata_df["center_x"] = spatial[:, top2[0]]
            metadata_df["center_y"] = spatial[:, top2[1]]
            remaining = [i for i in range(spatial.shape[1]) if i not in top2]
            if remaining:
                metadata_df["center_z"] = spatial[:, remaining[0]]

    return metadata_df, cbg_df


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


def run_mapping(validated_h5ad_path, output_dir, reference_dir, species,
                n_processors, flatten=False, method="hierarchical"):
    """Dispatch to the chosen mapping method. Writes mapping_output.csv."""
    if method == "hierarchical":
        return run_mapping_hierarchical(
            validated_h5ad_path, output_dir, reference_dir, species,
            n_processors, flatten=flatten,
        )
    elif method == "correlation":
        return run_mapping_correlation(
            validated_h5ad_path, output_dir, reference_dir, species
        )
    else:
        raise ValueError(f"Unknown mapping method: {method}")


def run_mapping_hierarchical(validated_h5ad_path, output_dir, reference_dir,
                              species, n_processors, flatten=False):
    from cell_type_mapper.cli.from_specified_markers import FromSpecifiedMarkersRunner

    config = TAXONOMY_CONFIG[species]
    output_dir = Path(output_dir)
    csv_result_path = output_dir / "mapping_output.csv"
    # extended_result_path is required by the mapper but we don't need it
    json_tmp = tempfile.NamedTemporaryFile(suffix=".json", delete=False)
    json_tmp.close()

    # Use SLURM allocation if available, otherwise detect physical RAM
    slurm_mem = os.environ.get("SLURM_MEM_PER_NODE")
    if slurm_mem:
        total_ram_gb = int(slurm_mem) // 1024  # SLURM reports in MB
    else:
        total_ram_gb = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES") // (1024 ** 3)
    max_gb = max(1, total_ram_gb // 2)
    logger.info("System: %d cores, %dGB RAM → using n_processors=%d, max_gb=%d",
                multiprocessing.cpu_count(), total_ram_gb, n_processors, max_gb)

    if flatten:
        bootstrap_iteration = 1
        bootstrap_factor = 1.0
        logger.info("Flatten mode: skipping hierarchy, bootstrap_iteration=1")
    else:
        bootstrap_iteration = 100
        bootstrap_factor = 0.9

    try:
        config_dict = {
            "query_path": str(validated_h5ad_path),
            "extended_result_path": json_tmp.name,
            "csv_result_path": str(csv_result_path),
            "drop_level": config["drop_level"],
            "flatten": flatten,
            "tmp_dir": os.environ.get("TMPDIR", "/tmp"),
            "max_gb": max_gb,
            "precomputed_stats": {"path": str(reference_dir / config["precomputed_stats"])},
            "query_markers": {"serialized_lookup": str(reference_dir / config["markers"])},
            "type_assignment": {
                "normalization": config["normalization"],
                "n_processors": n_processors,
                "chunk_size": 2000,
                "bootstrap_iteration": bootstrap_iteration,
                "bootstrap_factor": bootstrap_factor,
            },
        }
        runner = FromSpecifiedMarkersRunner(args=[], input_data=config_dict)
        runner.run()
    finally:
        os.unlink(json_tmp.name)
    return csv_result_path


def _read_scalar_json(dset):
    """Allen precomputed_stats h5 stores col_names / cluster_to_row / taxonomy_tree
    as JSON-encoded scalar datasets. Return the parsed Python object."""
    raw = dset[()]
    if isinstance(raw, bytes):
        raw = raw.decode()
    return json.loads(raw)


def _load_reference_cluster_means(precomputed_stats_path):
    """Read cluster-mean expression + taxonomy from Allen precomputed_stats h5.

    Returns (ref_gene_ids, cluster_ids, cluster_means, taxonomy_tree).
    cluster_means is (n_clusters, n_genes) of raw counts averaged per cluster.
    cluster_ids are the raw taxonomy IDs (e.g. 'CS20230722_CLUS_0001').
    """
    import h5py

    with h5py.File(precomputed_stats_path, "r") as f:
        gene_key = "col_names" if "col_names" in f else "gene_names"
        ref_genes = _read_scalar_json(f[gene_key])

        mapping = _read_scalar_json(f["cluster_to_row"])
        row_to_cluster = {int(v): k for k, v in mapping.items()}
        cluster_ids = [row_to_cluster[i] for i in range(len(row_to_cluster))]

        sums = f["sum"][:].astype(np.float64)
        n_cells = f["n_cells"][:].astype(np.float64)
        cluster_means = sums / np.maximum(n_cells[:, None], 1.0)

        taxonomy_tree = None
        if "taxonomy_tree" in f:
            try:
                taxonomy_tree = _read_scalar_json(f["taxonomy_tree"])
            except Exception:
                taxonomy_tree = None

    return ref_genes, cluster_ids, cluster_means, taxonomy_tree


def _resolve_name(name_mapper, level_code, node_id):
    """Translate an Allen taxonomy ID to a readable name via name_mapper;
    fall back to the raw ID if no mapping is available."""
    if not name_mapper or not level_code:
        return node_id
    entry = (name_mapper.get(level_code) or {}).get(node_id)
    if isinstance(entry, dict):
        return entry.get("name", node_id)
    if isinstance(entry, str):
        return entry
    return node_id


def _leaf_to_parent_maps(taxonomy_tree):
    """From Allen taxonomy tree, build {friendly_level_name: {leaf_id: readable_label}}.

    Uses 'hierarchy_mapper' to rename level codes like CCN20230722_CLAS → 'class',
    and 'name_mapper' to turn parent IDs into readable names.
    """
    if not isinstance(taxonomy_tree, dict):
        return {}
    hierarchy = taxonomy_tree.get("hierarchy")
    if not hierarchy or len(hierarchy) < 2:
        return {}

    hierarchy_mapper = taxonomy_tree.get("hierarchy_mapper", {}) or {}
    name_mapper = taxonomy_tree.get("name_mapper", {}) or {}
    level_children = {lvl: taxonomy_tree.get(lvl, {}) for lvl in hierarchy[:-1]}

    leaf_to_ancestors = {}  # leaf_id -> {level_code: parent_id}

    def recurse(level_idx, node, ancestors):
        if level_idx == len(hierarchy) - 1:
            leaf_to_ancestors[node] = dict(ancestors)
            return
        current_level = hierarchy[level_idx]
        next_ancestors = dict(ancestors)
        next_ancestors[current_level] = node
        for child in level_children.get(current_level, {}).get(node, []):
            recurse(level_idx + 1, child, next_ancestors)

    top_level = hierarchy[0]
    for top_node in level_children.get(top_level, {}).keys():
        recurse(0, top_node, {})

    per_level = {}
    for level_code in hierarchy[:-1]:
        friendly = hierarchy_mapper.get(level_code, level_code.lower())
        per_level[friendly] = {
            leaf: _resolve_name(name_mapper, level_code, anc.get(level_code, ""))
            for leaf, anc in leaf_to_ancestors.items()
        }
    return per_level


def run_mapping_correlation(validated_h5ad_path, output_dir, reference_dir, species):
    """Assign each cell to the reference cluster with the highest Pearson
    correlation of log1p-normalized expression profiles.

    Output CSV schema matches the hierarchical runner's downstream requirements:
    cell_id, cluster_name, cluster_label, plus parent taxonomy levels we can
    derive (e.g. class_name), plus correlation_coefficient in place of
    bootstrapping_probability.
    """
    import scipy.sparse as sp

    config = TAXONOMY_CONFIG[species]
    output_dir = Path(output_dir)
    csv_result_path = output_dir / "mapping_output.csv"
    precomputed_path = reference_dir / config["precomputed_stats"]

    logger.info("Correlation mapping: loading reference stats from %s", precomputed_path)
    ref_genes, cluster_ids, cluster_means, taxonomy_tree = (
        _load_reference_cluster_means(precomputed_path)
    )
    logger.info("Reference: %d clusters × %d genes", len(cluster_ids), len(ref_genes))

    # Resolve readable cluster names via name_mapper (fallback to ID).
    hierarchy = (taxonomy_tree or {}).get("hierarchy", [])
    leaf_level_code = hierarchy[-1] if hierarchy else None
    name_mapper = (taxonomy_tree or {}).get("name_mapper", {}) or {}
    cluster_names = [
        _resolve_name(name_mapper, leaf_level_code, cid) for cid in cluster_ids
    ]

    adata = anndata.read_h5ad(validated_h5ad_path)
    query_genes = list(adata.var.index)
    X = adata.X
    if sp.issparse(X):
        X = X.toarray()
    X = X.astype(np.float64, copy=False)

    ref_gene_to_idx = {g: i for i, g in enumerate(ref_genes)}
    shared = [(qi, ref_gene_to_idx[g]) for qi, g in enumerate(query_genes)
              if g in ref_gene_to_idx]
    if len(shared) < 10:
        raise ValueError(
            f"Only {len(shared)} genes overlap between query and reference. "
            "Correlation mapping requires enough shared genes to be meaningful."
        )
    logger.info("Correlation mapping: %d/%d query genes align with reference.",
                len(shared), len(query_genes))
    q_idx = np.array([s[0] for s in shared])
    r_idx = np.array([s[1] for s in shared])

    def _lognorm(M):
        totals = np.maximum(M.sum(axis=1, keepdims=True), 1.0)
        return np.log1p(M / totals * 1e4)

    def _zscore_rows(M):
        mu = M.mean(axis=1, keepdims=True)
        sd = M.std(axis=1, keepdims=True)
        sd = np.where(sd == 0, 1.0, sd)
        return (M - mu) / sd

    Q = _lognorm(X[:, q_idx])
    R = _lognorm(cluster_means[:, r_idx])
    Qz = _zscore_rows(Q)
    Rz = _zscore_rows(R)

    logger.info("Computing %d × %d correlation matrix...", Qz.shape[0], Rz.shape[0])
    corr = (Qz @ Rz.T) / Qz.shape[1]
    best = corr.argmax(axis=1)
    best_r = corr[np.arange(len(best)), best]
    assigned_ids = [cluster_ids[i] for i in best]
    assigned_names = [cluster_names[i] for i in best]

    result = pd.DataFrame({
        "cell_id": adata.obs.index.astype(str),
        "cluster_label": assigned_ids,
        "cluster_name": assigned_names,
        "correlation_coefficient": best_r,
    })

    parent_maps = _leaf_to_parent_maps(taxonomy_tree)
    for friendly_level, leaf_to_label in parent_maps.items():
        result[f"{friendly_level}_name"] = [
            leaf_to_label.get(cid, "") for cid in assigned_ids
        ]

    with open(csv_result_path, "w") as fh:
        fh.write("# mapping method: correlation (Pearson, log1p CPT-10k)\n")
        fh.write(f"# reference: {precomputed_path.name}\n")
        fh.write(f"# shared_genes: {len(shared)}\n")
        result.to_csv(fh, index=False)

    logger.info("Correlation mapping complete → %s", csv_result_path)
    return csv_result_path


def join_results(mapping_csv_path, metadata_df):
    mapping_df = pd.read_csv(mapping_csv_path, comment='#')
    mapping_df["cell_id"] = mapping_df["cell_id"].astype(str)

    meta_id_col = next(c for c in ("EntityID", "id", "cell_id") if c in metadata_df.columns)

    # Build the same compound key used in build_h5ad
    metadata_df = metadata_df.copy()
    if "_sample_id" in metadata_df.columns:
        metadata_df["_join_key"] = metadata_df["_sample_id"].astype(str) + "_" + metadata_df[meta_id_col].astype(str)
    else:
        metadata_df["_join_key"] = metadata_df[meta_id_col].astype(str)

    merged = metadata_df.merge(mapping_df, left_on="_join_key", right_on="cell_id", how="left")
    merged = merged.drop(columns=["_join_key"])
    return merged

def plot_cell_types(merged_df, output_dir):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    df = merged_df

    # Mouse has class_name, human has supercluster_name. Correlation mode may
    # only populate cluster_name. Fall back through them in preference order.
    class_col = None
    for candidate in ("class_name", "supercluster_name", "subclass_name",
                      "cluster_name", "cluster_label"):
        if candidate in df.columns and not df[candidate].isna().all():
            class_col = candidate
            break
    if class_col is None:
        logger.warning("No cell type column found, skipping plots.")
        return
    logger.info("Plotting cell types using column '%s'", class_col)

    # Pick top 3 most common cell type classes to highlight
    top3 = df[class_col].value_counts().head(3).index.tolist()

    for cell_type in top3:
        plt.style.use("dark_background")
        fig, ax = plt.subplots(figsize=(14, 14))

        # Plot all cells in grey first
        ax.scatter(df["center_x"], df["center_y"],
                   c="dimgrey", s=0.3, alpha=0.3, rasterized=True)

        # Overlay the highlighted cell type in color
        mask = df[class_col] == cell_type
        ax.scatter(df.loc[mask, "center_x"], df.loc[mask, "center_y"],
                   c="crimson", s=0.8, alpha=0.8, rasterized=True,
                   label=f"{cell_type} (n={mask.sum():,})")

        ax.set_title(f"{cell_type} — {mask.sum():,} of {len(df):,} cells")
        ax.set_aspect("equal")
        ax.legend(loc='upper right', fontsize=8, framealpha=0.3)

        plt.tight_layout()
        safe_name = str(cell_type).replace(" ", "_").replace("/", "-")
        out_path = Path(output_dir) / f"check_celltype_{safe_name}.png"
        plt.savefig(out_path, dpi=150)
        plt.close()
        logger.info("Saved %s", out_path)


def plot_correlation_diagnostic(merged_df, output_dir, bins=50):
    """Histogram of per-cell best-match correlation coefficients.

    Shows how confident the correlation mapper was across the dataset — a
    left-skewed distribution concentrated near 1 means most cells had a clear
    match; a flat or near-zero distribution means many cells are weakly
    differentiated against the reference (often because of sparse expression
    on the shared gene panel).
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    if "correlation_coefficient" not in merged_df.columns:
        logger.warning("No correlation_coefficient column; skipping diagnostic.")
        return

    r = pd.to_numeric(merged_df["correlation_coefficient"], errors="coerce").dropna()
    if r.empty:
        logger.warning("No numeric correlation values; skipping diagnostic.")
        return

    median = float(r.median())
    mean = float(r.mean())
    n_zero = int((r <= 0).sum())

    plt.style.use("default")
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(r, bins=bins, color="steelblue", edgecolor="black", linewidth=0.3)
    ax.axvline(median, color="crimson", linestyle="--", linewidth=1.2,
               label=f"median = {median:.3f}")
    ax.axvline(mean, color="darkorange", linestyle=":", linewidth=1.2,
               label=f"mean = {mean:.3f}")
    ax.set_xlabel("Pearson r (best-match reference cluster)")
    ax.set_ylabel("number of cells")
    ax.set_title(f"Correlation-method confidence — {len(r):,} cells "
                 f"({n_zero:,} with r ≤ 0)")
    ax.legend(loc="upper left", fontsize=9, framealpha=0.9)
    ax.grid(axis="y", linestyle=":", alpha=0.4)
    plt.tight_layout()
    out_path = Path(output_dir) / "correlation_diagnostic.png"
    plt.savefig(out_path, dpi=150)
    plt.close()
    logger.info("Saved %s", out_path)

def print_summary(merged_df, query_h5ad_path, original_cbg_gene_count):
    validated_adata = anndata.read_h5ad(query_h5ad_path, backed="r")
    n_genes = len(validated_adata.var.index)
    gene_overlap_pct = n_genes / original_cbg_gene_count * 100

    print(f"Gene overlap: {n_genes}/{original_cbg_gene_count} ({gene_overlap_pct:.1f}%)")
    print(f"Total cells mapped: {len(merged_df)}")

    cluster_col = "cluster_label" if "cluster_label" in merged_df.columns else merged_df.columns[-1]
    print(f"\nTop 5 assigned cluster labels ({cluster_col}):")
    print(merged_df[cluster_col].value_counts().head(5).to_string())


if __name__ == "__main__":
    args = parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        stream=sys.stdout,
    )

    input_dir = Path(args.input_dir).expanduser()

    reference_dir = resolve_reference_dir(args.reference_dir)
    output_dir = Path(args.output_dir).expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)

    if input_dir.is_file() and input_dir.suffix == ".h5ad":
        metadata_df, cbg_df = load_h5ad(input_dir)
    else:
        metadata_path = input_dir / "cell_metadata.csv"
        cbg_path = input_dir / "cell_by_gene.csv"

        if not metadata_path.exists():
            raise FileNotFoundError(f"cell_metadata.csv not found in {input_dir}")
        if not cbg_path.exists():
            raise FileNotFoundError(f"cell_by_gene.csv not found in {input_dir}")

        metadata_df, cbg_df = load_inputs(metadata_path, cbg_path)
    original_gene_count = len(cbg_df.columns) - 1

    gene_mapping_path = reference_dir / TAXONOMY_CONFIG[args.species]["gene_mapping"]
    cbg_df = translate_genes_to_ensembl(cbg_df, gene_mapping_path)

    query_h5ad_path = build_h5ad(cbg_df, metadata_df, output_dir)
    mapping_csv_path = run_mapping(
        query_h5ad_path, output_dir, reference_dir, args.species, args.n_processors,
        flatten=args.flatten, method=args.method,
    )
    merged_df = join_results(mapping_csv_path, metadata_df)
    plot_cell_types(merged_df, output_dir)
    if args.method == "correlation":
        plot_correlation_diagnostic(merged_df, output_dir)
    print_summary(merged_df, query_h5ad_path, original_gene_count)

    # Correlation mode: keep only mapping_output.csv + PNGs.
    if args.method == "correlation":
        try:
            Path(query_h5ad_path).unlink()
            logger.info("Removed intermediate %s", query_h5ad_path)
        except FileNotFoundError:
            pass
