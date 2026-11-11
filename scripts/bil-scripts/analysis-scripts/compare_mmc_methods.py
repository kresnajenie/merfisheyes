#!/usr/bin/env python3
"""
compare_mmc_methods.py

Compare hierarchical (mmc_output/) vs correlation (mmc_output_corr/) MapMyCells
results across one or many datasets.

Inputs (per sample, on the BIL cluster by default):
  ${MEYES_BASE}/<sample>/mmc_output/mapping_output.csv       (hierarchical)
  ${MEYES_BASE}/<sample>/mmc_output_corr/mapping_output.csv  (correlation)

Outputs (in --out-dir):
  comparison.tsv               one row per dataset (default mode)
  comparison.png               summary plot for the dataset-level TSV
  comparison_per_slice.tsv     one row per slice, with masked/unmasked metrics
                                (when --per-slice is set; see below)
  comparison_per_slice.png     summary plot for the per-slice TSV
  per_sample/<sample>.png      per-cell histograms + class-level confusion
                                matrix (only when --per-sample-plots is set,
                                or when run on a single sample in default mode)

Per-slice mode (--per-slice):
  - Slices come from combined_output/cell_metadata.csv `_sample_id` column.
  - Mapping CSV cell_ids are joined as "{_sample_id}_{cell_id}" (built by
    map_my_cell.build_h5ad).
  - Artifact mask is read from combined_output/<--mask-name>{,.csv}; cells
    flagged as artifacts are subtracted to produce the *_masked metrics.

Usage:
  # all rows in samples.csv (per dataset)
  python compare_mmc_methods.py samples.csv

  # all rows in samples.csv (per slice, with masked metrics)
  python compare_mmc_methods.py samples.csv --per-slice

  # one sample (also writes per_sample/<sample>.png in default mode)
  python compare_mmc_methods.py --sample ace-dip-use

  # custom roots
  python compare_mmc_methods.py samples.csv \
      --meyes-base /bil/data/meyes \
      --out-dir   /bil/data/meyes/_mmc_comparison
"""

import argparse
import csv
import logging
import os
import re
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

logger = logging.getLogger("compare_mmc")

DEFAULT_MEYES_BASE = "/bil/data/meyes"
HIER_DIR = "mmc_output"
CORR_DIR = "mmc_output_corr"
MAPPING_FILE = "mapping_output.csv"
COMBINED_DIR = "combined_output"
CELL_METADATA_FILE = "cell_metadata.csv"
MASK_NAME_RE = re.compile(r"^artifact_mask_p(\d+)(?:\.csv)?$", re.IGNORECASE)

LEVEL_COLS = ["class_name", "subclass_name", "supertype_name",
              "cluster_name", "cluster_label"]

# Allen's FromSpecifiedMarkersRunner writes one bootstrapping_probability per
# taxonomy level, e.g. cluster_bootstrapping_probability. The leaf level is
# the natural match for correlation_coefficient (also leaf-level).
HIER_PROB_LEAF_PRIORITY = [
    "cluster_bootstrapping_probability",
    "supertype_bootstrapping_probability",
    "subclass_bootstrapping_probability",
    "class_bootstrapping_probability",
    "bootstrapping_probability",  # fallback if upstream changed
]


def _hier_prob_column(df: pd.DataFrame) -> str | None:
    """Pick the deepest available bootstrapping_probability column."""
    for col in HIER_PROB_LEAF_PRIORITY:
        if col in df.columns:
            return col
    return None


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("samples_csv", nargs="?",
                   help="CSV with rows of `sample_name,input_path` (input_path "
                        "is ignored here). If omitted, --sample is required.")
    p.add_argument("--sample", action="append", default=[],
                   help="Process only this sample. Repeat to pass multiple.")
    p.add_argument("--meyes-base", default=DEFAULT_MEYES_BASE,
                   help=f"Root containing per-sample dirs (default: {DEFAULT_MEYES_BASE}).")
    p.add_argument("--out-dir", default=None,
                   help="Where to write comparison.tsv and plots "
                        "(default: <meyes-base>/_mmc_comparison).")
    p.add_argument("--species-map", default=None,
                   help="Optional CSV with columns sample_name,species. If "
                        "omitted, species is left blank (or read from a "
                        "third column of samples_csv if present).")
    p.add_argument("--per-sample-plots", action="store_true",
                   help="Also emit per_sample/<sample>.png for every sample.")
    p.add_argument("--low-conf-threshold", type=float, default=0.5,
                   help="Threshold for hier_low_conf_frac (default 0.5).")
    p.add_argument("--per-slice", action="store_true",
                   help="Emit one row per slice (uses _sample_id from "
                        "combined_output/cell_metadata.csv) instead of one "
                        "row per dataset. Writes comparison_per_slice.tsv "
                        "and applies the artifact mask to add *_masked "
                        "metric columns alongside the unmasked ones.")
    p.add_argument("--mask-name", default=None,
                   help="Override the artifact mask filename inside "
                        "combined_output/ (e.g. artifact_mask_p50). "
                        "Default: auto-pick the highest percentile found "
                        "(artifact_mask_p<N> with the largest N).")
    return p.parse_args()


def read_samples(args):
    """Returns list of (sample_name, species_or_None, input_path_or_None)."""
    rows = []
    species_lookup = {}
    path_lookup: dict[str, str] = {}

    if args.species_map:
        sm = pd.read_csv(args.species_map)
        if "sample_name" not in sm.columns or "species" not in sm.columns:
            raise ValueError("--species-map needs columns sample_name,species")
        species_lookup = dict(zip(sm["sample_name"].astype(str),
                                  sm["species"].astype(str)))

    # If samples_csv exists, read it for input_path lookup even when --sample
    # filters down to a subset (so we can still resolve publication dates).
    if args.samples_csv and Path(args.samples_csv).exists():
        with open(args.samples_csv) as fh:
            for r in csv.reader(fh):
                if not r:
                    continue
                n = r[0].strip()
                if not n or n.startswith("#"):
                    continue
                if len(r) >= 2 and r[1].strip():
                    path_lookup[n] = r[1].strip()

    if args.sample:
        for s in args.sample:
            rows.append((s, species_lookup.get(s), path_lookup.get(s)))
        return rows

    if not args.samples_csv:
        raise SystemExit("Must provide samples_csv or --sample")

    with open(args.samples_csv) as fh:
        reader = csv.reader(fh)
        for r in reader:
            if not r:
                continue
            name = r[0].strip()
            if not name or name.startswith("#"):
                continue
            sp = None
            if len(r) >= 3 and r[2].strip():
                sp = r[2].strip()
            sp = species_lookup.get(name, sp)
            input_path = r[1].strip() if len(r) >= 2 and r[1].strip() else None
            rows.append((name, sp, input_path))
    return rows


def path_publication_date(path_str: str | None) -> str:
    """Return YYYY-MM-DD mtime of the BIL data path (the date `ls -l` shows),
    or '' if the path is missing or unreadable."""
    if not path_str:
        return ""
    p = Path(path_str.rstrip("/"))
    try:
        st = os.stat(p)
    except (OSError, FileNotFoundError) as e:
        logger.warning("stat failed for %s: %s", p, e)
        return ""
    return datetime.fromtimestamp(st.st_mtime).strftime("%Y-%m-%d")


def load_mapping_csv(path: Path) -> pd.DataFrame | None:
    if not path.exists():
        return None
    try:
        return pd.read_csv(path, comment="#", low_memory=False)
    except Exception as e:
        logger.warning("Failed to read %s: %s", path, e)
        return None


def agreement(hier: pd.DataFrame, corr: pd.DataFrame, col: str) -> float | None:
    if col not in hier.columns or col not in corr.columns:
        return None
    merged = hier[["cell_id", col]].merge(
        corr[["cell_id", col]], on="cell_id", suffixes=("_h", "_c"))
    if merged.empty:
        return None
    return float((merged[f"{col}_h"].astype(str) ==
                  merged[f"{col}_c"].astype(str)).mean())


def _meta_id_col(df: pd.DataFrame) -> str | None:
    for c in ("EntityID", "id", "cell_id"):
        if c in df.columns:
            return c
    return None


def load_cell_metadata(metadata_path: Path) -> pd.DataFrame | None:
    if not metadata_path.exists():
        return None
    try:
        return pd.read_csv(metadata_path, low_memory=False)
    except Exception as e:
        logger.warning("Failed to read %s: %s", metadata_path, e)
        return None


def pick_mask_name(combined_dir: Path, override: str | None) -> str | None:
    """Pick the artifact mask filename to use for one dataset.

    If --mask-name was passed, return it. Otherwise scan combined_output/
    for entries matching artifact_mask_p<N>(.csv)? and return the name with
    the highest N. Returns None when nothing matches."""
    if override:
        return override
    if not combined_dir.is_dir():
        return None
    best: tuple[int, str] | None = None
    for entry in combined_dir.iterdir():
        m = MASK_NAME_RE.match(entry.name)
        if not m:
            continue
        n = int(m.group(1))
        if best is None or n > best[0]:
            best = (n, entry.name)
    return best[1] if best else None


def load_artifact_keys(mask_dir: Path, mask_name: str,
                       cell_metadata: pd.DataFrame | None) -> set | None:
    """Return the set of compound keys ({_sample_id}_{cell_id}) that are
    artifacts according to the mask CSV. Returns None when the mask cannot
    be located (so the caller can distinguish 'no mask' from 'empty mask').

    The mask file is looked up at:
      <mask_dir>/<mask_name>            (path as given)
      <mask_dir>/<mask_name>.csv        (with .csv extension)
      <mask_dir>/<mask_name>/*.csv      (if name is a directory)
    """
    candidates = [mask_dir / mask_name, mask_dir / f"{mask_name}.csv"]
    if (mask_dir / mask_name).is_dir():
        candidates.extend(sorted((mask_dir / mask_name).glob("*.csv")))
    mask_path = next((c for c in candidates if c.is_file()), None)
    if mask_path is None:
        return None

    try:
        mask_df = pd.read_csv(mask_path, low_memory=False)
    except Exception as e:
        logger.warning("Failed to read mask %s: %s", mask_path, e)
        return None

    id_col = _meta_id_col(mask_df)
    sample_col = "_sample_id" if "_sample_id" in mask_df.columns else None
    artifact_col = next(
        (c for c in mask_df.columns
         if c.lower() in ("is_artifact", "artifact", "artifact_mask",
                          "mask", "p25", "is_mask")
         or c.lower().startswith("artifact_mask")),
        None,
    )

    # Path A: mask has explicit ID columns.
    if id_col and artifact_col:
        flagged = mask_df[mask_df[artifact_col].astype(bool)]
        if sample_col:
            return {
                f"{s}_{c}" for s, c in zip(
                    flagged[sample_col].astype(str),
                    flagged[id_col].astype(str),
                )
            }
        if cell_metadata is not None and "_sample_id" in cell_metadata.columns:
            meta_id = _meta_id_col(cell_metadata)
            if meta_id is None:
                return set()
            lookup = dict(zip(
                cell_metadata[meta_id].astype(str),
                cell_metadata["_sample_id"].astype(str),
            ))
            keys = set()
            for cid in flagged[id_col].astype(str):
                sid = lookup.get(cid)
                if sid is not None:
                    keys.add(f"{sid}_{cid}")
            return keys
        return {str(c) for c in flagged[id_col].astype(str)}

    # Path B: mask is row-aligned to cell_metadata (no IDs in mask).
    if cell_metadata is not None and len(mask_df) == len(cell_metadata):
        meta_id = _meta_id_col(cell_metadata)
        if meta_id is None:
            return set()
        # Pick the first numeric/boolean-looking column.
        col = next(
            (c for c in mask_df.columns
             if pd.api.types.is_bool_dtype(mask_df[c])
             or pd.api.types.is_numeric_dtype(mask_df[c])),
            mask_df.columns[0],
        )
        flag = mask_df[col].astype(bool).to_numpy()
        sids = (cell_metadata["_sample_id"].astype(str).to_numpy()
                if "_sample_id" in cell_metadata.columns
                else np.array([""] * len(cell_metadata)))
        cids = cell_metadata[meta_id].astype(str).to_numpy()
        return {f"{s}_{c}" for s, c, f in zip(sids, cids, flag) if f}

    logger.warning("Mask %s has unrecognized schema; skipping.", mask_path)
    return None


def _filter_mapping_to_keys(df: pd.DataFrame | None, keys: set) -> pd.DataFrame | None:
    if df is None or "cell_id" not in df.columns:
        return df
    return df[df["cell_id"].astype(str).isin(keys)]


def compute_metrics(hier: pd.DataFrame | None, corr: pd.DataFrame | None,
                    low_conf_threshold: float, suffix: str = "") -> dict:
    """All summary metrics for a hier/corr pair, with column names suffixed."""
    out: dict = {}
    if hier is not None and len(hier):
        prob_col = _hier_prob_column(hier)
        if prob_col:
            p = pd.to_numeric(hier[prob_col], errors="coerce").dropna()
            if not p.empty:
                out[f"hier_mean_prob{suffix}"] = float(p.mean())
                out[f"hier_median_prob{suffix}"] = float(p.median())
                out[f"hier_low_conf_frac{suffix}"] = float((p < low_conf_threshold).mean())
    if corr is not None and len(corr) and "correlation_coefficient" in corr.columns:
        r = pd.to_numeric(corr["correlation_coefficient"], errors="coerce").dropna()
        if not r.empty:
            out[f"corr_mean_r{suffix}"] = float(r.mean())
            out[f"corr_median_r{suffix}"] = float(r.median())
            out[f"corr_neg_frac{suffix}"] = float((r <= 0).mean())
    if hier is not None and corr is not None and len(hier) and len(corr):
        for level_col in LEVEL_COLS:
            agr = agreement(hier, corr, level_col)
            if agr is not None:
                friendly = level_col.replace("_name", "").replace("_label", "_label")
                out[f"agreement_{friendly}{suffix}"] = agr
    return out


def summarize_dataset_per_slice(sample: str, species: str | None,
                                 meyes_base: Path, low_conf_threshold: float,
                                 input_path: str | None,
                                 mask_override: str | None) -> list[dict]:
    """Return one dict per slice for this sample. If cell_metadata.csv is
    missing or has no _sample_id column, returns a single row with
    slice_id='' (degenerate single-slice case)."""
    sample_dir = meyes_base / sample
    metadata_path = sample_dir / COMBINED_DIR / CELL_METADATA_FILE
    hier_full = load_mapping_csv(sample_dir / HIER_DIR / MAPPING_FILE)
    corr_full = load_mapping_csv(sample_dir / CORR_DIR / MAPPING_FILE)
    cell_meta = load_cell_metadata(metadata_path)

    mask_name = pick_mask_name(sample_dir / COMBINED_DIR, mask_override)
    if mask_name is None:
        artifact_keys = None
        logger.info("[%s] no artifact_mask_p* file in combined_output/ — "
                    "masked metrics will be blank", sample)
    else:
        artifact_keys = load_artifact_keys(sample_dir / COMBINED_DIR,
                                            mask_name, cell_meta)
        if artifact_keys is None:
            logger.info("[%s] mask %s could not be parsed — masked metrics "
                        "will be blank", sample, mask_name)
        else:
            logger.info("[%s] using mask %s (%d artifact cells flagged)",
                        sample, mask_name, len(artifact_keys))

    pub_date = path_publication_date(input_path)
    base = {
        "sample": sample,
        "species": species or "",
        "publication_date": pub_date,
        "input_path": input_path or "",
        "mask_name": mask_name or "",
    }

    rows: list[dict] = []

    if cell_meta is None:
        logger.warning("[%s] no %s; emitting whole-dataset row with "
                       "slice_id=''", sample, metadata_path)
        rows.append(_build_slice_row(
            base, slice_id="", slice_keys=None,
            hier_full=hier_full, corr_full=corr_full,
            artifact_keys=artifact_keys,
            low_conf_threshold=low_conf_threshold,
        ))
        return rows

    meta_id = _meta_id_col(cell_meta)
    if meta_id is None:
        logger.warning("[%s] cell_metadata.csv has no recognized id column; "
                       "skipping", sample)
        return []

    if "_sample_id" not in cell_meta.columns:
        logger.info("[%s] no _sample_id column; treating as a single slice",
                    sample)
        slice_iter: list[tuple[str, pd.DataFrame]] = [("", cell_meta)]
    else:
        cell_meta["_sample_id"] = cell_meta["_sample_id"].astype(str)
        slice_iter = [(sid, group) for sid, group
                      in cell_meta.groupby("_sample_id", sort=True)]

    for sid, group in slice_iter:
        slice_keys = {
            f"{sid}_{cid}" for cid in group[meta_id].astype(str)
        }
        rows.append(_build_slice_row(
            base, slice_id=sid, slice_keys=slice_keys,
            hier_full=hier_full, corr_full=corr_full,
            artifact_keys=artifact_keys,
            low_conf_threshold=low_conf_threshold,
        ))
    return rows


def _build_slice_row(base: dict, slice_id: str, slice_keys: set | None,
                     hier_full: pd.DataFrame | None,
                     corr_full: pd.DataFrame | None,
                     artifact_keys: set | None,
                     low_conf_threshold: float) -> dict:
    row = dict(base)
    row["slice_id"] = slice_id

    if slice_keys is None:
        # Degenerate: no cell_metadata to define a slice — use the entire
        # mapping CSVs and treat the artifact set as global.
        hier_slice = hier_full
        corr_slice = corr_full
        n_total = max(
            len(hier_full) if hier_full is not None else 0,
            len(corr_full) if corr_full is not None else 0,
        )
        in_mask_keys: set | None = None
    else:
        hier_slice = _filter_mapping_to_keys(hier_full, slice_keys)
        corr_slice = _filter_mapping_to_keys(corr_full, slice_keys)
        n_total = len(slice_keys)
        in_mask_keys = (slice_keys - artifact_keys) if artifact_keys is not None else None

    if artifact_keys is not None and slice_keys is not None:
        n_artifact = len(slice_keys & artifact_keys)
    else:
        n_artifact = np.nan

    n_in_mask = (n_total - n_artifact) if isinstance(n_artifact, int) else np.nan
    frac_in_mask = (n_in_mask / n_total) if isinstance(n_in_mask, int) and n_total else np.nan

    row["n_cells_total"] = n_total
    row["n_cells_artifact"] = n_artifact
    row["n_cells_in_mask"] = n_in_mask
    row["frac_in_mask"] = frac_in_mask
    row["n_cells_hier"] = (int(len(hier_slice)) if hier_slice is not None
                            else np.nan)
    row["n_cells_corr"] = (int(len(corr_slice)) if corr_slice is not None
                            else np.nan)
    row["hier_prob_col"] = (_hier_prob_column(hier_slice) or ""
                             if hier_slice is not None and len(hier_slice) else "")

    row.update(compute_metrics(hier_slice, corr_slice, low_conf_threshold))

    if in_mask_keys is not None:
        hier_masked = _filter_mapping_to_keys(hier_slice, in_mask_keys)
        corr_masked = _filter_mapping_to_keys(corr_slice, in_mask_keys)
        row["n_cells_hier_masked"] = (int(len(hier_masked))
                                       if hier_masked is not None else np.nan)
        row["n_cells_corr_masked"] = (int(len(corr_masked))
                                       if corr_masked is not None else np.nan)
        row.update(compute_metrics(hier_masked, corr_masked,
                                   low_conf_threshold, suffix="_masked"))

    return row


def summarize_sample(sample: str, species: str | None, meyes_base: Path,
                     low_conf_threshold: float,
                     input_path: str | None = None) -> dict:
    sample_dir = meyes_base / sample
    hier_path = sample_dir / HIER_DIR / MAPPING_FILE
    corr_path = sample_dir / CORR_DIR / MAPPING_FILE

    hier = load_mapping_csv(hier_path)
    corr = load_mapping_csv(corr_path)

    row: dict = {
        "sample": sample,
        "species": species or "",
        "publication_date": path_publication_date(input_path),
        "input_path": input_path or "",
    }

    # Pre-seed numeric columns so the DataFrame has them even when no sample
    # populates them (avoids KeyError in plotting).
    for k in ("n_cells", "n_cells_hier", "n_cells_corr",
              "hier_mean_prob", "hier_median_prob", "hier_low_conf_frac",
              "corr_mean_r", "corr_median_r", "corr_neg_frac"):
        row[k] = np.nan

    if hier is not None:
        row["hier_path"] = str(hier_path)
        row["n_cells_hier"] = len(hier)
        prob_col = _hier_prob_column(hier)
        row["hier_prob_col"] = prob_col or ""
        if prob_col:
            p = pd.to_numeric(hier[prob_col], errors="coerce").dropna()
            if not p.empty:
                row["hier_mean_prob"] = float(p.mean())
                row["hier_median_prob"] = float(p.median())
                row["hier_low_conf_frac"] = float((p < low_conf_threshold).mean())
    else:
        row["hier_path"] = ""

    if corr is not None:
        row["corr_path"] = str(corr_path)
        row["n_cells_corr"] = len(corr)
        if "correlation_coefficient" in corr.columns:
            r = pd.to_numeric(corr["correlation_coefficient"], errors="coerce").dropna()
            if not r.empty:
                row["corr_mean_r"] = float(r.mean())
                row["corr_median_r"] = float(r.median())
                row["corr_neg_frac"] = float((r <= 0).mean())
    else:
        row["corr_path"] = ""

    if hier is not None and corr is not None:
        row["n_cells"] = int(max(len(hier), len(corr)))
        for level_col in LEVEL_COLS:
            agr = agreement(hier, corr, level_col)
            if agr is not None:
                friendly = level_col.replace("_name", "").replace("_label", "_label")
                row[f"agreement_{friendly}"] = agr

        h_genes = _maybe_read_gene_count(hier_path)
        c_genes = _maybe_read_gene_count(corr_path)
        if h_genes is not None:
            row["hier_shared_genes"] = h_genes
        if c_genes is not None:
            row["corr_shared_genes"] = c_genes
    elif hier is not None:
        row["n_cells"] = len(hier)
    elif corr is not None:
        row["n_cells"] = len(corr)

    return row


def _maybe_read_gene_count(csv_path: Path) -> int | None:
    """Pull `# shared_genes: N` out of header comments if the writer added it."""
    try:
        with open(csv_path) as fh:
            for line in fh:
                if not line.startswith("#"):
                    return None
                if "shared_genes" in line:
                    parts = line.strip().split(":")
                    if len(parts) == 2:
                        try:
                            return int(parts[1].strip())
                        except ValueError:
                            return None
    except OSError:
        pass
    return None


def make_summary_plot(df: pd.DataFrame, out_path: Path, low_conf_threshold: float):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # 1. scatter: mean prob vs mean r
    ax = axes[0, 0]
    needed = [c for c in ("hier_mean_prob", "corr_mean_r") if c in df.columns]
    sub = df.dropna(subset=needed) if len(needed) == 2 else df.iloc[0:0]
    if not sub.empty:
        species_vals = sub["species"].fillna("").replace("", "unknown")
        size_col = ("n_cells_total" if "n_cells_total" in sub.columns
                    else "n_cells")
        for sp, group in sub.groupby(species_vals):
            sizes = np.clip(group[size_col].fillna(1).astype(float) / 200,
                            10, 400)
            ax.scatter(group["hier_mean_prob"], group["corr_mean_r"],
                       s=sizes, alpha=0.7, label=str(sp), edgecolor="black",
                       linewidth=0.4)
        ax.set_xlabel("hier mean bootstrapping_probability")
        ax.set_ylabel("corr mean Pearson r")
        ax.set_title(f"confidence (size ~ {size_col})")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8, framealpha=0.8)
    else:
        ax.text(0.5, 0.5, "no datasets with both methods",
                ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()

    # 2. bar chart: per-row class-level agreement (or histogram if too many)
    ax = axes[0, 1]
    agr_col = next((c for c in ["agreement_class", "agreement_subclass",
                                "agreement_cluster"] if c in df.columns), None)
    if agr_col and df[agr_col].notna().any():
        sub = df.dropna(subset=[agr_col]).sort_values(agr_col)
        if len(sub) > 80:
            ax.hist(sub[agr_col], bins=30, color="steelblue", edgecolor="black",
                    linewidth=0.3)
            ax.set_xlabel(agr_col)
            ax.set_ylabel("count")
            ax.set_xlim(0, 1)
            ax.set_title(f"distribution of {agr_col}")
            ax.grid(axis="y", alpha=0.3)
        else:
            label_col = ("sample" if "slice_id" not in sub.columns
                         else sub["sample"].astype(str) + "/" +
                              sub["slice_id"].astype(str))
            labels = sub[label_col] if isinstance(label_col, str) else label_col
            ax.barh(labels, sub[agr_col], color="steelblue",
                    edgecolor="black", linewidth=0.3)
            ax.set_xlabel(f"{agr_col} (fraction of cells)")
            ax.set_xlim(0, 1)
            ax.set_title(f"hier vs corr — {agr_col}")
            ax.tick_params(axis="y", labelsize=6)
            ax.grid(axis="x", alpha=0.3)
    else:
        ax.text(0.5, 0.5, "no agreement column available",
                ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()

    # 3. distribution of hier mean prob across datasets
    ax = axes[1, 0]
    if "hier_mean_prob" in df.columns and df["hier_mean_prob"].notna().any():
        ax.hist(df["hier_mean_prob"].dropna(), bins=20, color="darkorange",
                edgecolor="black", linewidth=0.3)
        ax.axvline(low_conf_threshold, color="red", linestyle="--",
                   label=f"threshold={low_conf_threshold}")
        ax.set_xlabel("hier_mean_prob")
        ax.set_ylabel("number of datasets")
        ax.set_title("hierarchical confidence across datasets")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
    else:
        ax.set_axis_off()

    # 4. distribution of corr mean r across datasets
    ax = axes[1, 1]
    if "corr_mean_r" in df.columns and df["corr_mean_r"].notna().any():
        ax.hist(df["corr_mean_r"].dropna(), bins=20, color="seagreen",
                edgecolor="black", linewidth=0.3)
        ax.set_xlabel("corr_mean_r")
        ax.set_ylabel("number of datasets")
        ax.set_title("correlation confidence across datasets")
        ax.grid(alpha=0.3)
    else:
        ax.set_axis_off()

    unit = "slices" if "slice_id" in df.columns else "datasets"
    fig.suptitle(f"MMC method comparison — {len(df)} {unit}", fontsize=14)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    logger.info("Wrote %s", out_path)


def make_per_sample_plot(sample: str, meyes_base: Path, out_path: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sample_dir = meyes_base / sample
    hier = load_mapping_csv(sample_dir / HIER_DIR / MAPPING_FILE)
    corr = load_mapping_csv(sample_dir / CORR_DIR / MAPPING_FILE)
    if hier is None and corr is None:
        logger.warning("No mapping CSVs for %s, skipping per-sample plot.", sample)
        return

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    ax = axes[0]
    hier_prob_col = _hier_prob_column(hier) if hier is not None else None
    if hier is not None and hier_prob_col:
        p = pd.to_numeric(hier[hier_prob_col], errors="coerce").dropna()
        ax.hist(p, bins=40, color="darkorange", edgecolor="black", linewidth=0.3)
        ax.axvline(p.median(), color="red", linestyle="--",
                   label=f"median={p.median():.3f}")
        ax.set_title(f"hierarchical — {len(p):,} cells ({hier_prob_col})")
        ax.set_xlabel(hier_prob_col)
        ax.set_ylabel("cells")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
    else:
        ax.set_axis_off()

    ax = axes[1]
    if corr is not None and "correlation_coefficient" in corr.columns:
        r = pd.to_numeric(corr["correlation_coefficient"], errors="coerce").dropna()
        ax.hist(r, bins=40, color="seagreen", edgecolor="black", linewidth=0.3)
        ax.axvline(r.median(), color="red", linestyle="--",
                   label=f"median={r.median():.3f}")
        ax.set_title(f"correlation — {len(r):,} cells")
        ax.set_xlabel("Pearson r")
        ax.set_ylabel("cells")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
    else:
        ax.set_axis_off()

    ax = axes[2]
    confusion_col = next((c for c in ("class_name", "subclass_name", "cluster_name")
                          if hier is not None and corr is not None
                          and c in hier.columns and c in corr.columns), None)
    if confusion_col:
        merged = hier[["cell_id", confusion_col]].merge(
            corr[["cell_id", confusion_col]], on="cell_id", suffixes=("_h", "_c"))
        ct = pd.crosstab(merged[f"{confusion_col}_h"], merged[f"{confusion_col}_c"])
        # Trim to top-15 most common labels for readability.
        top_h = merged[f"{confusion_col}_h"].value_counts().head(15).index
        top_c = merged[f"{confusion_col}_c"].value_counts().head(15).index
        ct = ct.reindex(index=top_h, columns=top_c, fill_value=0)
        norm = ct.div(ct.sum(axis=1).replace(0, 1), axis=0)
        im = ax.imshow(norm.values, aspect="auto", cmap="viridis", vmin=0, vmax=1)
        ax.set_xticks(range(len(ct.columns)))
        ax.set_xticklabels(ct.columns, rotation=90, fontsize=6)
        ax.set_yticks(range(len(ct.index)))
        ax.set_yticklabels(ct.index, fontsize=6)
        ax.set_xlabel(f"correlation {confusion_col}")
        ax.set_ylabel(f"hierarchical {confusion_col}")
        ax.set_title(f"row-normalized confusion ({confusion_col}, top 15)")
        fig.colorbar(im, ax=ax, fraction=0.04)
    else:
        ax.set_axis_off()

    fig.suptitle(f"{sample} — MMC method comparison", fontsize=13)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    logger.info("Wrote %s", out_path)


def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s",
                        stream=sys.stdout)

    meyes_base = Path(args.meyes_base).expanduser()
    out_dir = Path(args.out_dir).expanduser() if args.out_dir else \
        meyes_base / "_mmc_comparison"
    out_dir.mkdir(parents=True, exist_ok=True)

    samples = read_samples(args)
    if not samples:
        raise SystemExit("No samples to process.")

    rows: list[dict] = []
    if args.per_slice:
        for name, species, input_path in samples:
            logger.info("Summarizing %s (per slice)", name)
            rows.extend(summarize_dataset_per_slice(
                name, species, meyes_base, args.low_conf_threshold,
                input_path=input_path, mask_override=args.mask_name,
            ))
        tsv_name = "comparison_per_slice.tsv"
        png_name = "comparison_per_slice.png"
    else:
        for name, species, input_path in samples:
            logger.info("Summarizing %s", name)
            rows.append(summarize_sample(name, species, meyes_base,
                                         args.low_conf_threshold,
                                         input_path=input_path))
        tsv_name = "comparison.tsv"
        png_name = "comparison.png"

    df = pd.DataFrame(rows)
    tsv_path = out_dir / tsv_name
    df.to_csv(tsv_path, sep="\t", index=False)
    logger.info("Wrote %s (%d rows)", tsv_path, len(df))

    make_summary_plot(df, out_dir / png_name, args.low_conf_threshold)

    if args.per_sample_plots or (len(samples) == 1 and not args.per_slice):
        per_dir = out_dir / "per_sample"
        per_dir.mkdir(parents=True, exist_ok=True)
        for name, _, _ in samples:
            make_per_sample_plot(name, meyes_base, per_dir / f"{name}.png")


if __name__ == "__main__":
    main()
