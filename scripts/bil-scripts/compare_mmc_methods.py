#!/usr/bin/env python3
"""
compare_mmc_methods.py

Compare hierarchical (mmc_output/) vs correlation (mmc_output_corr/) MapMyCells
results across one or many datasets.

Inputs (per sample, on the BIL cluster by default):
  ${MEYES_BASE}/<sample>/mmc_output/mapping_output.csv       (hierarchical)
  ${MEYES_BASE}/<sample>/mmc_output_corr/mapping_output.csv  (correlation)

Outputs (in --out-dir):
  comparison.tsv          one row per dataset with summary metrics
  comparison.png          multi-panel summary across datasets
  per_sample/<sample>.png per-cell histograms + class-level confusion matrix
                          (only when --per-sample-plots is set, or when run
                          on a single sample)

Usage:
  # all rows in samples.csv
  python compare_mmc_methods.py samples.csv

  # one sample (also writes per_sample/<sample>.png)
  python compare_mmc_methods.py --sample ace-dip-use

  # custom roots
  python compare_mmc_methods.py samples.csv \
      --meyes-base /bil/data/meyes \
      --out-dir   /bil/data/meyes/_mmc_comparison
"""

import argparse
import csv
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

logger = logging.getLogger("compare_mmc")

DEFAULT_MEYES_BASE = "/bil/data/meyes"
HIER_DIR = "mmc_output"
CORR_DIR = "mmc_output_corr"
MAPPING_FILE = "mapping_output.csv"

LEVEL_COLS = ["class_name", "subclass_name", "supertype_name",
              "cluster_name", "cluster_label"]


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
    return p.parse_args()


def read_samples(args):
    """Returns list of (sample_name, species_or_None)."""
    rows = []
    species_lookup = {}

    if args.species_map:
        sm = pd.read_csv(args.species_map)
        if "sample_name" not in sm.columns or "species" not in sm.columns:
            raise ValueError("--species-map needs columns sample_name,species")
        species_lookup = dict(zip(sm["sample_name"].astype(str),
                                  sm["species"].astype(str)))

    if args.sample:
        for s in args.sample:
            rows.append((s, species_lookup.get(s)))
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
            rows.append((name, sp))
    return rows


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


def summarize_sample(sample: str, species: str | None, meyes_base: Path,
                     low_conf_threshold: float) -> dict:
    sample_dir = meyes_base / sample
    hier_path = sample_dir / HIER_DIR / MAPPING_FILE
    corr_path = sample_dir / CORR_DIR / MAPPING_FILE

    hier = load_mapping_csv(hier_path)
    corr = load_mapping_csv(corr_path)

    row: dict = {"sample": sample, "species": species or ""}

    if hier is not None:
        row["hier_path"] = str(hier_path)
        row["n_cells_hier"] = len(hier)
        if "bootstrapping_probability" in hier.columns:
            p = pd.to_numeric(hier["bootstrapping_probability"], errors="coerce").dropna()
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
        row["n_cells"] = int(min(len(hier), len(corr)))
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
    sub = df.dropna(subset=["hier_mean_prob", "corr_mean_r"])
    if not sub.empty:
        species_vals = sub["species"].fillna("").replace("", "unknown")
        for sp, group in sub.groupby(species_vals):
            sizes = np.clip(group["n_cells"].fillna(1).astype(float) / 200, 10, 400)
            ax.scatter(group["hier_mean_prob"], group["corr_mean_r"],
                       s=sizes, alpha=0.7, label=str(sp), edgecolor="black",
                       linewidth=0.4)
        ax.set_xlabel("hier mean bootstrapping_probability")
        ax.set_ylabel("corr mean Pearson r")
        ax.set_title("per-dataset confidence (size ~ n_cells)")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8, framealpha=0.8)
    else:
        ax.text(0.5, 0.5, "no datasets with both methods",
                ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()

    # 2. bar chart: per-dataset class-level agreement
    ax = axes[0, 1]
    agr_col = next((c for c in ["agreement_class", "agreement_subclass",
                                "agreement_cluster"] if c in df.columns), None)
    if agr_col and df[agr_col].notna().any():
        sub = df.dropna(subset=[agr_col]).sort_values(agr_col)
        ax.barh(sub["sample"], sub[agr_col], color="steelblue", edgecolor="black",
                linewidth=0.3)
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
    if df["hier_mean_prob"].notna().any():
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
    if df["corr_mean_r"].notna().any():
        ax.hist(df["corr_mean_r"].dropna(), bins=20, color="seagreen",
                edgecolor="black", linewidth=0.3)
        ax.set_xlabel("corr_mean_r")
        ax.set_ylabel("number of datasets")
        ax.set_title("correlation confidence across datasets")
        ax.grid(alpha=0.3)
    else:
        ax.set_axis_off()

    fig.suptitle(f"MMC method comparison — {len(df)} datasets", fontsize=14)
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
    if hier is not None and "bootstrapping_probability" in hier.columns:
        p = pd.to_numeric(hier["bootstrapping_probability"], errors="coerce").dropna()
        ax.hist(p, bins=40, color="darkorange", edgecolor="black", linewidth=0.3)
        ax.axvline(p.median(), color="red", linestyle="--",
                   label=f"median={p.median():.3f}")
        ax.set_title(f"hierarchical — {len(p):,} cells")
        ax.set_xlabel("bootstrapping_probability")
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

    rows = []
    for name, species in samples:
        logger.info("Summarizing %s", name)
        rows.append(summarize_sample(name, species, meyes_base,
                                     args.low_conf_threshold))

    df = pd.DataFrame(rows)
    tsv_path = out_dir / "comparison.tsv"
    df.to_csv(tsv_path, sep="\t", index=False)
    logger.info("Wrote %s (%d rows)", tsv_path, len(df))

    make_summary_plot(df, out_dir / "comparison.png", args.low_conf_threshold)

    if args.per_sample_plots or len(samples) == 1:
        per_dir = out_dir / "per_sample"
        per_dir.mkdir(parents=True, exist_ok=True)
        for name, _ in samples:
            make_per_sample_plot(name, meyes_base, per_dir / f"{name}.png")


if __name__ == "__main__":
    main()
