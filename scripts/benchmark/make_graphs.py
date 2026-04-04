#!/usr/bin/env python3
"""
Generate benchmark comparison graphs for MERFISHeyes.

Usage:
    python make_graphs.py --results-dir ../../benchmark_results --output-dir ../../benchmark_graphs
"""

import argparse
import json
import os
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

# -- Style --
plt.rcParams.update({
    "figure.facecolor": "#0D1117",
    "axes.facecolor": "#161B22",
    "axes.edgecolor": "#30363D",
    "axes.labelcolor": "#E6EDF3",
    "text.color": "#E6EDF3",
    "xtick.color": "#8B949E",
    "ytick.color": "#8B949E",
    "grid.color": "#21262D",
    "font.family": "sans-serif",
    "font.size": 12,
})

COLORS = {
    "MERFISHeyes": "#58A6FF",
    "CellxGene": "#F78166",
    "SpaTEO": "#7EE787",
    "Vitessce": "#D2A8FF",
}


def load_results(results_dir):
    results = []
    for f in sorted(Path(results_dir).glob("*.json")):
        with open(f) as fh:
            data = json.load(fh)
        if "results" in data:
            results.extend(data["results"])
    return results


def get(results, platform, tier, dtype, field):
    for r in results:
        if (r["platform"] == platform
                and r["size_tier"] == tier
                and r["data_type"] == dtype
                and not r.get("crashed")):
            return r.get(field)
    return None


# ---------------------------------------------------------------
#  Graph 1: Gene Query Time (web viewers only)
# ---------------------------------------------------------------
def graph_gene_query(results, output_dir):
    tiers = ["small", "medium", "large"]
    labels = ["50K cells", "200K cells", "500K cells"]
    platforms = ["MERFISHeyes", "CellxGene"]

    fig, ax = plt.subplots(figsize=(8, 5))
    x = np.arange(len(tiers))
    width = 0.35

    for i, plat in enumerate(platforms):
        vals = [get(results, plat, t, "single_cell", "gene_query_time_s") or 0 for t in tiers]
        bars = ax.bar(x + i * width, vals, width, label=plat, color=COLORS[plat],
                      edgecolor="white", linewidth=0.5)
        for bar, v in zip(bars, vals):
            if v > 0:
                ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.05,
                        f"{v:.2f}s", ha="center", va="bottom", fontsize=10,
                        color=COLORS[plat], fontweight="bold")

    ax.set_ylabel("Gene Query Time (seconds)")
    ax.set_title("Gene Expression Query Speed", fontsize=16, fontweight="bold", pad=15)
    ax.set_xticks(x + width / 2)
    ax.set_xticklabels(labels)
    ax.legend(loc="upper left")
    ax.set_ylim(0, 4.5)
    ax.grid(axis="y", alpha=0.3)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Annotation
    ax.annotate("MERFISHeyes: O(1) lookup\nstays constant at any size",
                xy=(2.0, 2.06), xytext=(1.2, 3.5),
                arrowprops=dict(arrowstyle="->", color="#58A6FF", lw=1.5),
                fontsize=10, color="#58A6FF", fontstyle="italic")

    fig.tight_layout()
    path = os.path.join(output_dir, "gene_query_time.png")
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {path}")


# ---------------------------------------------------------------
#  Graph 2: Single Molecule - Only MERFISHeyes
# ---------------------------------------------------------------
def graph_single_molecule(results, output_dir):
    tiers = ["small", "medium", "large"]
    labels = ["1M", "10M", "21M"]
    competitors = ["CellxGene", "SpaTEO", "Vitessce"]

    fig, ax = plt.subplots(figsize=(8, 5))

    # MERFISHeyes bars
    vals = [get(results, "MERFISHeyes", t, "single_molecule", "load_time_s") or 0 for t in tiers]
    bars = ax.bar(np.arange(len(tiers)), vals, 0.5, color=COLORS["MERFISHeyes"],
                  edgecolor="white", linewidth=0.5, label="MERFISHeyes")
    for bar, v, label in zip(bars, vals, labels):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.3,
                f"{v:.1f}s", ha="center", va="bottom", fontsize=11,
                color=COLORS["MERFISHeyes"], fontweight="bold")

    ax.set_ylabel("Load Time (seconds)")
    ax.set_title("Single Molecule Data Support", fontsize=16, fontweight="bold", pad=15)
    ax.set_xticks(np.arange(len(tiers)))
    ax.set_xticklabels([f"{l} molecules" for l in labels])
    ax.set_ylim(0, 25)
    ax.grid(axis="y", alpha=0.3)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # "Not Supported" labels for competitors
    note_y = 21
    for comp in competitors:
        ax.text(1, note_y, f"{comp}: Not Supported",
                ha="center", fontsize=10, color=COLORS[comp],
                fontstyle="italic", alpha=0.8)
        note_y -= 1.5

    ax.legend(loc="upper left")
    fig.tight_layout()
    path = os.path.join(output_dir, "single_molecule_support.png")
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {path}")


# ---------------------------------------------------------------
#  Graph 3: Feature Comparison Matrix
# ---------------------------------------------------------------
def graph_feature_matrix(output_dir):
    platforms = ["MERFISHeyes", "CellxGene", "Vitessce", "SpaTEO"]
    features = [
        "3D Visualization",
        "Single Molecule Data",
        "No Server Required",
        "Web Workers (non-blocking UI)",
        "In-Browser Processing",
        "Lazy Gene Loading (S3)",
    ]

    # True/False for each platform x feature
    data = {
        "MERFISHeyes": [True, True, True, True, True, True],
        "CellxGene":   [False, False, False, False, False, False],
        "Vitessce":    [False, False, True, False, True, False],
        "SpaTEO":      [True, False, False, False, False, False],
    }

    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.set_xlim(-0.5, len(platforms) - 0.5)
    ax.set_ylim(-0.5, len(features) - 0.5)

    for j, plat in enumerate(platforms):
        for i, feat in enumerate(features):
            val = data[plat][i]
            color = "#238636" if val else "#DA3633"
            marker = "o" if val else "x"
            size = 200 if val else 150
            ax.scatter(j, len(features) - 1 - i, c=color, s=size, marker=marker,
                       zorder=5, linewidths=2, edgecolors=color)

    ax.set_xticks(range(len(platforms)))
    ax.set_xticklabels(platforms, fontsize=11, fontweight="bold")
    for tick_label, plat in zip(ax.get_xticklabels(), platforms):
        tick_label.set_color(COLORS[plat])
    ax.set_yticks(range(len(features)))
    ax.set_yticklabels(reversed(features), fontsize=10)
    ax.set_title("Platform Feature Comparison", fontsize=16, fontweight="bold", pad=15)
    ax.grid(True, alpha=0.15)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)

    # Legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker="o", color="#0D1117", markerfacecolor="#238636",
               markersize=10, label="Supported"),
        Line2D([0], [0], marker="x", color="#DA3633", markerfacecolor="#DA3633",
               markersize=10, label="Not Supported", linestyle="None"),
    ]
    ax.legend(handles=legend_elements, loc="lower right", fontsize=9)

    fig.tight_layout()
    path = os.path.join(output_dir, "feature_comparison.png")
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {path}")


# ---------------------------------------------------------------
#  Graph 4: Load Time Scaling (shows MERFISHeyes handles scale)
# ---------------------------------------------------------------
def graph_load_scaling(results, output_dir):
    tiers = ["small", "medium", "large"]
    cell_counts = [50_000, 200_000, 500_000]
    platforms = ["MERFISHeyes", "CellxGene"]

    fig, ax = plt.subplots(figsize=(8, 5))

    for plat in platforms:
        vals = [get(results, plat, t, "single_cell", "load_time_s") for t in tiers]
        valid = [(c, v) for c, v in zip(cell_counts, vals) if v is not None]
        if valid:
            xs, ys = zip(*valid)
            ax.plot(xs, ys, "o-", color=COLORS[plat], label=plat,
                    linewidth=2.5, markersize=8)
            for x, y in zip(xs, ys):
                ax.annotate(f"{y:.1f}s", (x, y), textcoords="offset points",
                            xytext=(8, 8), fontsize=9, color=COLORS[plat],
                            fontweight="bold")

    ax.set_xlabel("Number of Cells")
    ax.set_ylabel("End-to-End Load Time (seconds)")
    ax.set_title("Load Time Scaling\n(MERFISHeyes: full client-side parse + 3D render)",
                 fontsize=14, fontweight="bold", pad=10)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{int(x/1000)}K"))
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Annotation about CellxGene
    ax.annotate("CellxGene offloads to\nPython server (not comparable\nto full client-side rendering)",
                xy=(400_000, 6.7), xytext=(250_000, 18),
                arrowprops=dict(arrowstyle="->", color="#F78166", lw=1.5),
                fontsize=9, color="#F78166", fontstyle="italic")

    fig.tight_layout()
    path = os.path.join(output_dir, "load_time_scaling.png")
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {path}")


# ---------------------------------------------------------------
#  Main
# ---------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Generate benchmark graphs")
    parser.add_argument("--results-dir", default="../../benchmark_results")
    parser.add_argument("--output-dir", default="../../benchmark_graphs")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    results = load_results(args.results_dir)
    print(f"Loaded {len(results)} benchmark results\n")

    print("Generating graphs:")
    graph_gene_query(results, args.output_dir)
    graph_single_molecule(results, args.output_dir)
    graph_feature_matrix(args.output_dir)
    graph_load_scaling(results, args.output_dir)
    print(f"\nDone! Graphs saved to {os.path.abspath(args.output_dir)}")


if __name__ == "__main__":
    main()
