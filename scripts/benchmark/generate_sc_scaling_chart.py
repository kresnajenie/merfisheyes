#!/usr/bin/env python3
"""Generate scaling chart: cells on x-axis, load time on y-axis, one line per format."""

import csv
import os

import matplotlib.pyplot as plt
import numpy as np

RESULTS_CSV = os.path.join(os.path.dirname(__file__), "../../benchmark_results/sc_load_benchmark.csv")
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "../../benchmark_results/charts")


def load_data():
    rows = []
    with open(RESULTS_CSV) as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["crashed"] == "True":
                continue
            rows.append(row)
    return rows


def setup_style():
    plt.rcParams.update({
        "figure.facecolor": "#1a1a2e",
        "axes.facecolor": "#16213e",
        "axes.edgecolor": "#e0e0e0",
        "axes.labelcolor": "#e0e0e0",
        "text.color": "#e0e0e0",
        "xtick.color": "#e0e0e0",
        "ytick.color": "#e0e0e0",
        "grid.color": "#2a2a4a",
        "grid.alpha": 0.5,
        "font.family": "sans-serif",
        "font.size": 12,
    })


COLORS = {
    "h5ad": "#3b82f6",      # blue
    "xenium": "#10b981",    # green
    "merscope": "#f59e0b",  # amber
}
MARKERS = {
    "h5ad": "o",
    "xenium": "s",
    "merscope": "D",
}
LABELS = {
    "h5ad": "H5AD",
    "xenium": "Xenium",
    "merscope": "MERSCOPE",
}


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    setup_style()
    rows = load_data()

    # Group by filetype
    by_type = {}
    for row in rows:
        ft = row["filetype"]
        if ft not in by_type:
            by_type[ft] = []
        by_type[ft].append(row)

    # Sort each group by n_cells
    for ft in by_type:
        by_type[ft].sort(key=lambda r: int(r["n_cells"]))

    # --- Chart 1: Cells vs Load Time ---
    fig, ax = plt.subplots(figsize=(12, 7))

    for ft in ["h5ad", "xenium", "merscope"]:
        if ft not in by_type:
            continue
        data = by_type[ft]
        cells = [int(r["n_cells"]) for r in data]
        times = [float(r["load_time_s"]) for r in data]

        ax.plot(cells, times,
                color=COLORS[ft], marker=MARKERS[ft], markersize=8,
                linewidth=2.5, label=LABELS[ft], alpha=0.9)

        # Add value labels
        for c, t in zip(cells, times):
            ax.annotate(f"{t:.1f}s", (c, t),
                        textcoords="offset points", xytext=(0, 12),
                        ha="center", fontsize=8, color=COLORS[ft])

    ax.set_xlabel("Number of Cells", fontsize=14)
    ax.set_ylabel("Load Time (seconds)", fontsize=14)
    ax.set_title("MERFISHeyes Single Cell Load Time vs Dataset Size\n(500 genes, local browser upload)",
                 fontsize=16, fontweight="bold", pad=15)
    ax.legend(fontsize=12, framealpha=0.8, facecolor="#1a1a2e", edgecolor="#444")
    ax.grid(alpha=0.3)

    # Format x-axis with K/M labels
    def cell_formatter(x, pos):
        if x >= 1_000_000:
            return f"{x/1_000_000:.0f}M"
        if x >= 1_000:
            return f"{x/1_000:.0f}K"
        return str(int(x))

    ax.xaxis.set_major_formatter(plt.FuncFormatter(cell_formatter))

    fig.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, "sc_load_scaling.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("  -> sc_load_scaling.png")

    # --- Chart 2: File Size vs Load Time ---
    fig, ax = plt.subplots(figsize=(12, 7))

    for ft in ["h5ad", "xenium", "merscope"]:
        if ft not in by_type:
            continue
        data = by_type[ft]
        sizes = [float(r["file_size_mb"]) for r in data]
        times = [float(r["load_time_s"]) for r in data]

        ax.plot(sizes, times,
                color=COLORS[ft], marker=MARKERS[ft], markersize=8,
                linewidth=2.5, label=LABELS[ft], alpha=0.9)

        for s, t in zip(sizes, times):
            label = f"{s:.0f}MB" if s < 1024 else f"{s/1024:.1f}GB"
            ax.annotate(f"{t:.1f}s", (s, t),
                        textcoords="offset points", xytext=(0, 12),
                        ha="center", fontsize=8, color=COLORS[ft])

    ax.set_xlabel("File Size (MB)", fontsize=14)
    ax.set_ylabel("Load Time (seconds)", fontsize=14)
    ax.set_title("MERFISHeyes Single Cell Load Time vs File Size\n(500 genes, local browser upload)",
                 fontsize=16, fontweight="bold", pad=15)
    ax.legend(fontsize=12, framealpha=0.8, facecolor="#1a1a2e", edgecolor="#444")
    ax.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, "sc_filesize_scaling.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("  -> sc_filesize_scaling.png")

    print(f"\nCharts saved to {OUTPUT_DIR}/")


if __name__ == "__main__":
    main()
