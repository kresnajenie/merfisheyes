#!/usr/bin/env python3
"""Generate chart for S3 single molecule lazy-loading benchmark."""

import csv
import os

import matplotlib.pyplot as plt
import numpy as np

RESULTS_CSV = os.path.join(os.path.dirname(__file__), "../../benchmark_results/s3_sm_benchmark.csv")
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "../../benchmark_results/charts")


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


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    setup_style()

    rows = []
    with open(RESULTS_CSV) as f:
        for row in csv.DictReader(f):
            if row["crashed"] == "True":
                continue
            rows.append(row)

    rows.sort(key=lambda r: int(r["molecules"]))

    molecules = [int(r["molecules"]) for r in rows]
    times = [float(r["load_time_s"]) for r in rows]
    names = [r["name"] for r in rows]

    def fmt_mol(n):
        if n >= 1_000_000_000:
            return f"{n/1e9:.1f}B"
        return f"{n/1e6:.0f}M"

    fig, ax = plt.subplots(figsize=(13, 7))

    # Plot line + markers
    color = "#7c3aed"
    ax.plot(molecules, times, "o-", color=color, linewidth=2.5, markersize=10, zorder=3)

    # Add labels
    for mol, t, name in zip(molecules, times, names):
        short = name.split(" (")[0][:20]
        ax.annotate(f"{t:.1f}s\n{short}",
                    (mol, t), textcoords="offset points",
                    xytext=(0, 18), ha="center", fontsize=8, color="#e0e0e0")

    # Constant line reference
    avg_time = np.mean(times)
    ax.axhline(y=avg_time, color="#22c55e", linestyle="--", alpha=0.5, linewidth=1.5)
    ax.text(molecules[-1] * 0.7, avg_time + 0.3,
            f"avg: {avg_time:.1f}s (manifest-only download)",
            fontsize=10, color="#22c55e", alpha=0.8)

    ax.set_xlabel("Total Molecules in Dataset", fontsize=14)
    ax.set_ylabel("Initial Load Time (seconds)", fontsize=14)
    ax.set_title(
        "MERFISHeyes S3 Lazy-Loading: Initial Load Time vs Dataset Size\n"
        "(Single Molecule, manifest download only — genes loaded on demand)",
        fontsize=15, fontweight="bold", pad=15,
    )
    ax.set_xscale("log")
    ax.grid(alpha=0.3)

    # Format x-axis
    def mol_formatter(x, pos):
        if x >= 1e9:
            return f"{x/1e9:.1f}B"
        if x >= 1e6:
            return f"{x/1e6:.0f}M"
        return f"{x/1e3:.0f}K"

    ax.xaxis.set_major_formatter(plt.FuncFormatter(mol_formatter))
    ax.set_ylim(0, max(times) * 1.6)

    # Add annotation box
    box_text = (
        "S3 Lazy Loading:\n"
        "  1. Download manifest.json.gz (~KB)\n"
        "  2. Parse gene list\n"
        "  3. Auto-select 3 genes\n"
        "  4. Download selected gene files\n"
        "\n"
        "Gene data loaded on-demand\n"
        "when user selects genes"
    )
    ax.text(0.02, 0.97, box_text, transform=ax.transAxes,
            fontsize=9, verticalalignment="top", fontfamily="monospace",
            bbox=dict(boxstyle="round,pad=0.5", facecolor="#16213e",
                      edgecolor="#4a4a6a", alpha=0.9))

    fig.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, "s3_sm_load_scaling.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("  -> s3_sm_load_scaling.png")


if __name__ == "__main__":
    main()
