#!/usr/bin/env python3
"""Chart the ingestion benchmark in dataset-results.csv.

Two panels, each answering a question the raw table makes you squint at:
  left  — where does wall-clock actually go? (upload vs processing)
  right — does processing time track input size, or format?

    python3 docs/server-side-ingestion/benchmarks/plot_results.py
"""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

HERE = Path(__file__).resolve().parent
FMT_COLOR = {"h5ad": "#5EA2EF", "merscope": "#FF1CF7"}


def load():
    rows = []
    with open(HERE / "dataset-results.csv") as f:
        for r in csv.DictReader(f):
            # The CSV-annotated re-run is the same input as its plain
            # counterpart; including it would double-plot one dataset.
            if r["annotation_csv"] == "yes":
                continue
            r["input_mb"] = float(r["input_mb"])
            r["e2e_upload_s"] = float(r["e2e_upload_s"])
            r["e2e_process_s"] = float(r["e2e_process_s"])
            r["local_process_s"] = float(r["local_process_s"])
            rows.append(r)
    return sorted(rows, key=lambda r: r["input_mb"])


def short(name):
    return name.replace(".h5ad", "").replace("_labelled", "").replace("_labeled", "")[:26]


def main():
    rows = load()
    labels = [f'{short(r["dataset"])}\n{r["input_mb"]:.0f} MB' for r in rows]
    y = range(len(rows))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 7))

    # ── Left: upload vs processing ───────────────────────────────────────
    up = [r["e2e_upload_s"] for r in rows]
    proc = [r["e2e_process_s"] for r in rows]

    ax1.barh(y, up, color="#94a3b8", label="Upload")
    ax1.barh(y, proc, left=up, color=[FMT_COLOR[r["format"]] for r in rows],
             label="Processing")
    for i, (u, p) in enumerate(zip(up, proc)):
        ax1.text(u + p + 12, i, f"{u + p:.0f}s", va="center", fontsize=9,
                 color="#475569")

    ax1.set_yticks(list(y))
    ax1.set_yticklabels(labels, fontsize=8)
    ax1.set_xlabel("seconds")
    ax1.set_title("Where the wall-clock goes\nUpload dominates everything except MERSCOPE",
                  fontsize=11, loc="left")
    ax1.legend(handles=[
        Patch(color="#94a3b8", label="Upload (~6 MB/s)"),
        Patch(color=FMT_COLOR["h5ad"], label="Processing — h5ad"),
        Patch(color=FMT_COLOR["merscope"], label="Processing — MERSCOPE"),
    ], loc="lower right", fontsize=9)
    ax1.grid(axis="x", alpha=0.25)
    ax1.set_axisbelow(True)

    # ── Right: processing time vs size, local vs Fargate ─────────────────
    for fmt in ("h5ad", "merscope"):
        sel = [r for r in rows if r["format"] == fmt]
        ax2.scatter([r["input_mb"] for r in sel], [r["e2e_process_s"] for r in sel],
                    s=90, color=FMT_COLOR[fmt], label=f"{fmt} — Fargate", zorder=3)
        ax2.scatter([r["input_mb"] for r in sel], [r["local_process_s"] for r in sel],
                    s=90, facecolors="none", edgecolors=FMT_COLOR[fmt],
                    linewidths=1.6, label=f"{fmt} — local", zorder=3)
        # Join each dataset's two measurements so the gap is legible.
        for r in sel:
            ax2.plot([r["input_mb"], r["input_mb"]],
                     [r["local_process_s"], r["e2e_process_s"]],
                     color=FMT_COLOR[fmt], alpha=0.35, lw=1.2, zorder=2)

    # Headroom so the callouts below don't run into the title.
    ax2.set_ylim(top=max(r["e2e_process_s"] for r in rows) * 1.28)

    worst = max(rows, key=lambda r: r["e2e_process_s"])
    ax2.annotate(
        f'{short(worst["dataset"])}\n{worst["e2e_process_s"]:.0f}s for only {worst["input_mb"]:.0f} MB',
        xy=(worst["input_mb"], worst["e2e_process_s"]),
        xytext=(worst["input_mb"] + 900, worst["e2e_process_s"] - 40),
        fontsize=9, color="#475569",
        arrowprops=dict(arrowstyle="->", color="#94a3b8"),
    )
    biggest = max(rows, key=lambda r: r["input_mb"])
    ax2.annotate(
        f'{biggest["input_mb"]:.0f} MB h5ad\nonly {biggest["e2e_process_s"]:.0f}s',
        xy=(biggest["input_mb"], biggest["e2e_process_s"]),
        xytext=(biggest["input_mb"] - 1450, biggest["e2e_process_s"] + 170),
        fontsize=9, color="#475569",
        arrowprops=dict(arrowstyle="->", color="#94a3b8"),
    )

    ax2.set_xlabel("input size (MB)")
    ax2.set_ylabel("processing time (s)")
    ax2.set_title("Processing time tracks format, not size\nFilled = Fargate (4 vCPU), hollow = local",
                  fontsize=11, loc="left")
    ax2.legend(fontsize=9)
    ax2.grid(alpha=0.25)
    ax2.set_axisbelow(True)

    fig.suptitle(
        "Server-side ingestion — 10 single-cell datasets, 135 MB to 4 GB, all COMPLETE",
        fontsize=13, y=0.99,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    out = HERE / "dataset-results.png"
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
