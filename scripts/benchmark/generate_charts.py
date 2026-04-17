#!/usr/bin/env python3
"""Generate benchmark charts from real_data_consolidated.json."""

import json
import os

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

RESULTS_PATH = os.path.join(os.path.dirname(__file__), "../../benchmark_results/real_data_consolidated.json")
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "../../benchmark_results/charts")


def load_data():
    with open(RESULTS_PATH) as f:
        return json.load(f)["results"]


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
        "font.size": 11,
    })


# ── Colors ──
SM_COLOR = "#7c3aed"      # purple
SC_MERSCOPE = "#3b82f6"   # blue
SC_XENIUM = "#10b981"     # green
CRASH_COLOR = "#ef4444"   # red
SKIP_COLOR = "#6b7280"    # grey


def chart_load_times(results):
    """Bar chart: load time by dataset, with crash/skip indicators."""
    fig, ax = plt.subplots(figsize=(14, 7))

    # Separate by category
    sm = [r for r in results if r["data_type"] == "single_molecule"]
    sc_m = [r for r in results if r["data_type"] == "single_cell" and "MERSCOPE" in r["dataset_name"]]
    sc_x = [r for r in results if r["data_type"] == "single_cell" and "Xenium" in r["dataset_name"]]

    all_ds = sm + sc_m + sc_x
    labels = []
    times = []
    colors = []
    hatches = []

    for r in all_ds:
        # Short label
        name = r["dataset_name"]
        size_gb = r["file_size_mb"] / 1024
        label = name.split("(")[0].strip()
        if size_gb >= 1:
            label += f"\n({size_gb:.1f} GB)"
        else:
            label += f"\n({r['file_size_mb']:.0f} MB)"
        labels.append(label)

        if r["crashed"]:
            # Show time before crash, or placeholder
            t = r["load_time_s"] if r["load_time_s"] else 60
            times.append(t)
            colors.append(CRASH_COLOR)
            hatches.append("//")
        elif r["load_time_s"] is None:
            times.append(0)
            colors.append(SKIP_COLOR)
            hatches.append("xx")
        else:
            times.append(r["load_time_s"])
            if r in sm:
                colors.append(SM_COLOR)
            elif r in sc_m:
                colors.append(SC_MERSCOPE)
            else:
                colors.append(SC_XENIUM)
            hatches.append("")

    x = np.arange(len(labels))

    # Compute y-axis cap before drawing
    reasonable_times = [t for t in times if t < 200]
    cap = max(reasonable_times) * 1.4 if reasonable_times else 120

    # Clamp bar heights to cap for display
    display_times = [min(t, cap * 0.95) for t in times]
    bars = ax.bar(x, display_times, color=colors, edgecolor="white", linewidth=0.5, width=0.7)

    # Add hatching for crash/skip
    for bar, h in zip(bars, hatches):
        if h:
            bar.set_hatch(h)
            bar.set_edgecolor("white")

    # Value labels on bars
    for i, (bar, t, r) in enumerate(zip(bars, times, all_ds)):
        bx = bar.get_x() + bar.get_width() / 2
        bh = bar.get_height()
        if r["crashed"] and r["load_time_s"] and r["load_time_s"] > 100:
            # Timeout - show inside truncated bar
            ax.text(bx, bh * 0.5, f"{t/60:.0f}min\nTIMEOUT",
                    ha="center", va="center", fontsize=8, color="white", fontweight="bold")
        elif r["crashed"]:
            ax.text(bx, bh + cap * 0.03, f"{t:.0f}s\nOOM",
                    ha="center", va="bottom", fontsize=8, color=CRASH_COLOR, fontweight="bold")
        elif r["load_time_s"] is None:
            ax.text(bx, cap * 0.05, "SKIPPED",
                    ha="center", va="bottom", fontsize=8, color="white", fontweight="bold")
        else:
            ax.text(bx, bh + cap * 0.02, f"{t:.1f}s",
                    ha="center", va="bottom", fontsize=9, color="white")

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8, rotation=0)
    ax.set_ylabel("Load Time (seconds)")
    ax.set_title("MERFISHeyes Load Time by Dataset", fontsize=16, fontweight="bold", pad=15)
    ax.set_ylim(0, cap)
    ax.grid(axis="y", alpha=0.3)


    # Legend
    legend_patches = [
        mpatches.Patch(facecolor=SM_COLOR, label="Single Molecule"),
        mpatches.Patch(facecolor=SC_MERSCOPE, label="Single Cell MERSCOPE"),
        mpatches.Patch(facecolor=SC_XENIUM, label="Single Cell Xenium"),
        mpatches.Patch(facecolor=CRASH_COLOR, hatch="//", edgecolor="white", label="Crashed / Timeout"),
        mpatches.Patch(facecolor=SKIP_COLOR, hatch="xx", edgecolor="white", label="Skipped"),
    ]
    ax.legend(handles=legend_patches, loc="upper left", framealpha=0.8,
              facecolor="#1a1a2e", edgecolor="#444")

    fig.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, "load_times.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("  -> load_times.png")


def chart_memory(results):
    """Bar chart: peak memory usage."""
    fig, ax = plt.subplots(figsize=(10, 6))

    # Only datasets with memory data
    mem_ds = [r for r in results if r["memory_peak_mb"] is not None]

    labels = []
    memory = []
    colors = []

    for r in mem_ds:
        name = r["dataset_name"].split("(")[0].strip()
        size_gb = r["file_size_mb"] / 1024
        if size_gb >= 1:
            labels.append(f"{name}\n({size_gb:.1f} GB)")
        else:
            labels.append(f"{name}\n({r['file_size_mb']:.0f} MB)")
        memory.append(r["memory_peak_mb"])

        if r["data_type"] == "single_molecule":
            colors.append(SM_COLOR)
        elif "MERSCOPE" in r["dataset_name"]:
            colors.append(SC_MERSCOPE)
        else:
            colors.append(SC_XENIUM)

    x = np.arange(len(labels))
    bars = ax.bar(x, memory, color=colors, edgecolor="white", linewidth=0.5, width=0.6)

    for bar, m in zip(bars, memory):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 10,
                f"{m:.0f} MB", ha="center", va="bottom", fontsize=9, color="white")

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("Peak JS Heap Memory (MB)")
    ax.set_title("MERFISHeyes Memory Usage by Dataset", fontsize=16, fontweight="bold", pad=15)
    ax.grid(axis="y", alpha=0.3)

    legend_patches = [
        mpatches.Patch(facecolor=SM_COLOR, label="Single Molecule"),
        mpatches.Patch(facecolor=SC_MERSCOPE, label="Single Cell MERSCOPE"),
        mpatches.Patch(facecolor=SC_XENIUM, label="Single Cell Xenium"),
    ]
    ax.legend(handles=legend_patches, loc="upper left", framealpha=0.8,
              facecolor="#1a1a2e", edgecolor="#444")

    fig.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, "memory_usage.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("  -> memory_usage.png")


def chart_fps(results):
    """Horizontal bar chart: FPS after load."""
    fig, ax = plt.subplots(figsize=(10, 5))

    fps_ds = [r for r in results if r["fps_after_load"] is not None]

    labels = []
    fps = []
    colors = []

    for r in fps_ds:
        name = r["dataset_name"].split("(")[0].strip()
        labels.append(name)
        fps.append(r["fps_after_load"])

        if r["data_type"] == "single_molecule":
            colors.append(SM_COLOR)
        elif "MERSCOPE" in r["dataset_name"]:
            colors.append(SC_MERSCOPE)
        else:
            colors.append(SC_XENIUM)

    y = np.arange(len(labels))
    bars = ax.barh(y, fps, color=colors, edgecolor="white", linewidth=0.5, height=0.5)

    for bar, f in zip(bars, fps):
        ax.text(bar.get_width() + 0.5, bar.get_y() + bar.get_height()/2,
                f"{f:.0f} FPS", ha="left", va="center", fontsize=10, color="white")

    # Reference lines
    ax.axvline(x=30, color="#22c55e", linestyle="--", alpha=0.5, linewidth=1)
    ax.axvline(x=10, color="#eab308", linestyle="--", alpha=0.5, linewidth=1)
    ax.text(30.5, len(labels) - 0.3, "30 FPS (smooth)", fontsize=7, color="#22c55e", alpha=0.7)
    ax.text(10.5, len(labels) - 0.3, "10 FPS (usable)", fontsize=7, color="#eab308", alpha=0.7)

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel("Frames Per Second")
    ax.set_title("MERFISHeyes Rendering Performance (FPS)", fontsize=16, fontweight="bold", pad=15)
    ax.set_xlim(0, max(fps) * 1.3)
    ax.grid(axis="x", alpha=0.3)
    ax.invert_yaxis()

    fig.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, "fps.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("  -> fps.png")


def chart_sm_scaling(results):
    """Line chart: single molecule scaling (load time & memory vs molecules)."""
    sm = [r for r in results if r["data_type"] == "single_molecule" and r["load_time_s"] is not None and not r["crashed"]]

    if len(sm) < 2:
        print("  -> sm_scaling.png SKIPPED (not enough data)")
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    molecules = [r["n_items"] / 1e6 for r in sm]
    load_times = [r["load_time_s"] for r in sm]
    memory = [r["memory_peak_mb"] for r in sm]

    # Load time scaling
    ax1.plot(molecules, load_times, "o-", color=SM_COLOR, linewidth=2, markersize=8)
    for m, t in zip(molecules, load_times):
        ax1.annotate(f"{t:.0f}s", (m, t), textcoords="offset points",
                     xytext=(0, 12), ha="center", fontsize=9, color="white")

    # Extrapolation line to 67M and 121M
    if len(sm) >= 2:
        # Linear extrapolation
        coeffs = np.polyfit(molecules, load_times, 1)
        extrap_x = [67.5, 120.7]
        extrap_y = [np.polyval(coeffs, x) for x in extrap_x]
        ax1.plot(extrap_x, extrap_y, "x", color=CRASH_COLOR, markersize=12, markeredgewidth=2)
        ax1.plot([molecules[-1], extrap_x[0]], [load_times[-1], extrap_y[0]],
                 "--", color=CRASH_COLOR, alpha=0.4)
        for ex, ey, label in zip(extrap_x, extrap_y, ["67M\n(timed out)", "121M\n(skipped)"]):
            ax1.annotate(f"~{ey/60:.0f}m\n{label}", (ex, ey), textcoords="offset points",
                         xytext=(0, 12), ha="center", fontsize=8, color=CRASH_COLOR)

    ax1.set_xlabel("Molecules (millions)")
    ax1.set_ylabel("Load Time (seconds)")
    ax1.set_title("Load Time Scaling", fontsize=13, fontweight="bold")
    ax1.grid(alpha=0.3)

    # Memory scaling
    ax2.plot(molecules, memory, "s-", color="#f59e0b", linewidth=2, markersize=8)
    for m, mem in zip(molecules, memory):
        ax2.annotate(f"{mem:.0f} MB", (m, mem), textcoords="offset points",
                     xytext=(0, 12), ha="center", fontsize=9, color="white")

    # Extrapolation
    if len(sm) >= 2:
        coeffs_mem = np.polyfit(molecules, memory, 1)
        extrap_mem = [np.polyval(coeffs_mem, x) for x in extrap_x]
        ax2.plot(extrap_x, extrap_mem, "x", color=CRASH_COLOR, markersize=12, markeredgewidth=2)
        ax2.plot([molecules[-1], extrap_x[0]], [memory[-1], extrap_mem[0]],
                 "--", color=CRASH_COLOR, alpha=0.4)
        for ex, ey in zip(extrap_x, extrap_mem):
            ax2.annotate(f"~{ey/1024:.1f} GB", (ex, ey), textcoords="offset points",
                         xytext=(0, 12), ha="center", fontsize=8, color=CRASH_COLOR)

    ax2.set_xlabel("Molecules (millions)")
    ax2.set_ylabel("Peak Memory (MB)")
    ax2.set_title("Memory Scaling", fontsize=13, fontweight="bold")
    ax2.grid(alpha=0.3)

    fig.suptitle("Single Molecule: Scaling with Dataset Size", fontsize=16, fontweight="bold", y=1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, "sm_scaling.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("  -> sm_scaling.png")


def chart_summary_dashboard(results):
    """Combined dashboard with all key metrics."""
    fig = plt.figure(figsize=(16, 10))
    fig.suptitle("MERFISHeyes Benchmark Dashboard — Real Test Data",
                 fontsize=18, fontweight="bold", y=0.98)

    # Datasets that loaded successfully
    ok = [r for r in results if r["load_time_s"] is not None and not r["crashed"]]

    # ── Top left: Load time bars ──
    ax1 = fig.add_subplot(2, 2, 1)
    names = [r["dataset_name"].split("(")[0].strip()[:20] for r in ok]
    times = [r["load_time_s"] for r in ok]
    clrs = []
    for r in ok:
        if r["data_type"] == "single_molecule":
            clrs.append(SM_COLOR)
        elif "MERSCOPE" in r["dataset_name"]:
            clrs.append(SC_MERSCOPE)
        else:
            clrs.append(SC_XENIUM)

    y = np.arange(len(ok))
    bars = ax1.barh(y, times, color=clrs, height=0.6)
    for bar, t in zip(bars, times):
        ax1.text(bar.get_width() + 1, bar.get_y() + bar.get_height()/2,
                 f"{t:.1f}s", va="center", fontsize=8, color="white")
    ax1.set_yticks(y)
    ax1.set_yticklabels(names, fontsize=7)
    ax1.set_xlabel("Load Time (s)")
    ax1.set_title("Load Time", fontsize=12, fontweight="bold")
    ax1.invert_yaxis()
    ax1.grid(axis="x", alpha=0.3)

    # ── Top right: Memory bars ──
    ax2 = fig.add_subplot(2, 2, 2)
    mem_ok = [r for r in ok if r["memory_peak_mb"] is not None]
    names_m = [r["dataset_name"].split("(")[0].strip()[:20] for r in mem_ok]
    mem = [r["memory_peak_mb"] for r in mem_ok]
    clrs_m = []
    for r in mem_ok:
        if r["data_type"] == "single_molecule":
            clrs_m.append(SM_COLOR)
        elif "MERSCOPE" in r["dataset_name"]:
            clrs_m.append(SC_MERSCOPE)
        else:
            clrs_m.append(SC_XENIUM)

    y2 = np.arange(len(mem_ok))
    bars2 = ax2.barh(y2, mem, color=clrs_m, height=0.6)
    for bar, m in zip(bars2, mem):
        ax2.text(bar.get_width() + 5, bar.get_y() + bar.get_height()/2,
                 f"{m:.0f} MB", va="center", fontsize=8, color="white")
    ax2.set_yticks(y2)
    ax2.set_yticklabels(names_m, fontsize=7)
    ax2.set_xlabel("Peak Memory (MB)")
    ax2.set_title("Memory Usage", fontsize=12, fontweight="bold")
    ax2.invert_yaxis()
    ax2.grid(axis="x", alpha=0.3)

    # ── Bottom left: FPS bars ──
    ax3 = fig.add_subplot(2, 2, 3)
    fps_ok = [r for r in ok if r["fps_after_load"] is not None]
    names_f = [r["dataset_name"].split("(")[0].strip()[:20] for r in fps_ok]
    fps = [r["fps_after_load"] for r in fps_ok]
    clrs_f = []
    for r in fps_ok:
        if r["data_type"] == "single_molecule":
            clrs_f.append(SM_COLOR)
        elif "MERSCOPE" in r["dataset_name"]:
            clrs_f.append(SC_MERSCOPE)
        else:
            clrs_f.append(SC_XENIUM)

    y3 = np.arange(len(fps_ok))
    bars3 = ax3.barh(y3, fps, color=clrs_f, height=0.6)
    ax3.axvline(x=30, color="#22c55e", linestyle="--", alpha=0.4)
    ax3.axvline(x=10, color="#eab308", linestyle="--", alpha=0.4)
    for bar, f in zip(bars3, fps):
        ax3.text(bar.get_width() + 0.5, bar.get_y() + bar.get_height()/2,
                 f"{f:.0f}", va="center", fontsize=8, color="white")
    ax3.set_yticks(y3)
    ax3.set_yticklabels(names_f, fontsize=7)
    ax3.set_xlabel("FPS")
    ax3.set_title("Rendering FPS", fontsize=12, fontweight="bold")
    ax3.invert_yaxis()
    ax3.grid(axis="x", alpha=0.3)

    # ── Bottom right: Status summary ──
    ax4 = fig.add_subplot(2, 2, 4)
    ax4.axis("off")

    total = len(results)
    succeeded = sum(1 for r in results if r["load_time_s"] is not None and not r["crashed"])
    crashed = sum(1 for r in results if r["crashed"])
    skipped = sum(1 for r in results if r["load_time_s"] is None and not r["crashed"])

    summary_text = (
        f"Total Datasets Tested: {total}\n\n"
        f"  Succeeded:  {succeeded}\n"
        f"  Crashed:    {crashed}  (OOM / timeout)\n"
        f"  Skipped:    {skipped}\n\n"
        f"Key Findings:\n"
        f"  - SM CSV: works up to ~32M molecules\n"
        f"    (67M+ times out after 93 min)\n"
        f"  - SC MERSCOPE: works up to ~150MB\n"
        f"    (800MB+ causes browser OOM)\n"
        f"  - SC Xenium: works up to 65GB!\n"
        f"    (only processes relevant files)\n"
        f"  - FPS degrades with cell count\n"
        f"    (43 FPS small -> 3 FPS very large)"
    )

    ax4.text(0.05, 0.95, summary_text, transform=ax4.transAxes,
             fontsize=10, verticalalignment="top", fontfamily="monospace",
             bbox=dict(boxstyle="round,pad=0.5", facecolor="#16213e",
                       edgecolor="#4a4a6a", alpha=0.9))

    # Legend at bottom
    legend_patches = [
        mpatches.Patch(facecolor=SM_COLOR, label="Single Molecule"),
        mpatches.Patch(facecolor=SC_MERSCOPE, label="SC MERSCOPE"),
        mpatches.Patch(facecolor=SC_XENIUM, label="SC Xenium"),
    ]
    fig.legend(handles=legend_patches, loc="lower center", ncol=3,
               framealpha=0.8, facecolor="#1a1a2e", edgecolor="#444",
               fontsize=10, bbox_to_anchor=(0.5, 0.01))

    fig.tight_layout(rect=[0, 0.04, 1, 0.95])
    fig.savefig(os.path.join(OUTPUT_DIR, "dashboard.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("  -> dashboard.png")


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    setup_style()
    results = load_data()

    print("Generating benchmark charts...")
    chart_load_times(results)
    chart_memory(results)
    chart_fps(results)
    chart_sm_scaling(results)
    chart_summary_dashboard(results)
    print(f"\nAll charts saved to {OUTPUT_DIR}/")


if __name__ == "__main__":
    main()
