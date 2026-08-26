"""Supplementary Fig 1 in the simple bar style: browser load time per format."""
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SC = "/tmp/claude-1002/-home-kjenie-merfisheyes/31c58ad3-3115-47ed-9650-60c98f62dd39/scratchpad"
ROWS = {r["id"]: r for r in json.load(open(f"{SC}/benchdata.json"))["rows"]}

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 7,
    "axes.titlesize": 7.5, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.6, "svg.fonttype": "none",
})

def t(id_):  # browser process time in s
    return ROWS[id_]["local"]["processMs"] / 1000

def cells(id_):
    return ROWS[id_]["shape"]["cells"]

def k(n):
    return f"{round(n/1000)}K"

PANELS = [
    ("a", "Xenium single-cell load time", "#5aade0", [
        (f"Spinal cord S8.2 ({k(cells('xenium__MouseSpinalCord_Run8_Slide2_Mus_M7_F7_ace-pie-cot'))})", t("xenium__MouseSpinalCord_Run8_Slide2_Mus_M7_F7_ace-pie-cot")),
        (f"Spinal cord S1.1 ({k(cells('xenium__MouseSpinalCord_Run1_Slide1_Mus_M1_F1_ace-pie-bin'))})", t("xenium__MouseSpinalCord_Run1_Slide1_Mus_M1_F1_ace-pie-bin")),
        (f"Spinal cord S4.2 ({k(cells('xenium__MouseSpinalCord_Run4_Slide2_Mus_M3_F3_ace-pie-bun'))})", t("xenium__MouseSpinalCord_Run4_Slide2_Mus_M3_F3_ace-pie-bun")),
        (f"Spinal cord S2.1 ({k(cells('xenium__MouseSpinalCord_Run2_Slide1_Mus_M2_F2_ace-pie-bog'))})", t("xenium__MouseSpinalCord_Run2_Slide1_Mus_M2_F2_ace-pie-bog")),
    ]),
    ("b", "MERSCOPE single-cell load time", "#e06c6c", [
        (f"Ovarian Cancer ({k(cells('merscope__Human_Ovarian_Cancer_150mb'))})", t("merscope__Human_Ovarian_Cancer_150mb")),
        (f"Mouse receptors r3.1 ({k(cells('merscope__Mouse_receptors_rep3_1_190mb'))})", t("merscope__Mouse_receptors_rep3_1_190mb")),
        (f"Mouse receptors r2.3 ({k(cells('merscope__Mouse_receptors_rep2_3_230mb'))})", t("merscope__Mouse_receptors_rep2_3_230mb")),
        (f"Mouse receptors r1.2 ({k(cells('merscope__Mouse_receptors_rep1_2_240mb'))})", t("merscope__Mouse_receptors_rep1_2_240mb")),
    ]),
    ("c", ".h5ad single-cell load time", "#8fbf7a", [
        (f"Cerebellum CBM2 ({k(cells('h5ad__140MB_cbm2_combined_neuropath'))})", t("h5ad__140MB_cbm2_combined_neuropath")),
        (f"Cerebellum CBM7 ({k(cells('h5ad__250MB_CBM7_labeled'))})", t("h5ad__250MB_CBM7_labeled")),
        (f"Cerebellum CBM5 ({k(cells('h5ad__cbm5_combined_neuropath'))})", t("h5ad__cbm5_combined_neuropath")),
        (f"MERFISH set3 ({k(cells('h5ad__scdata_2000gnMERFISH_DCBB_5_9_2024__set3'))})", t("h5ad__scdata_2000gnMERFISH_DCBB_5_9_2024__set3")),
        (f"Year-1 EdU labelled ({k(cells('h5ad__year1s_labelled_edu'))})", t("h5ad__year1s_labelled_edu")),
    ]),
    ("d", "single-molecule load time", "#f2a75c", [
        ("MERSCOPE 40 MB (0.5M mol)", t("single-molecule__merscope-800m-ace-dry-dog-detected_transcripts3-csv")),
        ("Xenium 0.4 GB (24M mol)", t("single-molecule__xenium-Human_Melanoma_10gb-parquet")),
        ("Xenium 0.3 GB (32M mol)", t("single-molecule__xenium-Human_Lung_Cancer_FFPE-parquet")),
    ]),
]

fig, axes = plt.subplots(1, 4, figsize=(11.5, 1.9))
fig.subplots_adjust(left=0.10, right=0.985, top=0.80, bottom=0.26, wspace=1.35)
for ax, (letter, title, color, bars) in zip(axes, PANELS):
    labels = [b[0] for b in bars][::-1]
    vals = [b[1] for b in bars][::-1]
    ax.barh(range(len(vals)), vals, color=color, edgecolor="none", height=0.62)
    ax.set_yticks(range(len(vals)), labels)
    for i, v in enumerate(vals):
        ax.text(v + max(vals) * 0.03, i, f"{v:.1f}s", va="center", fontsize=6)
    ax.set_xlim(0, max(vals) * 1.28)
    from matplotlib.ticker import MaxNLocator
    ax.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=5))
    ax.set_xlabel("load time (s)", fontsize=6.5)
    ax.set_title(title, fontsize=7, pad=4)
    ax.tick_params(length=2)
    ax.text(-0.02, 1.22, letter, transform=ax.transAxes, fontsize=10, fontweight="bold", ha="right")

for ext in ("svg", "pdf", "png"):
    fig.savefig(f"{SC}/supp_fig1_bars.{ext}", dpi=300, bbox_inches="tight")
print("saved supp_fig1_bars.{svg,pdf,png}")
