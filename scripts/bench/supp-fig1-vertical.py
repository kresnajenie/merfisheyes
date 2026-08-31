"""Supplementary figure: paired browser/server bars, all qualifying datasets."""
import json, re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

SC = "/tmp/claude-1002/-home-kjenie-merfisheyes/31c58ad3-3115-47ed-9650-60c98f62dd39/scratchpad"
ROWS = json.load(open(f"{SC}/benchdata.json"))["rows"]

BLUE, ORANGE, GREY = "#2a78d6", "#eb6834", "#8a887f"
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 12,
    "axes.titlesize": 15, "xtick.labelsize": 12, "ytick.labelsize": 12,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.6, "svg.fonttype": "none",
})

DUPES = {"single-molecule__xenium-Xenium_Prime_Mouse_Pup_FFPE-parquet",
         "single-molecule__xenium-Human_Melanoma_10gb-csv",
         "h5ad__retroB30_all_BRBB_Axon_Virus_AnnotBRBBTransfer"}

def nice_name(r):
    n = r["id"].split("__", 1)[1]
    m = re.match(r"MouseSpinalCord_Run(\d+)_Slide(\d+)", n)
    if m:
        return f"Spinal cord R{m.group(1)}S{m.group(2)}"
    m = re.search(r"Mouse_receptors_rep(\d)_(\d)", n)
    if m:
        suffix = " adata" if n.endswith("adata") else ""
        prefix = "Receptors" if suffix else "Mouse receptors"
        return f"{prefix} r{m.group(1)}.{m.group(2)}{suffix}"
    m = re.search(r"800m-ace-dry-dog-detected_transcripts(\d)", n)
    if m:
        return f"MERSCOPE 800MB #{m.group(1)}"
    RULES = [
        ("merscope-800m-ace-dry-dog", "MERSCOPE 800MB #1"),
        ("merscope-2g-ace-ear-nap-202204011403_2022040113532", "MERSCOPE 2GB"),
        ("merscope-5g-ace-dry-dog", "MERSCOPE 5GB"),
        ("merscope-ed_lein_ace-irk-sag-H22-30-001-CX-13-01-0", "ed lein H22"),
        ("merscope-ed_lein_ace-irk-sag", "ed lein"),
        ("xenium-Human_Lymph_Node", "Lymph Node"),
        ("xenium-Whole_mouse_pup_65gb", "Whole mouse pup"),
        ("xenium-Human_Renal_Carcinoma_50gb", "Renal Carcinoma"),
        ("xenium-Human_Melanoma_10gb", "Melanoma"),
        ("xenium-Human_Lung_Cancer_FFPE", "Lung Cancer FFPE"),
        ("WTA_Preview_FFPE_Breast_Cancer", "WTA Breast"),
        ("WTA_Preview_FFPE_Cervical_Cancer", "WTA Cervical"),
        ("Xenium_Prime_Mouse_Pup_FFPE", "Prime mouse pup"),
        ("Human_Ovarian_Cancer_150mb", "Ovarian Cancer"),
        ("Human_Uterine_Cancer_800mb", "Uterine Cancer"),
        ("Human_Colon_Cancer_1-2gb", "Colon Cancer"),
        ("hum_brain2_280mb", "Human brain 2"),
        ("140MB_cbm2_combined_neuropath", "Cerebellum CBM2"),
        ("250MB_CBM7_labeled", "Cerebellum CBM7"),
        ("cbm5_combined_neuropath", "Cerebellum CBM5"),
        ("retroB30_all_BRBB_Axon_Virus_AnnotBRBBTransfer", "retroB30"),
        ("year1s_labelled_edu", "Year-1 EdU"),
        ("scdata_2000gnMERFISH_DCBB_5_9_2024__set1_FULL_corrected__wit", "MERFISH set1 FULL"),
        ("scdata_2000gnMERFISH_DCBB_5_9_2024__set3", "MERFISH set3"),
        ("4GB_2kset2_labelled_clusters", "2K-gene set2"),
        ("180MB_HubcTMA_VMSC02511_region_R2", "HubcTMA R2"),
        ("190MB_HubcTMA_VMSC02511_region_R3", "HubcTMA R3"),
        ("270MB_HubcTMA_VMSC02511_region_R4", "HubcTMA R4"),
        ("290MB_HubcTMA_VMSC02511_region_R5", "HubcTMA R5"),
        ("340MB_HubcTMA_VMSC02511_region_R1", "HubcTMA R1"),
        ("430MB_MsBrain_VMSC02011_region_R1", "MsBrain R1"),
        ("440MB_mmc_MsBrain_VMSC02011_region_R1", "MsBrain R1 mmc"),
        ("570MB_HuBrain_VMSC16710_region_R1", "HuBrain R1"),
    ]
    for pat, out in RULES:
        if pat in n:
            return out
    return n.replace("-parquet", "").replace("-csv", "").replace("_", " ")[:26]

def size_tag(r):
    sh = r["shape"]
    if sh.get("molecules"):
        m = sh["molecules"]
        return f"({m/1e6:.1f}M mol)" if m < 1e7 else f"({round(m/1e6)}M mol)"
    return f"({round(sh['cells']/1000)}K)"

def fmt_t(sec, tier=None):
    lab = f"{sec:.1f}s" if sec < 100 else f"{sec/60:.1f}m"
    return lab + (f" @{tier} GB" if tier else "")

FAIL = {"oom": "browser: crashed", "timeout": "browser: froze",
        "failed": "browser: refused (>2 GB)", "aborted": "browser: aborted"}

KEEP_SPINAL = {"xenium__MouseSpinalCord_Run2_Slide1_Mus_M2_F2_ace-pie-bog",
               "xenium__MouseSpinalCord_Run2_Slide2_Mus_M2_F2_ace-pie-boo",
               "xenium__MouseSpinalCord_Run4_Slide2_Mus_M3_F3_ace-pie-bun"}
KEEP_RECEPTORS = {"merscope__Mouse_receptors_rep1_2_240mb",
                  "merscope__Mouse_receptors_rep2_3_230mb",
                  "merscope__Mouse_receptors_rep3_2_220mb"}

def build(fmt_filter):
    out = []
    for r in ROWS:
        if r["id"] in DUPES or not fmt_filter(r):
            continue
        if "MouseSpinalCord" in r["id"] and r["id"] not in KEEP_SPINAL:
            continue
        if "Mouse_receptors" in r["id"] and r["format"] == "merscope" and r["id"] not in KEEP_RECEPTORS:
            continue
        l, s = r.get("local"), r.get("server")
        e = (s or {}).get("escalation")
        l_ok = l and l["outcome"] == "ok" and (l.get("genes") or 0) > 1
        l_deg = l and l["outcome"] == "ok" and (l.get("genes") or 0) <= 1
        s_ok = s and s["outcome"] == "ok"
        e_ok = e and e.get("outcome") == "ok"
        if not (l_ok or s_ok or e_ok):
            continue
        sh = r["shape"]
        file_gb = (r.get("server") or {}).get("uploadGB") or r.get("relevantGB") or r["sizeGB"]
        if r["format"] == "single-molecule":
            ftype = ".parquet" if r["id"].endswith("-parquet") else ".csv"
        elif r["format"] == "h5ad":
            ftype = ".h5ad"
        else:
            ftype = r["format"] + " folder"
        if sh.get("molecules"):
            m = sh["molecules"]
            mol = f"{m/1e6:.1f}M" if m < 1e7 else f"{round(m/1e6)}M"
            label = f"{mol} molecules\n{ftype} \u00b7 in {file_gb:.2f} GB"
        else:
            dense = sh["cells"] * (sh.get("genes") or 0) * 4 / 1e9
            c = sh["cells"]
            kc = f"{c/1e6:.1f}M" if c >= 1e6 else f"{round(c/1000)}K"
            label = f"{kc}\u00d7{sh.get('genes'):,}\nin {file_gb:.2f} GB \u00b7 dense {dense:.2f} GB"
        entry = {
            "label": label,
            "bt": l["processMs"] / 1000 if l_ok else None,
            "bnote": "browser: matrix exceeds tab memory" if l_deg else (FAIL.get(l["outcome"], "browser: failed") if l else None),
            "st": (e["processMs"] / 1000 if e_ok and not s_ok else (s["processMs"] / 1000 if s_ok else None)),
            "tier": (32 if e_ok and not s_ok and e["tier"] == "32768" else 64 if e_ok and not s_ok else None),
            "scale": r["shape"].get("molecules") or (r["shape"]["cells"] * (r["shape"].get("genes") or 1)),
        }
        out.append(entry)
    # Sort by processing time (server when present, else browser) — which
    # tracks the dense in-memory size almost exactly.
    out.sort(key=lambda x: x["st"] if x["st"] is not None else (x["bt"] or 0))
    return out

PANELS = [
    ("a", "Xenium single cell", lambda r: r["format"] == "xenium"),
    ("b", "MERSCOPE single cell", lambda r: r["format"] == "merscope"),
    ("c", ".h5ad single cell", lambda r: r["format"] == "h5ad"),
    ("d", "single molecule", lambda r: r["format"] == "single-molecule"),
]
data = {ltr: build(f) for ltr, _, f in PANELS}
for ltr in data:
    print(ltr, len(data[ltr]), "rows")


# ---- one horizontal-bar figure per panel -------------------------------------
from matplotlib.lines import Line2D

for ltr, title, _ in PANELS:
    entries = data[ltr]
    n = len(entries)
    fig, ax = plt.subplots(figsize=(10.5, 0.72 * n + 1.9))
    fig.subplots_adjust(left=0.30, right=0.97, top=1 - 0.85 / (0.72 * n + 1.9), bottom=0.85 / (0.72 * n + 1.9))
    xmax = max([e["st"] or 0 for e in entries] + [e["bt"] or 0 for e in entries]) / 60
    for i, e in enumerate(entries):
        y = n - 1 - i  # fastest at top
        if e["bt"] is not None:
            ax.barh(y + 0.20, e["bt"] / 60, height=0.36, color=BLUE)
            ax.text(e["bt"] / 60 + xmax * 0.012, y + 0.20, fmt_t(e["bt"]),
                    va="center", fontsize=11, color=BLUE)
        elif e["bnote"]:
            ax.text(xmax * 0.012, y + 0.20, e["bnote"], va="center",
                    fontsize=10, color=GREY, style="italic")
        if e["st"] is not None:
            ax.barh(y - 0.20, e["st"] / 60, height=0.36, color=ORANGE)
            ax.text(e["st"] / 60 + xmax * 0.012, y - 0.20, fmt_t(e["st"], e["tier"]),
                    va="center", fontsize=11, color=ORANGE)
    ax.set_yticks([n - 1 - i for i in range(n)], [e["label"] for e in entries], fontsize=12)
    ax.set_ylim(-0.7, n - 0.3)
    ax.set_xlim(0, xmax * 1.24)
    ax.set_xlabel("processing time (min)", fontsize=13)
    ax.set_title(title, fontsize=15, pad=10)
    ax.tick_params(length=2)
    handles = [Line2D([], [], marker="s", ls="", color=BLUE, ms=6, label="browser (local)"),
               Line2D([], [], marker="s", ls="", color=ORANGE, ms=6, label="server")]
    ax.legend(handles=handles, frameon=False, loc="upper right", fontsize=12)
    for ext in ("png", "svg", "pdf"):
        fig.savefig(f"{SC}/supp_fig1_{ltr}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"saved supp_fig1_{ltr}.png ({n} datasets)")
