"""Supplementary Figure 1: browser vs server processing benchmarks."""
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

SC = "/tmp/claude-1002/-home-kjenie-merfisheyes/31c58ad3-3115-47ed-9650-60c98f62dd39/scratchpad"
DATA = json.load(open(f"{SC}/benchdata.json"))["rows"]

BLUE, ORANGE, GREY = "#2a78d6", "#eb6834", "#8a887f"
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 7, "axes.labelsize": 7.5, "axes.titlesize": 8,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "svg.fonttype": "none",
})

sc_rows = [r for r in DATA if r["format"] != "single-molecule" and r["shape"].get("cells") and r["shape"].get("genes")]
sm_rows = [r for r in DATA if r["format"] == "single-molecule" and r["shape"].get("molecules")]
x_sc = lambda r: r["shape"]["cells"] * r["shape"]["genes"]
x_sm = lambda r: r["shape"]["molecules"]
deg = lambda m: m and m.get("outcome") == "ok" and (m.get("genes") is not None and m["genes"] <= 1)

fig = plt.figure(figsize=(7.1, 6.4))
gs = fig.add_gridspec(2, 2, hspace=0.48, wspace=0.38,
                      left=0.10, right=0.965, top=0.885, bottom=0.07)

def split_axes(cell):
    sub = cell.subgridspec(1, 2, wspace=0.08)
    a1 = fig.add_subplot(sub[0]); a2 = fig.add_subplot(sub[1])
    a2.tick_params(labelleft=False)
    return a1, a2

def scatter_panel(ax, rows, xf, value, fail_val, ylim):
    """One sub-axis: browser + server dots, hollow=degraded, x=failure."""
    for r in rows:
        x = xf(r)
        l, s, e = r.get("local"), r.get("server"), (r.get("server") or {}).get("escalation")
        lv, lf = value(l, "local"), fail_val(l, "local")
        if lv is not None:
            if deg(l):
                ax.scatter(x, lv, s=11, facecolors="none", edgecolors=BLUE, linewidths=0.9, zorder=3)
            else:
                ax.scatter(x, lv, s=11, color=BLUE, zorder=3)
        elif lf is not None:
            ax.scatter(x, lf, s=16, color=BLUE, marker="x", linewidths=0.9, zorder=3)
        sv, sf = value(s, "server"), fail_val(s, "server")
        ev = value_esc(e) if e else None
        if sv is not None:
            (ax.scatter(x, sv, s=11, facecolors="none", edgecolors=ORANGE, linewidths=0.9, zorder=3)
             if deg(s) else ax.scatter(x, sv, s=11, color=ORANGE, zorder=3))
        elif ev is not None:
            ax.scatter(x, ev, s=16, color=ORANGE, marker="^", zorder=3)
        elif sf is not None:
            ax.scatter(x, sf, s=16, color=ORANGE, marker="x", linewidths=0.9, zorder=3)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylim(*ylim)

# ---- panel a: processing time --------------------------------------------
ca, cb = gs[0, 0], gs[0, 1]
a1, a2 = split_axes(ca)
TL = (0.8, 9000)
value_esc = lambda e: e["processMs"] / 1000 if e.get("outcome") == "ok" and e.get("processMs") else None
tval = lambda m, side: (m["processMs"] / 1000) if m and m.get("outcome") == "ok" and m.get("processMs") else None
tfail = lambda m, side: TL[1] * 0.75 if m and m.get("outcome") not in (None, "ok") else None
scatter_panel(a1, sc_rows, x_sc, tval, tfail, TL)
scatter_panel(a2, sm_rows, x_sm, tval, tfail, TL)
a1.text(0.03, 0.955, "\u2715 did not process", transform=a1.transAxes, ha="left", va="top", fontsize=5.8, color=GREY)
a1.set_xlabel("cells × genes"); a2.set_xlabel("molecules")
a1.set_ylabel("processing time (s)")
a1.set_title("single cell", loc="left", fontsize=7, color=GREY)
a2.set_title("single molecule", loc="left", fontsize=7, color=GREY)
fig.text(0.012, 0.90, "a", fontsize=11, fontweight="bold")

# ---- panel b: peak memory -------------------------------------------------
b1, b2 = split_axes(cb)
ML = (0.1, 90)
value_esc = lambda e: e.get("peakRssGB") if e.get("outcome") == "ok" else None
mval = lambda m, side: m.get("peakRssGB") if m and m.get("outcome") == "ok" else None
def mfail(m, side):
    if not m or m.get("outcome") in (None, "ok"):
        return None
    return m.get("peakRssGB") if m.get("outcome") == "oom" and m.get("peakRssGB") else None
for ax in ():
    pass
scatter_panel(b1, sc_rows, x_sc, mval, mfail, ML)
scatter_panel(b2, sm_rows, x_sm, mval, mfail, ML)
for ax in (b1, b2):
    ax.axhline(3.5, ls="--", lw=0.7, color=BLUE, alpha=0.55)
    ax.axhline(16, ls="--", lw=0.7, color=ORANGE, alpha=0.55)
b2.text(0.97, 3.5, "3.5 GB tab heap", transform=b2.get_yaxis_transform(), fontsize=5.8, color=BLUE, va="top", ha="right")
b2.text(0.97, 16, "16 GB server tier", transform=b2.get_yaxis_transform(), fontsize=5.8, color=ORANGE, va="bottom", ha="right")
b1.set_xlabel("cells × genes"); b2.set_xlabel("molecules")
b1.set_ylabel("peak memory (GB)")
b1.set_title("single cell", loc="left", fontsize=7, color=GREY)
b2.set_title("single molecule", loc="left", fontsize=7, color=GREY)
fig.text(0.525, 0.90, "b", fontsize=11, fontweight="bold")

# ---- panel c: escalation ladder ------------------------------------------
cx = fig.add_subplot(gs[1, 0])
esc_rows = [r for r in DATA if (r.get("server") or {}).get("escalation")]
esc_rows.sort(key=lambda r: (r["server"]["escalation"].get("peakRssGB") or 99))
labels, ys = [], []
for i, r in enumerate(esc_rows):
    e = r["server"]["escalation"]
    name = r["id"].split("__", 1)[1]
    name = (name.replace("WTA_Preview_FFPE_", "WTA ").replace("xenium-", "")
                .replace("Human_", "").replace("_65gb", "").replace("_Cancer", "").replace("Whole_mouse_pup", "Mouse_pup").replace("_", " "))
    kind = " (molecules)" if r["format"] == "single-molecule" else " (cells)"
    name = name.replace("-parquet", "") + kind
    labels.append(name); ys.append(i)
    if e["outcome"] == "ok" and e.get("peakRssGB"):
        cx.scatter(e["peakRssGB"], i, s=22, color=ORANGE, zorder=3)
        cx.text(e["peakRssGB"] - 2.4, i, f'{e["peakRssGB"]:.1f} GB', ha="right", va="center", fontsize=6.2, color=ORANGE)
    else:
        cx.scatter(62.5, i, s=22, color=GREY, marker="x", linewidths=1.0, zorder=3)
        cx.text(60.5, i, "exceeds 64 GB", ha="right", va="center", fontsize=6.2, color=GREY)
for tier in (16, 32, 64):
    cx.axvline(tier, ls="--", lw=0.7, color=GREY, alpha=0.6)
    cx.text(tier, len(esc_rows) - 0.3, f"{tier} GB", ha="center", fontsize=6, color=GREY)
cx.set_yticks(ys, labels)
cx.set_xlim(0, 66); cx.set_ylim(-0.6, len(esc_rows) - 0.2)
cx.set_xlabel("peak memory (GB)")
cx.invert_yaxis()
box = cx.get_position()
cx.set_position([0.20, box.y0, box.x1 - 0.20, box.height])
fig.text(0.012, 0.455, "c", fontsize=11, fontweight="bold")

# ---- panel d: cost = billed time x tier ----------------------------------
dx = fig.add_subplot(gs[1, 1])
TIER_RATE = {16: 4 * 0.04048 + 16 * 0.004445, 32: 8 * 0.04048 + 32 * 0.004445, 64: 16 * 0.04048 + 64 * 0.004445}
TIER_COLOR = {16: ORANGE, 32: "#b34700", 64: "#6e2c00"}
for r in DATA:
    c = (r.get("server") or {}).get("cost")
    if not c:
        continue
    for run in c["runs"]:
        if not run.get("runtimeMs") or run["usd"] <= 0:
            continue
        tier = int(run["memGB"])
        dx.scatter(run["runtimeMs"] / 1000, run["usd"], s=10,
                   color=TIER_COLOR.get(tier, GREY),
                   marker="o" if run["status"] == "SUCCEEDED" else "x", linewidths=0.8, zorder=3)
t = np.logspace(1.3, 3.9, 50)
for tier, rate in TIER_RATE.items():
    dx.plot(t, t * rate / 3600, lw=0.7, color=TIER_COLOR[tier], alpha=0.55)
    dx.text(t[-1], t[-1] * rate / 3600, f" {tier} GB", fontsize=6, color=TIER_COLOR[tier], va="center")
dx.set_xscale("log"); dx.set_yscale("log")
dx.set_xlabel("billed runtime (s)"); dx.set_ylabel("compute cost (USD)")
fig.text(0.525, 0.455, "d", fontsize=11, fontweight="bold")

# ---- legend ---------------------------------------------------------------
from matplotlib.lines import Line2D
handles = [
    Line2D([], [], marker="o", ls="", color=BLUE, ms=4, label="browser, complete"),
    Line2D([], [], marker="o", ls="", mfc="none", mec=BLUE, ms=4, label="browser, expression skipped"),
    Line2D([], [], marker="x", ls="", color=BLUE, ms=4.5, label="browser, failed"),
    Line2D([], [], marker="o", ls="", color=ORANGE, ms=4, label="server (16 GB)"),
    Line2D([], [], marker="^", ls="", color=ORANGE, ms=4.5, label="server, escalated tier"),
    Line2D([], [], marker="x", ls="", color=ORANGE, ms=4.5, label="server, out of memory"),
]
fig.legend(handles=handles, loc="upper center", ncol=3, frameon=False,
           bbox_to_anchor=(0.5, 1.005), columnspacing=1.6, handletextpad=0.4)

for ext in ("svg", "pdf", "png"):
    fig.savefig(f"{SC}/supp_fig1.{ext}", dpi=300)
print("saved supp_fig1.{svg,pdf,png}")
