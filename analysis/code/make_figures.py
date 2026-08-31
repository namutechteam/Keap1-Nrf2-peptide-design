#!/usr/bin/env python3
"""Two figures for the Keap1-Nrf2 manuscript.

Figure2_cascade_composition  Figure 4: composition of the surviving population
                             across the design cascade.
Figure3_computed_not_used    Quantities the workflow computed but did not use
                             for selection, against measured inhibition.

Usage: python make_figures.py <ANALYSIS_DIR> <DATASET_CSV> <OUTDIR>
"""
import csv, json, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ANA, DATASET, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7.5,
    "axes.linewidth": 0.7, "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "xtick.major.size": 3, "ytick.major.size": 3, "legend.frameon": False,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})

E79, T80, E82 = "#0072B2", "#999999", "#D55E00"      # colour-blind safe
HIT, DEAD, REF = "#009E73", "#BBBBBB", "#CC79A7"


def panel(ax, letter, dx=-0.17, dy=1.06):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontweight="bold",
            fontsize=9, va="top", ha="left")


def tidy(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def save(fig, name):
    for ext, dpi in (("png", 300), ("tif", 600)):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"), dpi=dpi,
                    pil_kwargs={"compression": "tiff_lzw"} if ext == "tif" else None)
    plt.close(fig)
    print("  wrote", name + ".png /", name + ".tif")


# ---------------------------------------------------------------- Figure 4 --
comp = list(csv.DictReader(open(os.path.join(ANA, "table_composition_by_stage.csv"))))
brk = list(csv.DictReader(open(os.path.join(ANA, "table_etge_break_position.csv"))))
hadd = json.load(open(os.path.join(ANA, "si_table_S4_haddock.json")))

SHORT = ["Generated", "Unique", "Filters 1–5", "Filter 6", "Filter 7", "Selected"]
NLAB = [f"{int(r['n']):,}" for r in comp]

fig = plt.figure(figsize=(7.0, 2.5))
gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1.0, 0.85], wspace=0.42)

# (a) native-residue retention across the cascade
ax = fig.add_subplot(gs[0])
x = range(len(comp))
for key, col, lab, mk in (("E79_percent", E79, "Glu79", "o"),
                          ("T80_percent", T80, "Thr80", "s"),
                          ("E82_percent", E82, "Glu82", "^")):
    ax.plot(x, [float(r[key]) for r in comp], marker=mk, ms=4, lw=1.4,
            color=col, label=lab, clip_on=False)
ax.axvspan(3.5, 4.5, color="#F0E442", alpha=0.30, lw=0)
ax.text(4, 111, "ETGE\nexclusion", ha="center", va="bottom", fontsize=6.5, color="#7a6a00")
ax.set_xticks(list(x))
ax.set_xticklabels([a + "\n" + b for a, b in zip(SHORT, NLAB)], rotation=34,
                   ha="right", fontsize=6.6, linespacing=1.3)
ax.set_ylim(0, 110); ax.set_yticks([0, 25, 50, 75, 100])
ax.set_ylabel("Sequences retaining the\nnative residue (%)")
ax.legend(loc="center", bbox_to_anchor=(0.63, 0.34), handlelength=1.4, fontsize=7)
tidy(ax); panel(ax, "a", dx=-0.30)

# (b) where the ETGE pattern was broken
ax = fig.add_subplot(gs[1])
keys = ["E79", "T80", "E82"]
lbl = ["Glu79", "Thr80", "Glu82"]
share = [float(next(r for r in brk if r["deviating_position(s)"] == k)["percent_of_5118"])
         for k in keys]
avail = [8.6, 2.2, 1.8]                      # non-native frequency entering filter 7
w = 0.36
xs = range(3)
ax.bar([i - w / 2 for i in xs], share, w, color=[E79, T80, E82], label="_")
ax.bar([i + w / 2 for i in xs], avail, w, color="white", edgecolor="#444444",
       hatch="////", lw=0.7)
for i, (s, a) in enumerate(zip(share, avail)):
    ax.text(i - w / 2, s + 1.6, f"{s:.1f}", ha="center", fontsize=6.6)
    ax.text(i + w / 2, a + 1.6, f"{a:.1f}", ha="center", fontsize=6.6)
ax.set_xticks(list(xs)); ax.set_xticklabels(lbl)
ax.set_ylim(0, 100); ax.set_ylabel("Percent")
ax.set_xlabel("Position at which ETGE was broken")
ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, fc="#777777", ec="none"),
                   plt.Rectangle((0, 0), 1, 1, fc="white", ec="#444444", hatch="////", lw=0.7)],
          labels=["Share of the 5,118\nretained sequences",
                  "Non-native frequency\nentering filter 7"],
          loc="upper right", handlelength=1.1, labelspacing=0.7, borderpad=0.2, fontsize=6.4)
tidy(ax); panel(ax, "b", dx=-0.26)

# (c) the structure-based steps, 10-15-residue classes
ax = fig.add_subplot(gs[2])
order = ["Sequence filter 7 (10–15 aa)", "HADDOCK3 checkpoint (10–15 aa)",
         "Structure-based selection (10–15 aa)"]
rows = {(r["stage"], r["position"]): r for r in hadd}
xs = range(3)
for pos, col, lab, mk in ((79, E79, "Glu79", "o"), (82, E82, "Glu82", "^")):
    ax.plot(xs, [rows[(s, pos)]["E" if pos in (79, 82) else "T"] for s in order],
            marker=mk, ms=4.5, lw=1.6, color=col, label=lab, clip_on=False)
for pos, col, dy in ((79, E79, 9), (82, E82, 9)):
    for i, s in enumerate(order):
        ax.annotate(f"{rows[(s, pos)]['E']:.1f}", (i, rows[(s, pos)]["E"]),
                    textcoords="offset points", xytext=(0, dy), ha="center",
                    fontsize=6.4, color=col)
ax.set_xticks(list(xs))
ax.set_xticklabels(["Filter 7\n(4,347)", "HADDOCK3\n(115)", "Selected\n(26)"], fontsize=7)
ax.set_ylim(-4, 120); ax.set_yticks([0, 25, 50, 75, 100])
ax.set_ylabel("Sequences retaining the\nnative residue (%)")
ax.legend(loc="center left", bbox_to_anchor=(0.06, 0.55), handlelength=1.4, fontsize=7)
tidy(ax); panel(ax, "c", dx=-0.34)

save(fig, "Figure2_cascade_composition")

# ------------------------------------------------- computed vs measured --
ds = [r for r in csv.DictReader(open(DATASET)) if r["round"] == "first-round"]
met = list(csv.DictReader(open(os.path.join(ANA, "table_metric_activity.csv"))))

fig = plt.figure(figsize=(7.0, 2.7))
gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.15], wspace=0.34)

# (a) interface confidence against measured inhibition
ax = fig.add_subplot(gs[0])
mark = {"Comp11": (HIT, "Comp11"), "Comp25": (REF, "Comp25"), "Comp27": ("#56B4E9", "Comp27")}
for r in ds:
    xi, yi = float(r["chai_iptm_mean"]), float(r["inhibition_1uM_pct"])
    if r["peptide_id"] in mark:
        c, _ = mark[r["peptide_id"]]
        ax.scatter(xi, yi, s=42, color=c, edgecolor="black", lw=0.6, zorder=3)
    else:
        ax.scatter(xi, yi, s=20, color=DEAD, edgecolor="none", zorder=2)
ann = {"Comp11": (-46, -4), "Comp25": (6, 9), "Comp27": (6, 9)}
for r in ds:
    if r["peptide_id"] in mark:
        ax.annotate(f"{r['peptide_id']}\n(ipTM {float(r['chai_iptm_mean']):.2f})",
                    (float(r["chai_iptm_mean"]), float(r["inhibition_1uM_pct"])),
                    textcoords="offset points", xytext=ann[r["peptide_id"]],
                    fontsize=6.6, color=mark[r["peptide_id"]][0], fontweight="bold")
ax.axhline(50, ls=":", lw=0.8, color="#888888")
ax.text(0.355, 52, "active (>50%)", fontsize=6.4, color="#888888")
ax.set_xlabel("Chai-1 interface pTM (mean of five models)")
ax.set_ylabel("Inhibition at 1 μM (%)")
ax.set_xlim(0.33, 0.92); ax.set_ylim(-6, 108)
tidy(ax); panel(ax, "a", dx=-0.20)

# (b) rank of the single active design under every computed quantity
ax = fig.add_subplot(gs[1])
def rank_val(r):
    return float(r["rank_of_active"].split("-")[0])
met_sorted = sorted(met, key=rank_val)
ys = range(len(met_sorted))
cols = ["#0072B2" if r["metric"].startswith("Chai-1") else "#E69F00" for r in met_sorted]
ax.barh(list(ys), [rank_val(r) for r in met_sorted], color=cols, height=0.68)
NICE = {"E_vdw": r"$E_\mathrm{vdw}$", "E_elec": r"$E_\mathrm{elec}$",
        "E_desolv": r"$E_\mathrm{desolv}$"}
ax.set_yticks(list(ys))
ax.set_yticklabels([NICE.get(r["metric"], r["metric"].replace("Chai-1 ", ""))
                    for r in met_sorted], fontsize=6.8)
ax.invert_yaxis()
ax.axvline(14.5, ls="--", lw=0.9, color="#555555")
ax.text(14.2, -0.9, "expected rank if uninformative", fontsize=6.4,
        color="#555555", va="bottom", ha="right")
for i, r in enumerate(met_sorted):
    ax.text(rank_val(r) + 0.4, i, r["rank_of_active"], va="center", fontsize=6.4)
ax.set_xlim(0, 29); ax.set_ylim(len(met_sorted) - 0.4, -1.6)
ax.set_xlabel("Rank of Comp11 among the 28 designs")
ax.legend(handles=[Line2D([], [], color="#0072B2", lw=5),
                   Line2D([], [], color="#E69F00", lw=5)],
          labels=["Chai-1 confidence", "HADDOCK3 metric"],
          loc="upper right", bbox_to_anchor=(1.0, 0.99), handlelength=1.0, fontsize=7)
tidy(ax); panel(ax, "b", dx=-0.40)

save(fig, "Figure3_computed_not_used")
print("figures written to", os.path.abspath(OUT))
