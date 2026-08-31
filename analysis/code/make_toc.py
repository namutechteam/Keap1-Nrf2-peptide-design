#!/usr/bin/env python3
"""ACS table-of-contents graphic (3.25 x 1.75 in, 600 dpi).

Usage: python make_toc.py <OUTDIR>
"""
import os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow, FancyBboxPatch

OUT = sys.argv[1] if len(sys.argv) > 1 else "."
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({"font.family": "sans-serif",
                     "font.sans-serif": ["Arial", "DejaVu Sans"]})

W, H = 3.25, 1.75
fig = plt.figure(figsize=(W, H))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

STEPS = [("ProteinMPNN", "530,299", "unique\nsequences", "#4C72B0"),
         ("sequence\nfilters", "5,118", "sequences", "#6A9AC4"),
         ("Chai-1 +\nHADDOCK3", "28", "peptides\nsynthesized", "#E8A33D"),
         ("FP assay", "1", "active", "#2E9E75")]

x0, w, gap = 3.0, 20.5, 2.2
for i, (tool, n, unit, col) in enumerate(STEPS):
    x = x0 + i * (w + gap)
    ax.add_patch(FancyBboxPatch((x, 46), w, 26, boxstyle="round,pad=0.4,rounding_size=2.0",
                                fc=col, ec="none"))
    ax.text(x + w / 2, 65.5, n, ha="center", va="center", fontsize=10.5,
            color="white", fontweight="bold")
    ax.text(x + w / 2, 53.5, unit, ha="center", va="center", fontsize=5.4,
            color="white", linespacing=1.15)
    ax.text(x + w / 2, 78.5, tool, ha="center", va="center", fontsize=6.2,
            color="#333333", linespacing=1.15)
    if i < len(STEPS) - 1:
        ax.add_patch(FancyArrow(x + w + 0.15, 59, gap - 0.5, 0, width=1.6,
                                head_width=5.2, head_length=1.6,
                                fc="#9A9A9A", ec="none", length_includes_head=True))

ax.text(50, 92, "A peptide design workflow built only from public tools,",
        ha="center", va="center", fontsize=6.6, color="#222222")
ax.text(50, 85.5, "with every prioritized design synthesized and measured",
        ha="center", va="center", fontsize=6.6, color="#222222")

# ---- outcome band -----------------------------------------------------------
ax.add_patch(FancyBboxPatch((3.0, 6), 94, 30, boxstyle="round,pad=0.4,rounding_size=2.0",
                            fc="#F4F4F2", ec="#DADAD6", lw=0.8))

def draw_runs(x, y, runs, size, family="monospace"):
    """Draw coloured text segments end to end, measuring each as it is placed."""
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    inv = ax.transData.inverted()
    for text, colour, weight in runs:
        t = ax.text(x, y, text, ha="left", va="center", fontsize=size,
                    family=family, color=colour, fontweight=weight)
        bb = t.get_window_extent(renderer=rend)
        x = inv.transform((bb.x1, bb.y0))[0]
    return x


ax.text(6.5, 27.5, "Hit", ha="left", va="center", fontsize=6.2, color="#666666")
end = draw_runs(14.0, 27.5,
                [("GLRLD", "#333333", "normal"),
                 ("PENGE", "#2E9E75", "bold"),
                 ("WN", "#333333", "normal")], 8.6)
ax.text(end + 6.0, 27.5, "96%", ha="left", va="center", fontsize=9.5,
        color="#2E9E75", fontweight="bold")
ax.text(end + 17.5, 26.9, "inhibition at 1 μM", ha="left", va="center",
        fontsize=5.6, color="#666666")

ax.text(6.5, 13.0, "No intact ETGE motif. Interface confidence, computed but not used",
        ha="left", va="center", fontsize=5.9, color="#444444")
ax.text(6.5, 7.6, "for selection, ranked this design 1 of 28.",
        ha="left", va="center", fontsize=5.9, color="#444444")

for ext, dpi in (("png", 600), ("tif", 600)):
    fig.savefig(os.path.join(OUT, f"TOC_re_write_3.25x1.75in.{ext}"), dpi=dpi,
                pil_kwargs={"compression": "tiff_lzw"} if ext == "tif" else None)
plt.close(fig)
print("wrote TOC_re_write_3.25x1.75in.png / .tif  ->", os.path.abspath(OUT))
