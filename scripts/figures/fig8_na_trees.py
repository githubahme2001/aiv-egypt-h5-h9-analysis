#!/usr/bin/env python3
"""Fig. 8 - neuraminidase phylogenies testing the origin of the N1 and N2 segments.

Panel A asks whether the N1 of post-2021 clade 2.3.4.4b H5N1 is the resident Egyptian N1.
Panel B asks whether the N2 of the 2024 H5N2 reassortants came from co-circulating H9N2.
"""
import json, itertools, random, statistics as st
from collections import Counter
from Bio import Phylo
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection

BLUE="#0072B2"; ORANGE="#E69F00"; GREEN="#009E73"; PINK="#CC79A7"; INK="#222222"
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":7,"text.color":INK})
random.seed(20260929)

def grp(n): return n.split("_")[1]

def layout(tree):
    """Assign x = root-to-node distance, y = tip order."""
    tips = tree.get_terminals()
    y = {}
    for i, t in enumerate(tips): y[t] = i
    def gety(cl):
        if cl in y: return y[cl]
        v = st.mean(gety(c) for c in cl.clades); y[cl] = v; return v
    gety(tree.root)
    x = {}
    def getx(cl, acc=0.0):
        x[cl] = acc
        for c in cl.clades: getx(c, acc + (c.branch_length or 0.0))
    getx(tree.root)
    return x, y, tips

def render(ax, tree_file, colmap, title, note, legend_items):
    t = Phylo.read(tree_file, "newick"); t.root_at_midpoint(); t.ladderize()
    x, y, tips = layout(t)
    segs = []
    for cl in t.find_clades(order="level"):
        for c in cl.clades:
            segs.append([(x[cl], y[c]), (x[c], y[c])])        # horizontal
        if cl.clades:
            ys = [y[c] for c in cl.clades]
            segs.append([(x[cl], min(ys)), (x[cl], max(ys))])  # vertical
    ax.add_collection(LineCollection(segs, colors="#9aa0a8", linewidths=0.55))
    for tp in tips:
        g = grp(tp.name)
        ax.plot(x[tp], y[tp], marker="o", ms=3.1, mfc=colmap.get(g, "#888"),
                mec="none", zorder=4)
    ax.set_xlim(-0.004, max(x.values()) * 1.06)
    ax.set_ylim(-1.5, len(tips) + 0.5)
    ax.set_title(title, loc="left", fontsize=9.5, fontweight="bold", pad=6)
    ax.set_xlabel("substitutions per site", fontsize=8)
    ax.set_yticks([]); ax.tick_params(labelsize=7.2)
    for sp in ("top", "right", "left"): ax.spines[sp].set_visible(False)
    ax.legend(handles=legend_items, loc="lower right", fontsize=7.3, frameon=False,
              handletextpad=0.3, labelspacing=0.32)
    ax.text(0.015, 0.985, note, transform=ax.transAxes, fontsize=7.1, va="top",
            color="#333", linespacing=1.5)
    return t, x, y, tips

fig, axes = plt.subplots(1, 2, figsize=(11.6, 6.6))

# ---------------- Panel A : N1 ----------------
legA = [Line2D([], [], marker='o', ls='', ms=5.5, mfc=BLUE, mec='none',
               label="H5N1, 2010–2015 (clade 2.2.1.x era)"),
        Line2D([], [], marker='o', ls='', ms=5.5, mfc=ORANGE, mec='none',
               label="H5N1, 2021–2026 (clade 2.3.4.4b era)")]
noteA = ("The two eras form reciprocally monophyletic groups — 45/45 tips\n"
         "at SH-like support 0.996 and 42/42 at 1.000. Mean pairwise\n"
         "nucleotide identity is 97.4% within the 2010–2015 group and\n"
         "98.5% within the 2021–2026 group, but only 86.8% between them.\n"
         "The post-2021 N1 is therefore not the resident Egyptian N1.")
tA, xA, yA, tipsA = render(axes[0], "work/na_N1.tree",
                           {"legacy": BLUE, "recent": ORANGE},
                           "A   Neuraminidase N1, Egyptian H5N1", noteA, legA)

# ---------------- Panel B : N2 ----------------
legB = [Line2D([], [], marker='o', ls='', ms=5.5, mfc=GREEN, mec='none',
               label="H9N2 neuraminidase"),
        Line2D([], [], marker='o', ls='', ms=5.5, mfc=PINK, mec='none',
               label="H5N2 neuraminidase")]
noteB = ("Five of the six Egyptian H5N2 neuraminidase sequences, all from\n"
         "the 2024 reassortants, fall within the Egyptian H9N2 N2 group;\n"
         "mean identity between the two sets is 92.6%. The sixth is the\n"
         "2016 EA-nonGsGD wild-bird H5N2 and falls outside. Co-circulating\n"
         "H9N2 is the plausible donor of the 2024 N2 segment.")
tB, xB, yB, tipsB = render(axes[1], "work/na_N2.tree",
                           {"H9N2": GREEN, "H5N2": PINK},
                           "B   Neuraminidase N2, Egyptian H5N2 versus H9N2", noteB, legB)

plt.tight_layout(w_pad=2.0)
plt.savefig("out/Fig8.png", dpi=400, bbox_inches="tight", facecolor="white")
plt.close()
print("wrote out/Fig8.png")
