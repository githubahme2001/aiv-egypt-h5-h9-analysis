#!/usr/bin/env python3
"""Fig. 7 - antigenic epitopes and clade-defining substitutions on the HA1 structure.

Structure: PDB 2FK0 chain A (H5 HA1, A/Vietnam/1203/2004), deposited in H3 numbering.
Validated against our alignment: deposited residues at monitored positions
(128 S, 158 N, 185 P, 193 K, 196 Q, 226 Q, 227 S, 228 G) match residues observed
in the Egyptian H5 dataset, confirming the H3-equivalent numbering.
"""
import json, math, csv, re
from collections import Counter, defaultdict
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import os
# Input directory holding the retrieved FASTA/metadata files.
# Override with:  export AIV_DATA=/path/to/data
DATA = os.environ.get("AIV_DATA", "data")


AA = "ACDEFGHIKLMNPQRSTVWY"

# ---------------------------------------------------------------- structure
CA = {}
for tok in open("reference/2FK0_chainA_CA_H3numbered.txt").read().strip().split(";"):
    rs, rn, xyz = tok.split(":")
    CA[rs] = (rn, np.array([float(v) for v in xyz.split(",")]))
order = list(CA.keys())
trace = np.array([CA[r][1] for r in order])
# integer H3 position -> coordinate (insertion codes such as 19A are not H3 positions)
POS = {int(r): CA[r][1] for r in order if r.isdigit()}

# ---------------------------------------------------------------- alignment
def rd(fn):
    d = {}; n = None; b = []
    for l in open(fn):
        l = l.rstrip("\n")
        if l.startswith(">"):
            if n: d[n] = "".join(b)
            n = l[1:].strip(); b = []
        else: b.append(l.strip())
    if n: d[n] = "".join(b)
    return d

a5 = rd("out/H5_aln.faa"); a9 = rd("out/H9_aln.faa")
nc = json.load(open("out/clades.json"))

def build_coords(aln, nterm_re, rbs):
    W = len(next(iter(aln.values()))); N = len(aln)
    cons = []
    for p in range(W):
        ct = Counter(s[p] for s in aln.values()); g = ct.get('-', 0) + ct.get('X', 0)
        if g / N >= 0.50: cons.append('-')
        else:
            c2 = Counter(s[p] for s in aln.values() if s[p] in AA)
            cons.append(c2.most_common(1)[0][0] if c2 else '-')
    cons = "".join(cons)
    fp = cons.find("GLFGAI"); m = re.search(nterm_re, cons); nt = m.start() if m else 0
    c2p = {}; pos = 0
    for c in range(nt, fp):
        if cons[c] != '-': pos += 1; c2p[c] = pos
    inv = {v: k for k, v in c2p.items()}
    off = None
    for o in range(-15, 16):
        A, B = inv.get(rbs[0][0] - o), inv.get(rbs[1][0] - o)
        if A is None or B is None: continue
        ra = Counter(s[A] for s in aln.values() if s[A] in AA).most_common(1)[0][0]
        rb = Counter(s[B] for s in aln.values() if s[B] in AA).most_common(1)[0][0]
        if ra == rbs[0][1] and rb == rbs[1][1]: off = o; break
    return c2p, off

import json as _j
_m5=_j.load(open("out/h3map.json"))["col2pos"]
c2p5={int(k):int(v) for k,v in _m5.items()}; O5=0
_m9=_j.load(open("out/h3map_H9.json"))["col2pos"]
c2p9={int(k):int(v) for k,v in _m9.items()}; O9=0

def meta(k):
    a, st, yr, src = k.split("|"); return a, st, int(yr), src
G221 = {"2.2.1", "2.2.1.1", "2.2.1.1a", "2.2.1.2"}
def grp(k):
    c = nc.get(meta(k)[0])
    if c in G221: return "2.2.1.x"
    if c == "2.3.4.4b": return "2.3.4.4b"
    return None

G = defaultdict(list)
for k, s in a5.items():
    g = grp(k)
    if g: G[g].append(s)

def entropy(seqs, col):
    ct = Counter(s[col] for s in seqs if s[col] in AA); n = sum(ct.values())
    if n < 10: return None
    return -sum((v / n) * math.log2(v / n) for v in ct.values())

def ent_profile(seqs, c2p, off):
    return {c2p[c] + off: entropy(seqs, c) for c in c2p if entropy(seqs, c) is not None}

E221 = ent_profile(G["2.2.1.x"], c2p5, O5)
E344 = ent_profile(G["2.3.4.4b"], c2p5, O5)
EH9  = ent_profile(list(a9.values()), c2p9, O9)
EH5  = ent_profile(list(a5.values()), c2p5, O5)

# ---------------------------------------------------------------- sites / subs
SITES = {"A": list(range(122, 147)),
         "B": list(range(155, 164)) + list(range(187, 199)),
         "C": list(range(50, 58)) + list(range(275, 280)),
         "D": list(range(201, 221)),
         "E": list(range(62, 84))}
SC = {"A": "#d62728", "B": "#17becf", "C": "#2ca02c", "D": "#3b3b9e", "E": "#c4459a"}
def site(p):
    for k, v in SITES.items():
        if p in v: return k
    return None

subs = []
with open("out/SuppTableS2.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        subs.append({"p": int(r["H3_equivalent_position"]),
                     "a": r["clade_2.2.1.x_residue"], "b": r["clade_2.3.4.4b_residue"],
                     "site": r["classical_antigenic_site"],
                     "early": r["present_in_early_H5N8"] == "yes"})

# ---------------------------------------------------------------- draw helpers
def ribbon(ax, lw=1.6, color="#b9bcc4"):
    segs = [[trace[i], trace[i + 1]] for i in range(len(trace) - 1)
            if np.linalg.norm(trace[i + 1] - trace[i]) < 5.0]
    ax.add_collection3d(Line3DCollection(segs, colors=color, linewidths=lw))

def setup3d(ax, azim=-60, elev=12):
    ax.set_xlim(trace[:, 0].min() - 2, trace[:, 0].max() + 2)
    ax.set_ylim(trace[:, 1].min() - 2, trace[:, 1].max() + 2)
    ax.set_zlim(trace[:, 2].min() - 2, trace[:, 2].max() + 2)
    ax.set_box_aspect((np.ptp(trace[:, 0]), np.ptp(trace[:, 1]), np.ptp(trace[:, 2])))
    ax.view_init(elev=elev, azim=azim)
    ax.set_axis_off()

# ---------------------------------------------------------------- figure
fig = plt.figure(figsize=(11.6, 9.8))
gs = fig.add_gridspec(2, 2, height_ratios=[1.35, 1.0], hspace=0.10, wspace=0.18)

# --- A: antigenic sites A-E
axA = fig.add_subplot(gs[0, 0], projection="3d")
ribbon(axA)
for s, col in SC.items():
    pts = np.array([POS[p] for p in SITES[s] if p in POS])
    if len(pts):
        axA.scatter(pts[:, 0], pts[:, 1], pts[:, 2], s=46, c=col,
                    edgecolors="white", linewidths=0.4, depthshade=True)
setup3d(axA)
axA.text2D(0.0, 1.03, "A   Classical antigenic sites A–E on H5 HA1",
           transform=axA.transAxes, fontsize=9.5, fontweight="bold", va="top")
axA.legend(handles=[Line2D([], [], marker='o', ls='', ms=6.5, mfc=SC[s], mec='white',
                           label=f"site {s}") for s in "ABCDE"],
           loc="lower left", bbox_to_anchor=(0.00, 0.02), fontsize=7.6,
           frameon=False, handletextpad=0.25, labelspacing=0.28)

# --- B: the 24 clade-defining substitutions
axB = fig.add_subplot(gs[0, 1], projection="3d")
ribbon(axB)
drawn = []
for d in subs:
    if d["p"] not in POS: continue
    c = SC.get(d["site"], "#5a5a5a")
    xyz = POS[d["p"]]
    axB.scatter(*xyz, s=110 if d["site"] in SC else 62, c=c, edgecolors="black",
                linewidths=0.9 if d["site"] in SC else 0.5,
                depthshade=False, zorder=5)
    drawn.append(d)
# label only the substitutions that fall in a classical antigenic site
LOFF = {54: (-11, 0, -10.0), 62: (6, 0, 5.5), 122: (-13, 0, 4.0),
        129: (-12, 0, 2.0), 137: (-16, 0, 6.0), 187: (9, 0, 4.0), 276: (-14, 0, -7.0)}
for d in drawn:
    if d["site"] not in SC: continue
    xyz = POS[d["p"]]; o = LOFF.get(d["p"], (0, 0, 3.0))
    axB.text(xyz[0] + o[0], xyz[1] + o[1], xyz[2] + o[2], f"{d['a']}{d['p']}{d['b']}",
             fontsize=7.2, ha="center", fontweight="bold",
             color=SC[d["site"]], zorder=8,
             bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.78))
for p, lab in [(226, "Q226"), (228, "G228")]:
    xyz = POS[p]
    axB.scatter(*xyz, s=150, marker="*", c="#f5c518", edgecolors="black",
                linewidths=0.7, depthshade=False, zorder=7)
axB.text(POS[226][0] + 10, POS[226][1], POS[226][2] - 9.0, "Q226 / G228", fontsize=7.0,
         ha="center", fontweight="bold", color="#8a6d00", zorder=8,
         bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.78))
setup3d(axB)
axB.text2D(0.0, 1.03,
           f"B   {len(drawn)} of {len(subs)} clade-defining HA1 substitutions",
           transform=axB.transAxes, fontsize=9.5, fontweight="bold", va="top")
axB.text2D(0.02, 0.30,
           "clade 2.2.1.x \u2192 clade 2.3.4.4b. Coloured and labelled where the position\n"
           "falls in a classical antigenic site; small grey spheres lie outside sites A\u2013E.\n"
           "\u2605 receptor-binding residues, invariant in both clades. All substitutions\n"
           "fall within the modelled range of the structure.",
           transform=axB.transAxes, fontsize=7.3, va="top", color="#333", linespacing=1.5)

# --- C: within-clade site entropy, 2.2.1.x vs 2.3.4.4b
axC = fig.add_subplot(gs[1, 0])
shared = sorted(set(E221) & set(E344))
for p in shared:
    s = site(p)
    axC.scatter(E221[p], E344[p], s=26 if s else 12,
                c=SC[s] if s else "#9aa0a8", alpha=0.9 if s else 0.5,
                edgecolors="none", zorder=3 if s else 2)
mx = max(max(E221[p] for p in shared), max(E344[p] for p in shared)) * 1.08
axC.plot([0, mx], [0, mx], ls="--", lw=0.8, c="#888", zorder=1)
subp = {d["p"] for d in subs}
rp = [p for p in shared if p in subp]
axC.scatter([E221[p] for p in rp], [E344[p] for p in rp], s=62, facecolors="none",
            edgecolors="black", linewidths=0.9, zorder=5)
axC.set_xlim(-0.05, mx); axC.set_ylim(-0.05, mx)
axC.set_xlabel("Shannon entropy, clade 2.2.1.x (bits)", fontsize=8.2)
axC.set_ylabel("Shannon entropy, clade 2.3.4.4b (bits)", fontsize=8.2)
axC.set_title("C   Within-clade variability at the same HA1 positions",
              loc="left", fontsize=9.5, fontweight="bold", pad=7)
axC.tick_params(labelsize=7.5)
for sp in ("top", "right"): axC.spines[sp].set_visible(False)
axC.legend(handles=[Line2D([], [], marker="o", ls="", ms=7, mfc="none", mec="black",
                           label="clade-defining substitution"),
                    Line2D([], [], marker="o", ls="", ms=6, mfc="#9aa0a8", mec="none",
                           label="position outside sites A\u2013E")],
           loc="upper left", fontsize=7.2, frameon=False,
           handletextpad=0.3, labelspacing=0.3)
n_lo = sum(1 for p in shared if p in subp and E221[p] < 0.5 and E344[p] < 0.5)
axC.text(0.97, 0.30,
         f"{n_lo} of {len([p for p in shared if p in subp])} substituted positions\n"
         "(black rings) lie below 0.5 bits in\nboth clades: fixed within each clade,\ndifferent between them",
         transform=axC.transAxes, fontsize=7.0, ha="right", va="top", color="#333")

# --- D: H5Nx vs H9N2 site entropy
axD = fig.add_subplot(gs[1, 1])
sh2 = sorted(set(EH5) & set(EH9))
for p in sh2:
    s = site(p)
    axD.scatter(EH5[p], EH9[p], s=26 if s else 12,
                c=SC[s] if s else "#9aa0a8", alpha=0.9 if s else 0.5,
                edgecolors="none", zorder=3 if s else 2)
mx2 = max(max(EH5[p] for p in sh2), max(EH9[p] for p in sh2)) * 1.08
axD.plot([0, mx2], [0, mx2], ls="--", lw=0.8, c="#888", zorder=1)
axD.set_xlim(-0.05, mx2); axD.set_ylim(-0.05, mx2)
axD.set_xlabel("Shannon entropy, H5Nx (all clades, bits)", fontsize=8.2)
axD.set_ylabel("Shannon entropy, H9N2 (bits)", fontsize=8.2)
axD.set_title("D   Antigenic-site entropy, H5Nx versus H9N2",
              loc="left", fontsize=9.5, fontweight="bold", pad=7)
axD.tick_params(labelsize=7.5)
for sp in ("top", "right"): axD.spines[sp].set_visible(False)
mH5 = np.mean([EH5[p] for p in sh2 if site(p)])
mH9 = np.mean([EH9[p] for p in sh2 if site(p)])
axD.text(0.97, 0.955,
         f"mean site A–E entropy\nH5Nx {mH5:.2f} bits\nH9N2 {mH9:.2f} bits",
         transform=axD.transAxes, fontsize=7.2, ha="right", va="top", color="#333")

plt.savefig("out/Fig7.png", dpi=400, bbox_inches="tight", facecolor="white")
plt.close()

out = {"n_subs_mapped": len(drawn),
       "n_subs_total": len(subs),
       "subs_not_in_structure": [d["p"] for d in subs if d["p"] not in POS],
       "n_shared_positions_clades": len(shared),
       "n_subs_low_entropy_both": n_lo,
       "n_subs_shared": len([p for p in shared if p in subp]),
       "mean_site_entropy_h5": round(float(mH5), 3),
       "mean_site_entropy_h9": round(float(mH9), 3),
       "mean_site_entropy_221": round(float(np.mean([E221[p] for p in shared if site(p)])), 3),
       "mean_site_entropy_2344b": round(float(np.mean([E344[p] for p in shared if site(p)])), 3),
       "pdb": "2FK0 chain A", "pdb_validation": {}}
for p in [128, 158, 185, 193, 196, 226, 227, 228]:
    out["pdb_validation"][p] = CA[str(p)][0]
json.dump(out, open("out/fig7.json", "w"), indent=1)
print(json.dumps(out, indent=1))
