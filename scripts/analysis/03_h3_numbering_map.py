#!/usr/bin/env python3
"""Build a per-position H5 -> H3 numbering map by alignment to PDB 2FK0 chain A.

The previous build used a single scalar offset calibrated at the receptor-binding
site. That is wrong away from the calibration point: the H5/H3 offset changes across
indel regions, so positions in the N-terminal and 130-loop regions were mis-assigned.

2FK0 chain A is an H5 HA1 deposited in H3 numbering, so aligning our alignment
consensus to it yields the H3 number of every column directly, with no scalar offset
and no assumption that the offset is constant.
"""
import re, json
from collections import Counter
from Bio import Align
import os
# Input directory holding the retrieved FASTA/metadata files.
# Override with:  export AIV_DATA=/path/to/data
DATA = os.environ.get("AIV_DATA", "data")


AA = "ACDEFGHIKLMNPQRSTVWY"

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

# ---- 2FK0 chain A: ordered (H3 label, residue) ----
pdb = []
for tok in open("reference/2FK0_chainA_CA_H3numbered.txt").read().strip().split(";"):
    rs, rn, _ = tok.split(":")
    pdb.append((rs, rn))
pdb_seq = "".join(r for _, r in pdb)

# ---- our alignment consensus, gap-aware ----
a5 = rd("out/H5_aln.faa")
W = len(next(iter(a5.values()))); N = len(a5)
cons = []
for p in range(W):
    ct = Counter(s[p] for s in a5.values()); g = ct.get('-', 0) + ct.get('X', 0)
    if g / N >= 0.50:
        cons.append('-')
    else:
        c2 = Counter(s[p] for s in a5.values() if s[p] in AA)
        cons.append(c2.most_common(1)[0][0] if c2 else '-')
cons = "".join(cons)

fp = cons.find("GLFGAI")                      # HA2 fusion peptide -> end of HA1
m = re.search(r"DQICIGYHA", cons)             # mature HA1 N-terminus
nt = m.start()
cols = [c for c in range(nt, fp) if cons[c] != '-']
ha1 = "".join(cons[c] for c in cols)
print(f"our HA1 consensus: {len(ha1)} residues; 2FK0 chain A: {len(pdb_seq)} residues")

# ---- global alignment, our HA1 vs 2FK0 ----
al = Align.PairwiseAligner()
al.mode = "global"
al.open_gap_score = -11
al.extend_gap_score = -1
al.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")
aln = al.align(ha1, pdb_seq)[0]
print(f"alignment score {aln.score:.0f}")

# map our HA1 index -> H3 label
col2h3 = {}
idx_a = idx_b = 0
A, B = aln.aligned
for (a0, a1), (b0, b1) in zip(A, B):
    for k in range(a1 - a0):
        our_i = a0 + k          # index into ha1
        pdb_i = b0 + k          # index into pdb_seq
        col2h3[cols[our_i]] = pdb[pdb_i][0]

ident = sum(1 for c, lab in col2h3.items() if cons[c] == pdb[[i for i, (l, _) in enumerate(pdb) if l == lab][0]][1])
print(f"mapped {len(col2h3)} of {len(cols)} HA1 columns to an H3 label; identity at mapped sites {ident}/{len(col2h3)} ({100*ident/len(col2h3):.1f}%)")

# integer H3 positions only (insertion codes such as 19A are not H3 positions)
col2pos = {c: int(l) for c, l in col2h3.items() if l.isdigit()}
pos2col = {v: k for k, v in col2pos.items()}

# ---- validate against the deposited residues ----
PDBRES = {int(l): r for l, r in pdb if l.isdigit()}
checks = [50, 58, 128, 138, 140, 155, 158, 185, 190, 193, 196, 222, 226, 227, 228]
print()
print("H3   2FK0  our dominant   observed(>5%)          ")
ok = 0; tot = 0
report = {}
for p in checks:
    c = pos2col.get(p)
    if c is None: continue
    ct = Counter(s[c] for s in a5.values() if s[c] in AA)
    tot_n = sum(ct.values())
    dom = ct.most_common(1)[0][0]
    obs = sorted(a for a, n2 in ct.items() if n2 / tot_n > 0.05)
    hit = PDBRES[p] in obs
    tot += 1; ok += hit
    report[p] = {"pdb": PDBRES[p], "dominant": dom, "observed": obs, "match": hit}
    print(f"{p:4d}   {PDBRES[p]}      {dom}            {obs}  {'OK' if hit else 'MISMATCH'}")
print()
print(f"deposited residue observed in the Egyptian data at {ok}/{tot} checked positions")

# implied offset per region, to document that it is not constant
offs = sorted({p - (list(col2pos.keys()).index(pos2col[p]) + 1) for p in pos2col if p in (20, 60, 100, 140, 180, 220, 260, 300)})
json.dump({"n_mapped": len(col2pos), "ha1_len": len(cols),
           "validation": report, "validated_ok": ok, "validated_total": tot,
           "col2pos": {str(k): v for k, v in col2pos.items()}},
          open("out/h3map.json", "w"), indent=1)
print("wrote out/h3map.json")
