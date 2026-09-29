#!/usr/bin/env python3
"""Re-test the post-2021 antigenic-site entropy decline with sampling depth controlled.

Three analyses:
  (1) the original year-level Wilcoxon, reproduced, plus the exact test;
  (2) rarefaction of every year to a common depth, bootstrapped;
  (3) the full breakpoint scan, so the 2021 choice is disclosed rather than implied.
"""
import json, math, re, random
from collections import Counter, defaultdict
import statistics as st
import numpy as np
from scipy.stats import ranksums, mannwhitneyu

AA = "ACDEFGHIKLMNPQRSTVWY"
random.seed(20260926); np.random.seed(20260926)

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

a5 = rd(__import__("os").environ.get("ALN","out/H5_aln.faa"))

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
_m=_j.load(open(__import__("os").environ.get("MAP","out/h3map.json")))["col2pos"]
c2p5={int(k):int(v) for k,v in _m.items()}; O5=0

SITES = {"A": list(range(122, 147)),
         "B": list(range(155, 164)) + list(range(187, 199)),
         "C": list(range(50, 58)) + list(range(275, 280)),
         "D": list(range(201, 221)),
         "E": list(range(62, 84))}
allsite = set().union(*SITES.values())
ag_cols = [c for c in sorted(c2p5) if (c2p5[c] + O5) in allsite]

by_year = defaultdict(list)
for k, s in a5.items():
    _, _, yr, _ = k.split("|")
    by_year[int(yr)].append(s)

def H(chars):
    ct = Counter(c for c in chars if c in AA); n = sum(ct.values())
    if n == 0: return None
    return -sum((v/n) * math.log2(v/n) for v in ct.values())

def year_entropy(seqs):
    vals = [H([s[c] for s in seqs]) for c in ag_cols]
    vals = [v for v in vals if v is not None]
    return st.mean(vals) if vals else None

years = sorted(y for y in by_year if len(by_year[y]) >= 10)
obs = {y: year_entropy(by_year[y]) for y in years}
depth = {y: len(by_year[y]) for y in years}

out = {"n_ag_cols": len(ag_cols), "years": years,
       "depth": depth, "observed_entropy": {y: round(obs[y], 4) for y in years}}

# ---------------------------------------------------------- (1) as published
def split(bp):
    return [obs[y] for y in years if y < bp], [obs[y] for y in years if y >= bp]

pre, post = split(2021)
out["published"] = {
    "pre_n": len(pre), "pre_mean": round(st.mean(pre), 4),
    "post_n": len(post), "post_mean": round(st.mean(post), 4),
    "ranksums_p": round(float(ranksums(pre, post).pvalue), 4),
    "exact_mwu_p": round(float(mannwhitneyu(pre, post, alternative="two-sided",
                                            method="exact").pvalue), 4)}

# ---------------------------------------------------------- (2) rarefaction
D = min(depth.values())
out["rarefaction_depth"] = D
B = 500
rare = {y: [] for y in years}
for y in years:
    seqs = by_year[y]
    for _ in range(B):
        sub = random.sample(seqs, D)
        rare[y].append(year_entropy(sub))
rare_mean = {y: st.mean(rare[y]) for y in years}
out["rarefied_entropy"] = {y: round(rare_mean[y], 4) for y in years}

# breakpoint test on rarefied point estimates
rpre = [rare_mean[y] for y in years if y < 2021]
rpost = [rare_mean[y] for y in years if y >= 2021]
out["rarefied_2021"] = {
    "pre_mean": round(st.mean(rpre), 4), "post_mean": round(st.mean(rpost), 4),
    "exact_mwu_p": round(float(mannwhitneyu(rpre, rpost, alternative="two-sided",
                                            method="exact").pvalue), 4)}

# bootstrap the test itself: one rarefied draw per year per replicate
ps = []
for b in range(B):
    pr = [rare[y][b] for y in years if y < 2021]
    po = [rare[y][b] for y in years if y >= 2021]
    ps.append(float(mannwhitneyu(pr, po, alternative="two-sided", method="exact").pvalue))
ps = np.array(ps)
out["rarefied_bootstrap"] = {
    "median_p": round(float(np.median(ps)), 4),
    "frac_p_below_0.05": round(float((ps < 0.05).mean()), 3),
    "p_2.5pct": round(float(np.percentile(ps, 2.5)), 4),
    "p_97.5pct": round(float(np.percentile(ps, 97.5)), 4)}

# Hodges-Lehmann shift + bootstrap CI on the rarefied difference
diffs = sorted(a - b for a in rpre for b in rpost)
hl = float(np.median(diffs))
bs = []
for b in range(B):
    pr = [rare[y][b] for y in years if y < 2021]
    po = [rare[y][b] for y in years if y >= 2021]
    bs.append(float(np.median(sorted(x - y2 for x in pr for y2 in po))))
out["rarefied_effect"] = {"hodges_lehmann_bits": round(hl, 4),
                          "ci95": [round(float(np.percentile(bs, 2.5)), 4),
                                   round(float(np.percentile(bs, 97.5)), 4)]}

# ---------------------------------------------------------- (3) breakpoint scan
scan = {}
for bp in years[2:-2]:
    a, b2 = split(bp)
    if len(a) >= 3 and len(b2) >= 3:
        scan[bp] = round(float(mannwhitneyu(a, b2, alternative="two-sided",
                                            method="exact").pvalue), 4)
out["breakpoint_scan_observed"] = scan
if scan:
    best = min(scan, key=scan.get)
    out["scan_min_year"] = best
    out["scan_min_p"] = scan[best]
    out["scan_n_tested"] = len(scan)
    out["scan_bonferroni_p"] = round(min(1.0, scan[best] * len(scan)), 4)

# same scan on rarefied estimates
rscan = {}
for bp in years[2:-2]:
    a = [rare_mean[y] for y in years if y < bp]
    b2 = [rare_mean[y] for y in years if y >= bp]
    if len(a) >= 3 and len(b2) >= 3:
        rscan[bp] = round(float(mannwhitneyu(a, b2, alternative="two-sided",
                                             method="exact").pvalue), 4)
out["breakpoint_scan_rarefied"] = rscan

json.dump(out, open(__import__("os").environ.get("OUT","out/entropy_rarefied.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
