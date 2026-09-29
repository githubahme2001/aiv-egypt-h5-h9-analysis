#!/usr/bin/env python3
"""Fig. 9 - does the post-2021 Egyptian H5N1 sit inside the global clade 2.3.4.4b radiation?

The Egypt-only HA tree cannot distinguish local descent from re-introduction, because both
hypotheses put Egyptian H5N8 and Egyptian post-2021 H5N1 in the same clade. Adding contextual
taxa from the flyway countries makes the question answerable:

  * if the Egyptian post-2021 sequences form an exclusive group that is sister to the Egyptian
    2016-2020 H5N8, local descent is supported;
  * if they fall among non-Egyptian sequences, separated from the Egyptian H5N8, a separate
    introduction is supported, which is what the neuraminidase result (Fig. 8A) already implies.
"""
import csv, re, json, random, subprocess, os
from collections import Counter, defaultdict
import statistics as st
import os
# Input directory holding the retrieved FASTA/metadata files.
# Override with:  export AIV_DATA=/path/to/data
DATA = os.environ.get("AIV_DATA", "data")

random.seed(20260929)
U = DATA.rstrip("/") + "/"
NC = "/tmp/nc"; DS = "/tmp/h5ds"

def read_fa(fn):
    d = {}; n = None; b = []
    for l in open(fn, encoding="utf-8", errors="replace"):
        l = l.rstrip("\n")
        if l.startswith(">"):
            if n: d[n] = "".join(b)
            n = l[1:].strip(); b = []
        else: b.append(l.strip())
    if n: d[n] = "".join(b)
    return d

def year_of(v):
    m = re.search(r'(19|20)\d{2}', str(v))
    return int(m.group(0)) if m else None

# ---------------- contextual sequences ----------------
# Everything needed is in the NCBI defline, e.g.
#   >QB062156.1 Influenza A virus (A/Turkey/SD/26G11288-002-original/2026(H5N1)) segment 4 hemagglutinin (HA) gene
# so no separate metadata file is required.
ctx = read_fa(U + "NCBI_CONTEXT_H5_HA.fasta")
keep = {}
for h, s in ctx.items():
    acc = h.split()[0]
    hl = h.lower()
    if "segment 4" not in hl or "emagglutinin" not in hl: continue
    if len(s) < 1500: continue
    m = re.search(r'\((A/[^()]+?)\((H\d+N\d+|H\d+)\)\)', h) or re.search(r'\((A/[^()]+?)\)', h)
    if not m: continue
    strain = m.group(1)
    sero = m.group(2) if m.lastindex and m.lastindex >= 2 else ""
    y = year_of(strain.split("/")[-1]) or year_of(strain)
    if not y or y < 2015: continue
    # The retrieval term matched host names ("turkey vulture") and US state codes as well as
    # country names, so the contextual set is broader than the 26 countries requested. That is
    # acceptable - wider global context only strengthens the test - but the locality string is
    # not reliable, so contextual tips are labelled by region class rather than by country.
    parts = [p.strip() for p in strain.split("/")]
    loc = " ".join(parts[1:3]).lower()
    NA_STATES = set("mn oh sd bc qc ca mi ia wi ne ks mo ut id wa or ny nj pa va nc sc ga fl tx ok ar ms al tn ky in il mt wy co nm az nv ak hi me nh vt ma ct ri de md wv".split())
    if any(k in loc for k in ("minnesota","dakota","ohio","iowa","wisconsin","california","michigan","indiana","missouri","kansas","nebraska","colorado","montana","washington","oregon","texas","canada","quebec","ontario","alberta","british columbia")) \
       or any(p.lower() in NA_STATES for p in parts[1:3]):
        region = "N-America"
    elif any(k in loc for k in ("nigeria","niger","ghana","cameroon","burkina","south africa","lesotho","senegal","mali","benin","togo")):
        region = "Africa"
    elif any(k in loc for k in ("israel","iraq","saudi","kuwait","iran","jordan","lebanon","syria")):
        region = "W-Asia"
    elif any(k in loc for k in ("netherlands","germany","united kingdom","england","scotland","france","italy","poland","hungary","romania","bulgaria","greece","belgium","denmark","sweden","norway","spain","ireland","czech","austria","switzerland")):
        region = "Europe"
    elif any(k in loc for k in ("russia","kazakhstan","ukraine","siberia","georgia","azerbaijan")):
        region = "Eurasia"
    else:
        region = "other"
    keep[acc] = (region, y, sero, s)
print(f"contextual HA retained: {len(keep)} of {len(ctx)}")
print("  by country:", Counter(v[0] for v in keep.values()).most_common(10))

# ---------------- clade-assign the contextual set ----------------
with open("work/ctx.fna", "w") as f:
    for acc, (c, y, sero, s) in keep.items():
        f.write(f">{acc}\n{s}\n")
subprocess.run(f"{NC} run --input-dataset {DS} --output-tsv work/ctx_nc.tsv --jobs 4 work/ctx.fna",
               shell=True, check=True, capture_output=True)
cc = {}
for r in csv.DictReader(open("work/ctx_nc.tsv"), delimiter="\t"):
    cc[r['seqName'].split()[0]] = r['clade'].strip()
print("  contextual clades:", Counter(cc.values()).most_common(8))

# ---------------- assemble the tree set ----------------
eg = read_fa("out/H5_consolidated.fasta")
egc = json.load(open("out/clades.json"))
G221 = {"2.2.1", "2.2.1.1", "2.2.1.1a", "2.2.1.2", "2.2"}

def nm(acc, grp, extra, y):
    return f"{acc.replace('_','-').replace('.','-')}_{grp}_{extra}_{y}"

recs = []
# Egyptian strata
strat = defaultdict(list)
for k, s in eg.items():
    acc, sero, y, src = k.split("|"); y = int(y)
    cl = egc.get(acc)
    if not cl: continue
    if cl == "2.3.4.4b":
        if y <= 2020: strat["EGY-2344b-H5N8-2016-20"].append((acc, sero, y, s))
        else:         strat["EGY-2344b-2021-26"].append((acc, sero, y, s))
    elif cl in G221:
        strat["EGY-221x"].append((acc, sero, y, s))
for g, v in strat.items():
    sel = v if len(v) <= 30 else random.sample(v, 30)
    for acc, sero, y, s in sel: recs.append((nm(acc, g, sero, y), s))
    print(f"  {g}: {len(v)} available, {len(sel)} used")
# contextual 2.3.4.4b
ctx244 = [(a, v) for a, v in keep.items() if cc.get(a) == "2.3.4.4b"]
bycountry = defaultdict(list)
for a, v in ctx244: bycountry[v[0]].append((a, v))
sel = []
for cn, lst in bycountry.items():
    sel.extend(lst if len(lst) <= 6 else random.sample(lst, 6))
if len(sel) > 70: sel = random.sample(sel, 70)
for a, (cn, y, sero, s) in sel:
    recs.append((nm(a, "GLOBAL-2344b", cn.replace(' ', '-'), y), s))
print(f"  GLOBAL-2344b: {len(ctx244)} available, {len(sel)} used from {len(bycountry)} countries")
# outgroup: an Egyptian EA-nonGsGD if present
for k, s in eg.items():
    acc = k.split("|")[0]
    if egc.get(acc) == "EA-nonGsGD":
        recs.append((nm(acc, "OUTGROUP", "EA-nonGsGD", k.split("|")[2]), s)); break

with open("work/glob.fna", "w") as f:
    for n, s in recs: f.write(f">{n}\n{s}\n")
print(f"tree set: {len(recs)} sequences")
subprocess.run("mafft --auto --thread 4 --quiet work/glob.fna > work/glob_aln.fna", shell=True, check=True)
subprocess.run("FastTree -nt -gtr -quiet work/glob_aln.fna > work/glob.tree", shell=True, check=True)

# ---------------- the test ----------------
from Bio import Phylo
t = Phylo.read("work/glob.tree", "newick")
tips = [x.name for x in t.get_terminals()]
og = [x for x in tips if "_OUTGROUP_" in x]
if og: t.root_with_outgroup({"name": og[0]})
else:   t.root_at_midpoint()
def g(n): return n.split("_")[1]
res = {}
for grp in ["EGY-2344b-2021-26", "EGY-2344b-H5N8-2016-20", "GLOBAL-2344b", "EGY-221x"]:
    mem = [x for x in tips if g(x) == grp]
    if len(mem) < 2: continue
    anc = t.common_ancestor(*[{"name": m} for m in mem])
    desc = [x.name for x in anc.get_terminals()]
    res[grp] = {"n": len(mem), "mrca_tips": len(desc),
                "purity": round(sum(1 for d in desc if g(d) == grp) / len(desc), 3),
                "monophyletic": len(desc) == len(mem),
                "support": anc.confidence,
                "mrca_composition": dict(Counter(g(d) for d in desc))}
# nearest neighbour of each Egyptian post-2021 tip
rec = [x for x in tips if g(x) == "EGY-2344b-2021-26"]
nn = Counter()
for r in rec:
    d = sorted(((t.distance(r, o), g(o)) for o in tips if o != r))[:1]
    if d: nn[d[0][1]] += 1
res["nearest_neighbour_of_post2021_egyptian"] = dict(nn)
json.dump(res, open("out/global_tree.json", "w"), indent=1)
print(json.dumps(res, indent=1))
