#!/usr/bin/env python3
"""Neuraminidase phylogenies: does the post-2021 Egyptian H5N1 carry the resident Egyptian N1?

Two questions the HA tree cannot answer:
  A. N1. If the N1 of post-2021 clade 2.3.4.4b H5N1 falls inside the N1 diversity of the
     pre-2016 Egyptian clade 2.2.1.x H5N1 population, reassortment with the resident lineage
     is supported. If the two form separate groups, the N1 came from elsewhere, which is what
     a wild-bird re-introduction predicts.
  B. N2. The Discussion elsewhere invokes reassortment pressure from co-circulating H9N2.
     H9N2 carries N2, so that hypothesis is testable directly: do the Egyptian H5N2 NA
     sequences fall among Egyptian H9N2 N2, or apart from them?
"""
import csv, re, random, subprocess, os, json
from collections import Counter, defaultdict
import os
# Input directory holding the retrieved FASTA/metadata files.
# Override with:  export AIV_DATA=/path/to/data
DATA = os.environ.get("AIV_DATA", "data")

random.seed(20260929)
U = DATA.rstrip("/") + "/"
MAXG = 45

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

def yr(r):
    m = re.search(r'(19|20)\d{2}', str(r.get('collection_date') or '') + " " + str(r.get('strain') or ''))
    return int(m.group(0)) if m else None

meta = {}
for tag in ["H5", "H9"]:
    for r in csv.DictReader(open(f"{U}NCBI_{tag}_NA_metadata.tsv", encoding="utf-8-sig"), delimiter="\t"):
        meta[r['accession']] = r
seqs = {}
for tag in ["H5", "H9"]:
    for h, s in read_fa(f"{U}NCBI_{tag}_NA.fasta").items():
        seqs[h.split()[0]] = s

def pick(pred, label):
    out = []
    for acc, s in seqs.items():
        r = meta.get(acc)
        if not r: continue
        if 'egypt' not in (r.get('country', '') or '').lower(): continue
        y = yr(r)
        if not y or not (2010 <= y <= 2026): continue
        if len(s) < 1200: continue
        # the retrieval query's bare "NA" term also matches segment-4 records;
        # keep only true neuraminidase sequences
        if (r.get('segment','') or '').strip() != '6': continue
        if 'neuraminidase' not in (r.get('title','') or '').lower(): continue
        if pred(r, y): out.append((acc, y, r.get('serotype', ''), s))
    random.shuffle(out)
    sel = out[:MAXG]
    print(f"  {label}: {len(out)} available, {len(sel)} used")
    return sel

def build(groups, tag, model_nt=True):
    recs = []
    for gname, sel in groups:
        for acc, y, st, s in sel:
            recs.append((f"{acc.replace('.','-')}_{gname}_{st}_{y}", s))
    fa = f"work/na_{tag}.fna"
    with open(fa, "w") as f:
        for n, s in recs: f.write(f">{n}\n{s}\n")
    aln = f"work/na_{tag}_aln.fna"
    subprocess.run(f"mafft --auto --thread 4 --quiet {fa} > {aln}", shell=True, check=True)
    tre = f"work/na_{tag}.tree"
    subprocess.run(f"FastTree -nt -gtr -quiet {aln} > {tre}", shell=True, check=True)
    print(f"  {tag}: {len(recs)} sequences -> {tre}")
    return tre

os.makedirs("work", exist_ok=True)
print("Panel A - N1:")
gA = [("legacy", pick(lambda r, y: r['serotype'] == 'H5N1' and y <= 2015, "H5N1 NA, 2010-2015 (clade 2.2.1.x era)")),
      ("recent", pick(lambda r, y: r['serotype'] == 'H5N1' and y >= 2021, "H5N1 NA, 2021-2026 (clade 2.3.4.4b era)"))]
tA = build(gA, "N1")

print("Panel B - N2:")
gB = [("H5N2", pick(lambda r, y: r['serotype'] == 'H5N2', "H5N2 NA")),
      ("H9N2", pick(lambda r, y: r['serotype'] == 'H9N2', "H9N2 NA"))]
tB = build(gB, "N2")

# ---------------- topology test ----------------
from Bio import Phylo
res = {}
for tag, tre, groups in [("N1", tA, ["legacy", "recent"]), ("N2", tB, ["H5N2", "H9N2"])]:
    t = Phylo.read(tre, "newick")
    t.root_at_midpoint()
    tips = [x.name for x in t.get_terminals()]
    def grp(n): return n.split("_")[1]
    out = {"n_tips": len(tips), "groups": dict(Counter(grp(x) for x in tips))}
    for g in groups:
        members = [x for x in tips if grp(x) == g]
        if len(members) < 2: continue
        anc = t.common_ancestor(*[{"name": m} for m in members])
        desc = [x.name for x in anc.get_terminals()]
        purity = sum(1 for d in desc if grp(d) == g) / len(desc)
        out[g] = {"n": len(members), "mrca_tips": len(desc),
                  "mrca_purity": round(purity, 3),
                  "monophyletic": len(desc) == len(members),
                  "support": anc.confidence}
    res[tag] = out
    print(f"\n{tag}: {json.dumps(out, indent=1)}")
json.dump(res, open("out/na_phylo.json", "w"), indent=1)
