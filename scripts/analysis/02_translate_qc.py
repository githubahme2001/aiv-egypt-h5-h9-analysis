#!/usr/bin/env python3
"""Translate + QC (manuscript section 2.2 method), then emit protein FASTA for MAFFT."""
import re, sys
from collections import Counter
TAB={}
_b="TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
_m="TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
_e="TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"
_a="FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for b,m,e,a in zip(_b,_m,_e,_a): TAB[b+m+e]=a
def tr(nt):
    nt=nt.upper().replace("U","T")
    return "".join(TAB.get(nt[i:i+3],"X") for i in range(0,len(nt)-2,3))
def km(s,k=6): return {s[i:i+k] for i in range(len(s)-k+1)}
def rf(fn):
    n=None;b=[]
    for l in open(fn):
        l=l.rstrip("\n")
        if l.startswith(">"):
            if n is not None: yield n,"".join(b)
            n=l[1:].strip();b=[]
        else: b.append(l.strip())
    if n is not None: yield n,"".join(b)
def longest(nt):
    best=""
    for f in range(3):
        for seg in tr(nt[f:]).split("*"):
            if len(seg)>len(best): best=seg
    return best

import os
TAGS=[("H5",1600)]+([("H9",1500)] if os.path.exists("out/H9_consolidated.fasta") else [])
for tag,minlen in TAGS:
    recs=list(rf(f"out/{tag}_consolidated.fasta"))
    prots=[(longest(s),h) for h,s in recs if len(s)>minlen]
    L=Counter(len(p) for p,_ in prots)
    modal=L.most_common(1)[0][0]
    REF=next(p for p,_ in prots if len(p)==modal); RK=km(REF)
    out=[];drop=Counter()
    for h,s in recs:
        c=[]
        for f in range(3):
            seg=max(tr(s[f:]).split("*"),key=len)
            if len(seg)<50: continue
            k=km(seg); c.append((len(k&RK)/max(1,len(k)),len(seg),seg))
        if not c: drop["no_orf"]+=1; continue
        c.sort(key=lambda x:(-x[0],-x[1])); cont,ln,prot=c[0]
        if ln<=150: drop["product<=150aa"]+=1; continue
        if cont<0.30: drop["low_containment"]+=1; continue
        out.append((h,prot))
    print(f"{tag}: reference {modal} aa (modal of {len(prots)} full-length) | "
          f"retained {len(out)}/{len(recs)} | drops {dict(drop)}")
    with open(f"out/{tag}_protein.faa","w") as f:
        for h,p in out: f.write(f">{h}\n{p}\n")
