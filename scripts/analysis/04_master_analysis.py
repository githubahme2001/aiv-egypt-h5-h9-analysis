#!/usr/bin/env python3
"""FINAL master analysis — Egyptian sequences only, formal Nextclade clades."""
import json, math, csv, re
from collections import Counter, defaultdict
import statistics as stat
from scipy.stats import ranksums
AA="ACDEFGHIKLMNPQRSTVWY"
def rd(fn):
    d={};n=None;b=[]
    for l in open(fn):
        l=l.rstrip("\n")
        if l.startswith(">"):
            if n: d[n]="".join(b)
            n=l[1:].strip();b=[]
        else: b.append(l.strip())
    if n: d[n]="".join(b)
    return d
import os
a5=rd("out/H5_aln.faa")
HAS9=os.path.exists("out/H9_aln.faa")
a9=rd("out/H9_aln.faa") if HAS9 else {}
nc=json.load(open("out/clades.json"))
def build_coords(aln,nterm_re,rbs):
    W=len(next(iter(aln.values()))); N=len(aln)
    cons=[]
    for p in range(W):
        ct=Counter(s[p] for s in aln.values()); g=ct.get('-',0)+ct.get('X',0)
        if g/N>=0.50: cons.append('-')
        else:
            c2=Counter(s[p] for s in aln.values() if s[p] in AA)
            cons.append(c2.most_common(1)[0][0] if c2 else '-')
    cons="".join(cons)
    fp=cons.find("GLFGAI"); m=re.search(nterm_re,cons); nt=m.start() if m else 0
    c2p={};pos=0
    for c in range(nt,fp):
        if cons[c]!='-': pos+=1; c2p[c]=pos
    inv={v:k for k,v in c2p.items()}
    off=None
    for o in range(-15,16):
        A,B=inv.get(rbs[0][0]-o),inv.get(rbs[1][0]-o)
        if A is None or B is None: continue
        ra=Counter(s[A] for s in aln.values() if s[A] in AA).most_common(1)[0][0]
        rb=Counter(s[B] for s in aln.values() if s[B] in AA).most_common(1)[0][0]
        if ra==rbs[0][1] and rb==rbs[1][1]: off=o; break
    return c2p,inv,off,pos
c2p5,inv5,O5,L5=build_coords(a5,r"DQICIGYHA",[(226,"Q"),(228,"G")])
# --- override the scalar offset with the alignment-based H5 -> H3 map (work/62_h3map.py).
# A single offset is wrong away from its calibration point because the H5/H3 offset
# changes across indel regions; the map is validated against PDB 2FK0 at 15/15 positions.
import json as _json
_m=_json.load(open("out/h3map.json"))["col2pos"]
c2p5={int(k):int(v) for k,v in _m.items()}
inv5={v:k for k,v in c2p5.items()}
O5=0
L5=max(c2p5.values())
c2p9,inv9,O9,L9=(build_coords(a9,r"[DN][KQ]ICIG",[(226,"L"),(228,"G")]) if HAS9 else ({},{},0,0))
if HAS9:
    _m9=_json.load(open("out/h3map_H9.json"))["col2pos"]
    c2p9={int(k):int(v) for k,v in _m9.items()}
    inv9={v:k for k,v in c2p9.items()}
    O9=0
    L9=max(c2p9.values())
print(f"H5 HA1={L5} offset=+{O5} | H9 HA1={L9} offset=+{O9}")
def meta(k):
    a,st,yr,src=k.split("|"); return a,st,int(yr),src
G221={"2.2.1","2.2.1.1","2.2.1.1a","2.2.1.2"}
def grp(k):
    c=nc.get(meta(k)[0])
    if c in G221: return "2.2.1.x"
    if c=="2.3.4.4b": return "2.3.4.4b"
    return "EA-nonGsGD" if c=="EA-nonGsGD" else None
G=defaultdict(list)
for k,s in a5.items():
    g=grp(k)
    if g: G[g].append((k,s))
R={}
R["dataset"]={"h5":len(a5),"h9":len(a9),"total":len(a5)+len(a9),
 "h5n1":sum(1 for k in a5 if meta(k)[1]=="H5N1"),"h5n8":sum(1 for k in a5 if meta(k)[1]=="H5N8"),
 "h5n2":sum(1 for k in a5 if meta(k)[1]=="H5N2"),
 "by_clade":dict(Counter(nc.get(meta(k)[0]) for k in a5)),
 "by_group":{g:len(v) for g,v in G.items()},
 "h5_ha1":L5,"h9_ha1":L9,"h5_off":O5,"h9_off":O9,
 "h5_src":dict(Counter(meta(k)[3] for k in a5)),"h9_src":dict(Counter(meta(k)[3] for k in a9))}
t=defaultdict(Counter)
for k in a5:
    c=nc.get(meta(k)[0])
    if c: t[meta(k)[2]][c]+=1
R["clade_by_year"]={str(y):dict(t[y]) for y in sorted(t)}
R["subtype_clade"]={st:dict(Counter(nc.get(meta(k)[0]) for k in a5 if meta(k)[1]==st)) for st in ["H5N1","H5N8","H5N2"]}
ny=defaultdict(Counter)
for k in a5:
    if grp(k)=="2.3.4.4b": ny[meta(k)[2]][meta(k)[1]]+=1
R["na_within_2344b"]={str(y):dict(ny[y]) for y in sorted(ny)}
def dom(g,c):
    ct=Counter(s[c] for _,s in g if s[c] in AA); n=sum(ct.values())
    if not n: return None,0
    a,x=ct.most_common(1)[0]; return a,x/n
SITES={"A":range(122,147),"B":list(range(155,164))+list(range(187,199)),
       "C":list(range(50,58))+list(range(275,280)),"D":range(201,221),"E":range(62,84)}
def site(p):
    for k,v in SITES.items():
        if p in v: return k
    return None
def subs(x,y,thr):
    o=[]
    for c in sorted(c2p5):
        a,fa=dom(x,c); b,fb=dom(y,c)
        if a and b and a!=b and fa>=thr and fb>=thr: o.append((c2p5[c]+O5,c,a,b,round(fa,3),round(fb,3)))
    return o
n8e=[(k,s) for k,s in a5.items() if meta(k)[1]=="H5N8" and meta(k)[2]<=2019]
s90=subs(G["2.2.1.x"],G["2.3.4.4b"],0.90); s50=subs(G["2.2.1.x"],G["2.3.4.4b"],0.50)
ok=sum(1 for p,c,a,b,_,_ in s90 if dom(n8e,c)[0]==b)
R["substitutions"]={"n_majority":len(s50),"n_090":len(s90),
 "positions":[p for p,*_ in s90],"list":[f"{p}{a}>{b}" for p,c,a,b,_,_ in s90],
 "at_sites":[[p,f"{a}>{b}",site(p)] for p,c,a,b,_,_ in s90 if site(p)],
 "n_at_sites":sum(1 for p,c,a,b,_,_ in s90 if site(p)),
 "site_breakdown":dict(Counter(site(p) for p,c,a,b,_,_ in s90 if site(p))),
 "already_in_early_H5N8":f"{ok}/{len(s90)}"}
with open("out/SuppTableS2.tsv","w",newline="") as f:
    w=csv.writer(f,delimiter="\t")
    w.writerow(["H3_equivalent_position","alignment_column","clade_2.2.1.x_residue","freq_2.2.1.x",
        "clade_2.3.4.4b_residue","freq_2.3.4.4b","classical_antigenic_site",
        "residue_in_H5N8_2016_2019","freq_H5N8_2016_2019","present_in_early_H5N8"])
    for p,c,a,b,fa,fb in s90:
        a8,f8=dom(n8e,c)
        w.writerow([p,c+1,a,fa,b,fb,site(p) or "-",a8 or "-",round(f8,3) if a8 else "-","yes" if a8==b else "no"])
def cons(g,cols):
    f=[];s=0;n=0
    for c in cols:
        ct=Counter(x[c] for _,x in g if x[c] in AA); t2=sum(ct.values())
        if t2<10: continue
        n+=1; fr=ct.most_common(1)[0][1]/t2; f.append(fr); s+= fr>=0.90
    return [round(stat.mean(f),4),round(s/n,3),n]
cols5=sorted(c2p5); cols9=sorted(c2p9)
R["conservation"]={g:cons(G[g],cols5) for g in ["2.2.1.x","2.3.4.4b"]}
R["conservation"]["h5_pooled"]=cons(G["2.2.1.x"]+G["2.3.4.4b"],cols5)
R["conservation"]["h9"]=cons(list(a9.items()),cols9) if HAS9 else None
def rbs(g,p,inv,off):
    c=inv[p-off]; ct=Counter(s[c] for _,s in g if s[c] in AA)
    a,n=ct.most_common(1)[0]; return [a,n,sum(ct.values()),round(n/sum(ct.values()),3)]
R["rbs"]={f"h5_{p}_{g}":rbs(G[g],p,inv5,O5) for p in (226,228) for g in ["2.2.1.x","2.3.4.4b"]}
R["rbs"]["h9_226"]=rbs(list(a9.items()),226,inv9,O9) if HAS9 else None
R["rbs"]["h9_228"]=rbs(list(a9.items()),228,inv9,O9) if HAS9 else None
P5=[83,128,138,140,156,158,160,163,165,167,182,185,193,195,196,197,204,207,208,217,226,227,228,263]
P9=[83,95,127,134,148,155,156,158,160,182,193,195,198,226,227,228]
def panel(aln,P,inv,off,minseq=10):
    by=defaultdict(list)
    for k,s in aln.items(): by[meta(k)[2]].append(s)
    yrs=[y for y in sorted(by) if len(by[y])>=minseq]
    stable=0; changed=[]; tot=0
    for p in P:
        c=inv.get(p-off)
        if c is None: continue
        tot+=1; ds=[]
        for y in yrs:
            ct=Counter(s[c] for s in by[y] if s[c] in AA); n=sum(ct.values())
            if n<minseq: continue
            a,k2=ct.most_common(1)[0]; ds.append((y,a,k2/n))
        res={d[1] for d in ds}
        if len(res)==1 and all(d[2]>=0.90 for d in ds): stable+=1
        else: changed.append([p,sorted(res),round(min(d[2] for d in ds),2),
                              [d[0] for i,d in enumerate(ds) if i and d[1]!=ds[i-1][1]]])
    return stable,tot,changed
s5,t5,c5c=panel(a5,P5,inv5,O5)
s9,t9,c9c=(panel(a9,P9,inv9,O9) if HAS9 else (0,0,[]))
R["panel"]={"h5_stable":s5,"h5_total":t5,"h5_changed":c5c,
            "h5_residue_changes":[c for c in c5c if len(c[1])>1],
            "h9_stable":s9,"h9_total":t9,"h9_changed":c9c,
            "h9_residue_changes":[c for c in c9c if len(c[1])>1]}
ag5=[c for c in cols5 if site(c2p5[c]+O5)]
ag9=[inv9[p-O9] for p in P9 if (p-O9) in inv9] if HAS9 else []
def ent(aln,cols,minseq=10):
    by=defaultdict(list)
    for k,s in aln.items(): by[meta(k)[2]].append(s)
    o={}
    for y,seqs in by.items():
        if len(seqs)<minseq: continue
        v=[]
        for c in cols:
            ct=Counter(s[c] for s in seqs if s[c] in AA); n=sum(ct.values())
            if n<minseq: continue
            v.append(-sum((q/n)*math.log2(q/n) for q in ct.values()))
        if v: o[y]=round(stat.mean(v),4)
    return o
e5=ent(a5,ag5)
e9=ent(a9,ag9) if HAS9 else {}
def wt(e,bp):
    pre=[v for y,v in e.items() if y<bp]; post=[v for y,v in e.items() if y>=bp]
    if len(pre)<2 or len(post)<2: return None
    return dict(n_pre=len(pre),mean_pre=round(stat.mean(pre),4),n_post=len(post),
                mean_post=round(stat.mean(post),4),p=round(float(ranksums(pre,post)[1]),4))
R["entropy"]={"h5_by_year":{str(k):v for k,v in sorted(e5.items())},
 "h5_bp2016":wt(e5,2016),"h5_bp2021":wt(e5,2021),"n_ag_cols":len(ag5),
 "h9_by_year":{str(k):v for k,v in sorted(e9.items())},"h9_bp2021":wt(e9,2021)}
json.dump(R,open("out/NUMBERS.json","w"),indent=1,default=str)
for k in ["dataset","subtype_clade","na_within_2344b","substitutions","conservation","rbs","panel","entropy"]:
    print(f"\n=== {k.upper()} ==="); print(json.dumps(R[k],indent=1,default=str)[:1300])
