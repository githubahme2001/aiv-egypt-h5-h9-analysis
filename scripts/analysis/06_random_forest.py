#!/usr/bin/env python3
"""Retrain the per-position Random Forest on the EXPANDED series, with persistence
and random baselines on identical splits. Replicates the published framework:
lag-2, 48 features, RF(200, depth 6, balanced), time-series LOO + prospective holdout."""
import json, math, warnings
from collections import Counter, defaultdict
import numpy as np
from sklearn.ensemble import RandomForestClassifier
warnings.filterwarnings("ignore")
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

def coords_h5_eg():
    import json as J, re as _re
    from collections import Counter as _C
    aln=rd("out/H5_aln.faa"); W=len(next(iter(aln.values()))); N=len(aln)
    cons=[]
    for p in range(W):
        ct=_C(s[p] for s in aln.values()); g=ct.get("-",0)+ct.get("X",0)
        if g/N>=0.50: cons.append("-")
        else:
            c2=_C(s[p] for s in aln.values() if s[p] in AA)
            cons.append(c2.most_common(1)[0][0] if c2 else "-")
    cons="".join(cons); fp=cons.find("GLFGAI"); m=_re.search(r"DQICIGYHA",cons); nt=m.start() if m else 0
    c2p={};pos=0
    for c in range(nt,fp):
        if cons[c]!="-": pos+=1; c2p[c]=pos
    off=J.load(open("out/NUMBERS.json"))["dataset"]["h5_off"]
    return {v:k for k,v in c2p.items()}, off

def coords_h9(aln):
    import re
    W=len(next(iter(aln.values())))
    cons="".join((Counter(s[p] for s in aln.values() if s[p] not in "-X").most_common(1) or [("-",0)])[0][0] for p in range(W))
    fp=cons.find("GLFGAI"); m=re.search(r"[DN][KQ]ICIG",cons); nt=m.start() if m else 0
    c2p={};pos=0
    for c in range(nt,fp):
        if cons[c]!="-": pos+=1; c2p[c]=pos
    inv={v:k for k,v in c2p.items()}
    off=None
    for o in range(-12,13):
        a,b=inv.get(226-o),inv.get(228-o)
        if a is None or b is None: continue
        ra=Counter(s[a] for s in aln.values() if s[a] not in "-X").most_common(1)[0][0]
        rb=Counter(s[b] for s in aln.values() if s[b] not in "-X").most_common(1)[0][0]
        if ra=="L" and rb=="G": off=o; break
    return inv, off

def year_profiles(aln, cols, filt, minseq=10):
    """per (position, year): AA frequency vector + summary stats"""
    by=defaultdict(list)
    for k,s in aln.items():
        st,yr=k.split("|")[1],int(k.split("|")[2])
        if filt(k): by[yr].append(s)
    prof={}
    for yr,seqs in sorted(by.items()):
        if len(seqs)<minseq: continue
        for c in cols:
            ct=Counter(s[c] for s in seqs if s[c] in AA); n=sum(ct.values())
            if n<minseq: continue
            f=np.array([ct.get(a,0)/n for a in AA])
            H=-sum(x*math.log2(x) for x in f if x>0)
            dom=AA[int(np.argmax(f))]
            prof[(c,yr)]=dict(freq=f,H=H,n=n,dom=dom,domfreq=float(f.max()))
    return prof

def build_xy(prof, cols):
    """lag-2 -> 48 features, label = dominant residue at t"""
    data=defaultdict(list)
    yrs=sorted({y for (_,y) in prof})
    for c in cols:
        cy=[y for y in yrs if (c,y) in prof]
        for i,y in enumerate(cy):
            y1,y2=y-1,y-2
            if (c,y1) not in prof or (c,y2) not in prof: continue
            p1,p2=prof[(c,y1)],prof[(c,y2)]
            x=np.concatenate([[p1["domfreq"],p1["H"],p1["n"],AA.index(p1["dom"])],
                              [p2["domfreq"],p2["H"],p2["n"],AA.index(p2["dom"])],
                              p1["freq"],p2["freq"]])
            data[c].append((y,x,prof[(c,y)]["dom"],p1["dom"]))
    return data

def evaluate(data, label):
    rows=[]
    for c,recs in data.items():
        if len(recs)<4: continue
        ys=[r[0] for r in recs]; X=np.array([r[1] for r in recs])
        Y=np.array([r[2] for r in recs]); PREV=np.array([r[3] for r in recs])
        classes=sorted(set(Y))
        # leave-one-out over time points
        mc=pc=rc=0
        for i in range(len(recs)):
            tr=[j for j in range(len(recs)) if j!=i]
            if len(set(Y[tr]))<2:
                pred=Y[tr][0]
            else:
                m=RandomForestClassifier(n_estimators=200,max_depth=6,class_weight="balanced",random_state=0)
                m.fit(X[tr],Y[tr]); pred=m.predict(X[i:i+1])[0]
            mc += int(pred==Y[i]); pc += int(PREV[i]==Y[i])
            rc += 1.0/max(1,len(set(Y[tr])))
        rows.append(dict(col=c,n=len(recs),model=mc/len(recs),persist=pc/len(recs),
                         rand=rc/len(recs),n_classes=len(classes)))
    if not rows: return None
    M=np.mean([r["model"] for r in rows]); P=np.mean([r["persist"] for r in rows]); R=np.mean([r["rand"] for r in rows])
    print(f"  {label:<44} positions={len(rows):>3}  model={M:.3f}  persistence={P:.3f}  random={R:.3f}  advantage={M-P:+.3f}")
    return dict(label=label,n_positions=len(rows),model=M,persistence=P,random=R,advantage=M-P,rows=rows)

# ================= H5 =================
a5=rd("out/H5_aln.faa")
import json as _j
_m5=_j.load(open("out/h3map.json"))["col2pos"]
inv5={int(v):int(k) for k,v in _m5.items()}; OFF5=0
P5=[83,128,138,140,156,158,160,163,165,167,182,185,193,195,196,197,204,207,208,217,226,227,228,263]
c5=[inv5[p-OFF5] for p in P5 if (p-OFF5) in inv5]
import json as _j
_nc=_j.load(open("out/clades.json"))
_G221={"2.2.1","2.2.1.1","2.2.1.1a","2.2.1.2"}
def _cl(k): return _nc.get(k.split("|")[0])
def h5_all(k): return _cl(k) in _G221|{"2.3.4.4b"}
def h5_344(k): return _cl(k)=="2.3.4.4b"
def h5_221(k): return _cl(k) in _G221
print("H5 (24 monitored positions, expanded series):")
res={}
res["H5_pooled"]=evaluate(build_xy(year_profiles(a5,c5,h5_all),c5),"H5Nx pooled across clade eras")
res["H5_221"]=evaluate(build_xy(year_profiles(a5,c5,h5_221),c5),"within clade 2.2.1.x only")
res["H5_344"]=evaluate(build_xy(year_profiles(a5,c5,h5_344),c5),"within clade 2.3.4.4b only")


# ================= H9 =================
a9=rd("out/H9_aln.faa")
_m9=_j.load(open("out/h3map_H9.json"))["col2pos"]
inv9={int(v):int(k) for k,v in _m9.items()}; OFF9=0
P9=[83,95,127,134,148,155,156,158,160,182,193,195,198,226,227,228]
c9=[inv9[p-OFF9] for p in P9 if (p-OFF9) in inv9]
print(f"\nH9N2 ({len(c9)} monitored positions, offset +{OFF9}, expanded series):")
res["H9"]=evaluate(build_xy(year_profiles(a9,c9,lambda k:True),c9),"H9N2 all years")

json.dump({k:{kk:vv for kk,vv in v.items() if kk!="rows"} for k,v in res.items() if v},
          open("out/ml.json","w"),indent=1)
print("\nwrote out/ml_results.json")

# ---------- prospective holdout: train <=2022, test 2023-2025 ----------
def holdout(data,label):
    mc=pc=rc=tot=0
    for c,recs in data.items():
        tr=[r for r in recs if r[0]<=2022]; te=[r for r in recs if 2023<=r[0]<=2025]
        if len(tr)<3 or not te: continue
        X=np.array([r[1] for r in tr]); Y=np.array([r[2] for r in tr])
        for y,x,lab,prev in te:
            if len(set(Y))<2: pred=Y[0]
            else:
                m=RandomForestClassifier(n_estimators=200,max_depth=6,class_weight="balanced",random_state=0)
                m.fit(X,Y); pred=m.predict(x.reshape(1,-1))[0]
            mc+=int(pred==lab); pc+=int(prev==lab); rc+=1.0/max(1,len(set(Y))); tot+=1
    if tot:
        print(f"  {label:<44} test_pts={tot:>4}  model={mc/tot:.3f}  persistence={pc/tot:.3f}  random={rc/tot:.3f}")
        return dict(n=tot,model=mc/tot,persistence=pc/tot,random=rc/tot)
print("\nPROSPECTIVE HOLDOUT (train <=2022, test 2023-2025):")
ho={}
ho["H5_pooled"]=holdout(build_xy(year_profiles(a5,c5,h5_all),c5),"H5Nx pooled")
ho["H9"]=holdout(build_xy(year_profiles(a9,c9,lambda k:True),c9),"H9N2")
json.dump(ho,open("out/ml_holdout.json","w"),indent=1)
