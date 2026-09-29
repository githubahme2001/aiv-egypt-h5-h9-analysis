#!/usr/bin/env python3
from Bio import Phylo
from collections import Counter, defaultdict
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
BLUE="#0072B2"; ORANGE="#E69F00"; PINK="#CC79A7"; INK="#222222"; MUT="#7a7a7a"
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":7,"text.color":INK})
t=Phylo.read("work/sub.tree","newick")
def cl(n): return n.split("_")[1].replace("-",".")
def na(n): return n.split("_")[2]
def yr(n): return n.split("_")[3]
tips=[x.name for x in t.get_terminals()]
out=[x for x in tips if cl(x).startswith("EA")]
if out: t.root_with_outgroup({"name":out[0]})
def mrca(names):
    c=t.common_ancestor(*[{"name":n} for n in names]); m=[x.name for x in c.get_terminals()]
    return c,m,Counter(cl(n) for n in m)
b44=[x for x in tips if cl(x)=="2.3.4.4b"]; g221=[x for x in tips if cl(x).startswith("2.2.1")]
cB,mB,ccB=mrca(b44); cG,mG,ccG=mrca(g221)
print(f"rooted on {out[0] if out else 'midpoint'}")
print(f"2.3.4.4b clade: {len(mB)} tips {dict(ccB)} support={cB.confidence}")
print(f"2.2.1.x clade:  {len(mG)} tips {dict(ccG)} support={cG.confidence}")
COL={"2.3.4.4b":ORANGE,"EA.nonGsGD":PINK}
def colour(n):
    c=cl(n); return COL.get(c,BLUE)
fig,ax=plt.subplots(figsize=(6.6,7.0))
Phylo.draw(t,axes=ax,do_show=False,show_confidence=False,
           label_func=lambda c: "" ,
           branch_labels=None)
# recolour tip markers manually
ycoord={}
def calc(clade,y=[0]):
    for c in clade.clades: calc(c)
    if clade.is_terminal(): y[0]+=1; ycoord[clade]=y[0]
    else: ycoord[clade]=sum(ycoord[c] for c in clade.clades)/len(clade.clades)
calc(t.root)
def depth(clade,d=0.0,acc=None):
    acc={} if acc is None else acc
    acc[clade]=d
    for c in clade.clades: depth(c,d+(c.branch_length or 0),acc)
    return acc
dep=depth(t.root)
for tip in t.get_terminals():
    ax.plot(dep[tip],ycoord[tip],"o",ms=2.6,color=colour(tip.name),
            mec="white",mew=.3,zorder=5,clip_on=False)
for c,lbl,col in [(cB,f"clade 2.3.4.4b\n(n={ccB['2.3.4.4b']}, support {cB.confidence:.2f})",ORANGE),
                  (cG,f"clade 2.2.1.x\n(n={sum(v for k,v in ccG.items() if k.startswith('2.2.1'))}, support {cG.confidence:.2f})" if cG.confidence else "clade 2.2.1.x",BLUE)]:
    ys=[ycoord[x] for x in c.get_terminals()]
    ax.annotate(lbl,xy=(dep[c],(min(ys)+max(ys))/2),xytext=(ax.get_xlim()[1]*0.78,(min(ys)+max(ys))/2),
      fontsize=7.5,color=col,fontweight="bold",va="center",
      arrowprops=dict(arrowstyle="-",lw=1.2,color=col))
ax.set_ylabel(""); ax.set_yticks([]); ax.set_xlabel("substitutions per site")
for s in ("top","right","left"): ax.spines[s].set_visible(False)
ax.legend(handles=[Line2D([],[],marker="o",ls="",color=BLUE,ms=5,label="clade 2.2.1.x"),
                   Line2D([],[],marker="o",ls="",color=ORANGE,ms=5,label="clade 2.3.4.4b"),
                   Line2D([],[],marker="o",ls="",color=PINK,ms=5,label="EA-nonGsGD (outgroup)")],
          frameon=False,fontsize=7.5,loc="lower right")
ax.set_title("Maximum-likelihood phylogeny of Egyptian H5 HA (n = 198 representatives)",
             loc="left",fontsize=9,fontweight="bold",pad=8)
plt.tight_layout(); plt.savefig("out/Fig5.png",dpi=400,bbox_inches="tight",facecolor="white")
print("wrote out/Fig5.png")
