#!/usr/bin/env python3
import json, csv
from collections import Counter, defaultdict
import numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
AA="ACDEFGHIKLMNPQRSTVWY"
BLUE="#0072B2"; ORANGE="#E69F00"; PINK="#CC79A7"; GREEN="#009E73"; INK="#222222"; MUT="#7a7a7a"; GRID="#dddddd"
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":8,"axes.edgecolor":MUT,
 "axes.linewidth":.7,"text.color":INK,"axes.labelcolor":INK,"xtick.color":INK,"ytick.color":INK})
SEQ=LinearSegmentedColormap.from_list("s",["#f7fbfd","#cfe3f0","#8ec4de","#3d8fc0","#0b4f78"])
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
a5=rd("out/H5_aln.faa"); a9=rd("out/H9_aln.faa")
nc=json.load(open("out/clades.json")); R=json.load(open("out/NUMBERS.json"))
import re
def coords(aln,ntre,rbs):
    W=len(next(iter(aln.values()))); N=len(aln); cons=[]
    for p in range(W):
        ct=Counter(s[p] for s in aln.values()); g=ct.get('-',0)+ct.get('X',0)
        if g/N>=.5: cons.append('-')
        else:
            c2=Counter(s[p] for s in aln.values() if s[p] in AA); cons.append(c2.most_common(1)[0][0] if c2 else '-')
    cons="".join(cons); fp=cons.find("GLFGAI"); m=re.search(ntre,cons); nt=m.start() if m else 0
    c2p={};pos=0
    for c in range(nt,fp):
        if cons[c]!='-': pos+=1; c2p[c]=pos
    inv={v:k for k,v in c2p.items()}
    for o in range(-15,16):
        A,B=inv.get(rbs[0][0]-o),inv.get(rbs[1][0]-o)
        if A is None or B is None: continue
        if Counter(s[A] for s in aln.values() if s[A] in AA).most_common(1)[0][0]==rbs[0][1] and \
           Counter(s[B] for s in aln.values() if s[B] in AA).most_common(1)[0][0]==rbs[1][1]: return c2p,inv,o
    return c2p,inv,0
import json as _j
_m5=_j.load(open("out/h3map.json"))["col2pos"]
c2p5={int(k):int(v) for k,v in _m5.items()}; inv5={v:k for k,v in c2p5.items()}; O5=0
_m9=_j.load(open("out/h3map_H9.json"))["col2pos"]
c2p9={int(k):int(v) for k,v in _m9.items()}; inv9={v:k for k,v in c2p9.items()}; O9=0
def meta(k):
    a,st,yr,src=k.split("|"); return a,st,int(yr),src
G221={"2.2.1","2.2.1.1","2.2.1.1a","2.2.1.2"}
def grp(k):
    c=nc.get(meta(k)[0])
    return "2.2.1.x" if c in G221 else ("2.3.4.4b" if c=="2.3.4.4b" else ("EA-nonGsGD" if c=="EA-nonGsGD" else None))
# ===== FIG 1 =====
cnt=defaultdict(Counter); ny=defaultdict(Counter)
for k in a5:
    g=grp(k)
    if not g: continue
    a,st,y,src=meta(k); cnt[y][g]+=1
    if g=="2.3.4.4b": ny[y][st]+=1
years=sorted(cnt)
fig,(ax1,ax2)=plt.subplots(2,1,figsize=(7.2,5.9),gridspec_kw={"height_ratios":[1.3,1]})
x=np.arange(len(years)); bot=np.zeros(len(years))
for g,col,lb in [("2.2.1.x",BLUE,"clade 2.2.1.x"),("2.3.4.4b",ORANGE,"clade 2.3.4.4b"),("EA-nonGsGD",PINK,"EA-nonGsGD")]:
    v=np.array([cnt[y][g] for y in years],float)
    if v.sum()==0: continue
    ax1.bar(x,v,.74,bottom=bot,color=col,label=lb,zorder=3,edgecolor="white",lw=.7); bot+=v
for xi,t in zip(x,bot):
    if t: ax1.text(xi,t+6,str(int(t)),ha="center",va="bottom",fontsize=6.4)
i16,i19=years.index(2016),years.index(2019)
ax1.axvspan(i16-.5,i19+.5,color="#000",alpha=.055,zorder=1)
ax1.annotate("2.3.4.4b arrives (H5N8, 2016);\ntwo clades co-circulate to 2019",
  xy=(i16-.28,cnt[2016]["2.2.1.x"]+cnt[2016]["2.3.4.4b"]/2),xytext=(i16-3.1,195),fontsize=7,ha="center",
  arrowprops=dict(arrowstyle="->",lw=.8,color=INK,connectionstyle="arc3,rad=-.26"))
ax1.set_xticks(x); ax1.set_xticklabels([str(y) for y in years],fontsize=7)
ax1.set_ylabel("QC-passed HA sequences (n)"); ax1.set_ylim(0,270); ax1.set_xlim(-.7,len(years)-.3)
ax1.yaxis.grid(True,color=GRID,lw=.6,zorder=0); ax1.set_axisbelow(True)
for s in ("top","right"): ax1.spines[s].set_visible(False)
ax1.legend(frameon=False,fontsize=7,loc="upper right",handlelength=1.4)
ax1.set_title(f"A   Egyptian H5Nx by year and WHO/WOAH/FAO clade (n = {R['dataset']['h5']:,})",
  loc="left",fontsize=8.5,fontweight="bold",pad=7)
yb=[y for y in years if sum(ny[y].values())]
xb=np.arange(len(yb))
n8=np.array([ny[y]["H5N8"] for y in yb],float); n1=np.array([ny[y]["H5N1"] for y in yb],float)
n2=np.array([ny[y]["H5N2"] for y in yb],float); n5=np.array([ny[y]["H5N5"] for y in yb],float)
tt=n8+n1+n2+n5
ax2.bar(xb,n8/tt*100,.7,color=ORANGE,hatch="///",edgecolor="white",lw=.8,label="N8",zorder=3)
ax2.bar(xb,n5/tt*100,.7,bottom=n8/tt*100,color=GREEN,edgecolor="white",lw=.8,label="N5",zorder=3)
ax2.bar(xb,n1/tt*100,.7,bottom=(n8+n5)/tt*100,color=ORANGE,edgecolor="white",lw=.8,label="N1",zorder=3)
ax2.bar(xb,n2/tt*100,.7,bottom=(n8+n5+n1)/tt*100,color=PINK,edgecolor="white",lw=.8,label="N2",zorder=3)
for xi,y in zip(xb,yb): ax2.text(xi,103,f"n={int(sum(ny[y].values()))}",ha="center",fontsize=6,color=MUT)
ax2.set_xticks(xb); ax2.set_xticklabels([str(y) for y in yb],fontsize=7)
ax2.set_ylabel("neuraminidase subtype\nwithin clade 2.3.4.4b (%)"); ax2.set_ylim(0,114); ax2.set_yticks([0,25,50,75,100])
ax2.yaxis.grid(True,color=GRID,lw=.6,zorder=0); ax2.set_axisbelow(True)
for s in ("top","right"): ax2.spines[s].set_visible(False)
ax2.legend(frameon=False,fontsize=7.5,loc="center left",bbox_to_anchor=(1.005,.5),handlelength=1.3)
ax2.set_title("B   The clade arrives on N8 (with N5) and progressively acquires N1 from 2021",
  loc="left",fontsize=8.5,fontweight="bold",pad=7)
plt.tight_layout(h_pad=1.5); plt.savefig("out/Fig1.png",dpi=400,bbox_inches="tight",facecolor="white"); plt.close()
# ===== FIG 2/3 =====
def land(aln,P,inv,off,minseq=10):
    by=defaultdict(list)
    for k,s in aln.items(): by[meta(k)[2]].append(s)
    yrs=[y for y in sorted(by) if len(by[y])>=minseq]
    G=np.full((len(P),len(yrs)),np.nan); L=[[""]*len(yrs) for _ in P]
    for i,p in enumerate(P):
        c=inv.get(p-off)
        if c is None: continue
        for j,y in enumerate(yrs):
            ct=Counter(s[c] for s in by[y] if s[c] in AA); n=sum(ct.values())
            if n<minseq: continue
            a,k2=ct.most_common(1)[0]; G[i,j]=k2/n; L[i][j]=a
    return G,L,yrs,[len(by[y]) for y in yrs]
def draw(G,L,Y,N,P,title,fn,h):
    fig,ax=plt.subplots(figsize=(min(13,1.05+.52*len(Y)),h))
    im=ax.imshow(G,cmap=SEQ,vmin=.3,vmax=1.,aspect="auto")
    for i in range(G.shape[0]):
        for j in range(G.shape[1]):
            if not np.isnan(G[i,j]):
                ax.text(j,i,L[i][j],ha="center",va="center",fontsize=6.6,
                        color="white" if G[i,j]>.72 else INK,fontweight="bold")
    ax.set_xticks(range(len(Y))); ax.set_xticklabels([f"{y}\n(n={n})" for y,n in zip(Y,N)],fontsize=6.2)
    ax.set_yticks(range(len(P))); ax.set_yticklabels(P,fontsize=7)
    ax.set_ylabel("HA1 position (H3-equivalent)")
    ax.set_title(title,loc="left",fontsize=9,fontweight="bold",pad=8)
    ax.set_xticks(np.arange(-.5,len(Y),1),minor=True); ax.set_yticks(np.arange(-.5,len(P),1),minor=True)
    ax.grid(which="minor",color="white",lw=1.1); ax.tick_params(which="minor",length=0)
    cb=fig.colorbar(im,ax=ax,fraction=.022,pad=.012); cb.set_label("dominant-residue frequency",fontsize=7)
    cb.ax.tick_params(labelsize=6.5); cb.outline.set_linewidth(.5)
    plt.tight_layout(); plt.savefig(fn,dpi=400,bbox_inches="tight",facecolor="white"); plt.close()
P5=[83,128,138,140,156,158,160,163,165,167,182,185,193,195,196,197,204,207,208,217,226,227,228,263]
P9=[83,95,127,134,148,155,156,158,160,182,193,195,198,226,227,228]
G,L,Y,N=land(a5,P5,inv5,O5); draw(G,L,Y,N,P5,f"H5Nx antigenic-site landscape (n = {R['dataset']['h5']:,}). Letter = dominant residue; shade = its frequency.","out/Fig2.png",6.4)
G,L,Y,N=land(a9,P9,inv9,O9); draw(G,L,Y,N,P9,f"H9N2 antigenic-site landscape (n = {R['dataset']['h9']:,}). Letter = dominant residue; shade = its frequency.","out/Fig3.png",4.8)
# ===== FIG 4 =====
Gg=defaultdict(list)
for k,s in a5.items():
    g=grp(k)
    if g: Gg[g].append((k,s))
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
subs=[]
for c in sorted(c2p5):
    a,fa=dom(Gg["2.2.1.x"],c); b,fb=dom(Gg["2.3.4.4b"],c)
    if a and b and a!=b and fa>=.9 and fb>=.9: subs.append((c2p5[c]+O5,c,a,b,fa,fb))
n8e=[(k,s) for k,s in a5.items() if meta(k)[1]=="H5N8" and meta(k)[2]<=2019]
rbsp=[(226,inv5[226-O5]),(228,inv5[228-O5])]
fig,ax=plt.subplots(figsize=(7.2,.29*(len(subs)+len(rbsp))+2.0))
for i,(p,c,a,b,fa,fb) in enumerate(subs):
    ax.barh(i,fa,.6,color=BLUE,zorder=3); ax.barh(i,-fb,.6,color=ORANGE,zorder=3)
    ax.text(fa+.03,i,f"{a} {fa:.2f}",va="center",fontsize=6.2)
    ax.text(-fb-.03,i,f"{fb:.2f} {b}",va="center",ha="right",fontsize=6.2)
    s=site(p); ax.text(0,i,f" {p}{'*' if s else ''} ",va="center",ha="center",fontsize=6.2,
        bbox=dict(boxstyle="round,pad=.1",fc="#fff3d6" if s else "white",ec="none"))
    if dom(n8e,c)[0]==b: ax.plot(-1.26,i,"D",ms=3.4,color="#9a6b00",zorder=5)
off=len(subs)+.7
for j,(p,c) in enumerate(rbsp):
    a,fa=dom(Gg["2.2.1.x"],c); b,fb=dom(Gg["2.3.4.4b"],c)
    ax.barh(off+j,fa,.6,color=BLUE,alpha=.4,zorder=3); ax.barh(off+j,-fb,.6,color=ORANGE,alpha=.4,zorder=3)
    ax.text(fa+.03,off+j,f"{a} {fa:.2f}",va="center",fontsize=6.2)
    ax.text(-fb-.03,off+j,f"{fb:.2f} {b}",va="center",ha="right",fontsize=6.2)
    ax.text(0,off+j,f" {p} ",va="center",ha="center",fontsize=6.2,bbox=dict(boxstyle="round,pad=.1",fc="white",ec="none"))
ax.axhline(len(subs)-.2,color=MUT,lw=.8,ls="--")
ax.text(0,len(subs)+.05,"receptor-binding positions (unchanged)",ha="center",fontsize=6.6,color=MUT)
ax.set_yticks([]); ax.set_xlim(-1.45,1.45); ax.set_ylim(-.9,off+len(rbsp))
ax.set_xticks([-1,-.5,0,.5,1]); ax.set_xticklabels(["1.0","0.5","0","0.5","1.0"],fontsize=7)
ax.set_xlabel("within-clade dominant-residue frequency")
ax.text(-.72,off+len(rbsp)-.25,"clade 2.3.4.4b",ha="center",fontsize=8,fontweight="bold",color="#9a6b00")
ax.text(.72,off+len(rbsp)-.25,"clade 2.2.1.x",ha="center",fontsize=8,fontweight="bold",color=BLUE)
for s in ("top","right","left"): ax.spines[s].set_visible(False)
ax.xaxis.grid(True,color=GRID,lw=.6,zorder=0); ax.set_axisbelow(True)
ax.legend(handles=[Line2D([],[],marker="D",ls="",color="#9a6b00",ms=4,
  label="residue already present in H5N8, 2016–2019")],frameon=False,fontsize=7,loc="lower left",bbox_to_anchor=(0,-.26))
ax.set_title(f"HA1 substitutions, clade 2.2.1.x vs 2.3.4.4b ({len(subs)} positions at ≥ 0.90 in both)\n* marks a classical antigenic site (A–E)",
  loc="left",fontsize=8.8,fontweight="bold",pad=10)
plt.tight_layout(); plt.savefig("out/Fig4.png",dpi=400,bbox_inches="tight",facecolor="white"); plt.close()
# ===== FIG 6 =====
M=json.load(open("out/ml.json"))
order=[("H5_pooled","H5Nx, pooled across clades"),("H5_221","within clade 2.2.1.x"),
       ("H5_344","within clade 2.3.4.4b"),("H9","H9N2")]
fig,ax=plt.subplots(figsize=(7.2,2.95))
for i,(k,lbl) in enumerate(order):
    r=M[k]
    ax.plot(r["random"],i,"o",ms=7,color="#bbbbbb",zorder=3)
    ax.plot(r["persistence"],i,"s",ms=7,color=BLUE,zorder=4)
    ax.plot(r["model"],i,"D",ms=7,color=ORANGE,zorder=5)
    ax.text(1.012,i,f"{r['advantage']:+.3f}",va="center",fontsize=7.2,color="#B00020" if r['advantage']<-.001 else INK)
ax.set_yticks(range(len(order))); ax.set_yticklabels([l for _,l in order],fontsize=8)
ax.set_xlim(.915,1.05); ax.set_xlabel("mean leave-one-out accuracy")
ax.text(1.012,len(order)-.40,"model − persistence",fontsize=7,color=MUT); ax.set_ylim(-.6,len(order)-.12)
ax.legend(handles=[Line2D([],[],marker="D",ls="",color=ORANGE,label="Random Forest"),
                   Line2D([],[],marker="s",ls="",color=BLUE,label="persistence baseline"),
                   Line2D([],[],marker="o",ls="",color="#bbbbbb",label="random baseline")],
  frameon=False,fontsize=7.5,loc="lower left",ncol=3,bbox_to_anchor=(0,-.40))
ax.xaxis.grid(True,color=GRID,lw=.6,zorder=0); ax.set_axisbelow(True)
for s in ("top","right","left"): ax.spines[s].set_visible(False)
ax.set_title("The Random Forest does not outperform a carry-forward baseline",loc="left",fontsize=9,fontweight="bold",pad=8)
plt.tight_layout(); plt.savefig("out/Fig6.png",dpi=400,bbox_inches="tight",facecolor="white"); plt.close()
print(f"wrote Fig1-4,6 | substitutions in Fig4: {len(subs)}")
