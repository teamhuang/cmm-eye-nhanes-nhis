"""Supplementary Figures S1-S3 (three-line tables companion)."""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, pandas as pd, numpy as np, os
OUT=r"H:\眼科论文\眼科2\analysis\output"; FIG=r"H:\眼科论文\眼科2\figures"
CM=1/2.54
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","DejaVu Sans"],
 "font.size":7,"axes.labelsize":7,"axes.titlesize":7.5,"xtick.labelsize":6.5,"ytick.labelsize":6.5,
 "axes.linewidth":0.6,"xtick.major.width":0.6,"ytick.major.width":0.6,"xtick.major.size":2,
 "axes.spines.top":False,"axes.spines.right":False,"savefig.dpi":600,"pdf.fonttype":42})
C=dict(blue="#2E5A88",orange="#C46A1B",green="#3A7D44",grey="#8A8A8A",red="#A63A3A",purple="#7E5A8C")

# S1 sample flow
steps=[("Examined, aged \u2265 40 y",7081),("Valid cataract-surgery data",7074),
       ("Valid AMD data",7005),("Complete CMM4",5773),("Complete CRP",5770),
       ("Complete covariates",5366)]
fig,ax=plt.subplots(figsize=(8*CM,9*CM)); ax.set_axis_off()
ax.set_xlim(0,1); ax.set_ylim(0,len(steps)+0.5)
for i,(lab,n) in enumerate(steps):
    y=len(steps)-i-0.4
    ax.add_patch(plt.Rectangle((0.18,y-0.28),0.64,0.56,facecolor="#EEF2F7",edgecolor=C["blue"],lw=0.7))
    ax.text(0.50,y+0.08,f"{lab}",ha="center",va="center",fontsize=6.8)
    ax.text(0.50,y-0.14,f"n = {n:,}",ha="center",va="center",fontsize=6.5,color="#333333")
    if i<len(steps)-1:
        d=steps[i+1][1]-n
        ax.annotate("",xy=(0.50,y-0.72),xytext=(0.50,y-0.30),arrowprops=dict(arrowstyle="-|>",lw=0.6,color=C["grey"]))
        if d<0: ax.text(0.53,y-0.51,f"\u2212{abs(d):,}",ha="left",va="center",fontsize=6,color="#555555")
fig.savefig(os.path.join(FIG,"FigS1_flow.png"),bbox_inches="tight",pad_inches=0.02)
fig.savefig(os.path.join(FIG,"FigS1_flow.pdf"),bbox_inches="tight",pad_inches=0.02); plt.close(fig)
print("S1 done")

# S2 CMM3 x dyslipidaemia
jd=pd.read_csv(os.path.join(OUT,"sr2_cmm3_dys_joint.csv"))
jd["CMM3"]=jd["CMM3"].astype(int); jd["dys"]=jd["dyslipidemia"].astype(int)
pct=jd.pivot(index="CMM3",columns="dys",values="pct_within_CMM3").fillna(0)
fig,ax=plt.subplots(figsize=(8*CM,5*CM))
ax.bar(pct.index,pct[1],width=0.6,color=C["purple"],edgecolor="none")
for x,v in zip(pct.index,pct[1]): ax.text(x,v+1,f"{v:.1f}%",ha="center",fontsize=6)
ax.set_xticks(pct.index); ax.set_xlabel("Number of other cardiometabolic conditions (CMM3)")
ax.set_ylabel("Dyslipidaemia present (%)"); ax.set_ylim(0,65)
fig.savefig(os.path.join(FIG,"FigS2_cmm3_dys.png"),bbox_inches="tight",pad_inches=0.02)
fig.savefig(os.path.join(FIG,"FigS2_cmm3_dys.pdf"),bbox_inches="tight",pad_inches=0.02); plt.close(fig)
print("S2 done")

# S3 bootstrap distributions
b=pd.read_csv(os.path.join(OUT,"sr2_mediation_boot.csv"))
ic=b["tot_cat"]-b["dir_cat"]; ia=b["tot_amd"]-b["dir_amd"]
fig,ax=plt.subplots(figsize=(8*CM,5*CM))
bins=np.linspace(-0.09,0.06,45)
ax.hist(ic,bins=bins,density=True,histtype="step",lw=1.0,color=C["blue"],label="Cataract")
ax.hist(ia,bins=bins,density=True,histtype="step",lw=1.0,color=C["orange"],label="AMD")
ax.axvline(0,color=C["grey"],lw=0.6,ls=(0,(3,2)))
ax.set_xlabel("Difference in log-odds (total minus CRP-adjusted)"); ax.set_ylabel("Density")
ax.legend(frameon=False,loc="upper left",fontsize=6)
fig.savefig(os.path.join(FIG,"FigS3_bootstrap.png"),bbox_inches="tight",pad_inches=0.02)
fig.savefig(os.path.join(FIG,"FigS3_bootstrap.pdf"),bbox_inches="tight",pad_inches=0.02); plt.close(fig)
print("S3 done")
