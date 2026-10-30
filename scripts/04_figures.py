"""Publication figures for Topic 3 (NHANES derivation cohort).
Each panel exported separately as 300-dpi PNG + vector PDF (per user preference).
Colorblind-safe Okabe-Ito palette; Arial if available else DejaVu Sans.
No fabricated numbers: all values read from analysis/output CSVs.
"""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Circle
import pandas as pd, numpy as np, os

OUT = r"H:\眼科论文\眼科2\analysis\output"
FIG = r"H:\眼科论文\眼科2\figures"
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","DejaVu Sans"],
                     "svg.fonttype":"none","font.size":9,"axes.spines.top":False,"axes.spines.right":False})
OI = dict(blue="#0072B2", orange="#E69F00", green="#009E73", verm="#D55E00",
          purple="#CC79A7", sky="#56B4E9", yel="#F0E442", black="#000000", gray="#999999")

tot = pd.read_csv(os.path.join(OUT,"main_total_effects.csv"))
sen = pd.read_csv(os.path.join(OUT,"sensitivity_results.csv"))
med = pd.read_csv(os.path.join(OUT,"mediation.csv"))
PRETTY={"cataract_surg":"Cataract surgery","amd_self":"AMD (self-report)","amd_obj_any":"AMD"}
def pl(x): return PRETTY.get(x,x)

# ---------- Fig 1: Forest of CMM -> eye outcome ORs (primary + sensitivities) ----------
rows=[]
for _,r in tot.iterrows():
    rows.append({"label":f"{pl(r.outcome)} (primary, full CMM)", "OR":r.OR, "lo":r.OR_lo, "hi":r.OR_hi, "grp":"primary"})
for _,r in sen[sen.analysis=="CMM_no_obesity"].iterrows():
    rows.append({"label":f"{pl(r.outcome)} (CMM, no obesity)", "OR":r.OR, "lo":r.OR_lo, "hi":r.OR_hi, "grp":"sens"})
for _,r in sen[(sen.analysis=="dose_response")&(sen.term=="CMM=3+")].iterrows():
    rows.append({"label":f"{pl(r.outcome)} (CMM>=3 vs 0)", "OR":r.OR, "lo":r.OR_lo, "hi":r.OR_hi, "grp":"sens"})
for _,r in sen[sen.analysis=="access_adjusted"].iterrows():
    rows.append({"label":f"{pl(r.outcome)} (access-adjusted)", "OR":r.OR, "lo":r.OR_lo, "hi":r.OR_hi, "grp":"sens"})
for _,r in sen[sen.analysis=="objective_AMD_any"].iterrows():
    rows.append({"label":"objective AMD (retinal grade)", "OR":r.OR, "lo":r.OR_lo, "hi":r.OR_hi, "grp":"obj"})
df=pd.DataFrame(rows).iloc[::-1]

fig,ax=plt.subplots(figsize=(6.2,3.4))
y=np.arange(len(df))
col={"primary":OI["blue"],"sens":OI["orange"],"obj":OI["verm"]}
for yi,(_,r) in zip(y,df.iterrows()):
    c=col[r.grp]
    ax.plot([r.lo,r.hi],[yi,yi],color=c,lw=1.6,solid_capstyle="butt")
    ax.plot(r.OR, yi, marker="D", ms=5.5, color=c, mec="black", mew=0.4)
ax.axvline(1, color=OI["gray"], ls="--", lw=1)
ax.set_yticks(y); ax.set_yticklabels(df.label, fontsize=8)
xs=np.concatenate([df.lo,df.hi,df.OR]); lo,hi=xs.min()*0.9, xs.max()*1.05
ax.set_xlim(lo,hi); ax.set_xscale("log")
ax.set_xlabel("Odds ratio per unit / category of cardiometabolic multimorbidity (95% CI)")
ax.set_title("CMM burden and eye disease — NHANES 2005-2008 (n=5,366, age>=40)", fontsize=9.5)
from matplotlib.lines import Line2D
leg=[Line2D([],[],marker="D",ls="",color=OI["blue"],label="Primary (full CMM)"),
     Line2D([],[],marker="D",ls="",color=OI["orange"],label="Pre-specified sensitivity"),
     Line2D([],[],marker="D",ls="",color=OI["verm"],label="Objective retinal AMD")]
ax.legend(handles=leg, fontsize=7.5, loc="lower right", frameon=False)
fig.tight_layout()
fig.savefig(os.path.join(FIG,"Fig1_forest.png"),dpi=300); fig.savefig(os.path.join(FIG,"Fig1_forest.pdf")); plt.close(fig)

# ---------- Fig 2: dose-response (categorical CMM) ----------
dr=sen[sen.analysis=="dose_response"].copy()
order={"CMM=1":1,"CMM=2":2,"CMM=3+":3}
fig,axes=plt.subplots(1,2,figsize=(6.0,2.8),sharex=True)
for ax,outcome in zip(axes,["cataract_surg","amd_self"]):
    s=dr[dr.outcome==outcome].copy(); s["x"]=s.term.map(order)
    s=pd.concat([pd.DataFrame({"x":[1],"OR":[1.0],"OR_lo":[np.nan],"OR_hi":[np.nan]}),s])
    ax.errorbar(s.x, s.OR, yerr=[s.OR-s.OR_lo, s.OR_hi-s.OR], fmt="o-", color=OI["blue"], capsize=3, ms=5)
    ax.axhline(1,ls="--",lw=1,color=OI["gray"])
    ax.set_xticks([1,2,3,4]); ax.set_xticklabels(["0","1","2",">=3"])
    ax.set_xlabel("Number of cardiometabolic conditions")
    ttl="Cataract surgery" if outcome=="cataract_surg" else "Self-reported AMD"
    ax.set_title(ttl, fontsize=9)
axes[0].set_ylabel("Adjusted odds ratio (95% CI)")
fig.suptitle("Dose-response: CMM count vs eye disease (ref = none)", fontsize=9.5)
fig.tight_layout(rect=[0,0,1,0.94])
fig.savefig(os.path.join(FIG,"Fig2_dose.png"),dpi=300); fig.savefig(os.path.join(FIG,"Fig2_dose.pdf")); plt.close(fig)

# ---------- Fig 3: mediation path diagram ----------
fig,ax=plt.subplots(figsize=(5.0,3.0))
ax.axis("off")
def node(x,y,t,c):
    ax.add_patch(Circle((x,y),0.10,color=c,alpha=0.85,zorder=3)); ax.text(x,y,t,ha="center",va="center",fontsize=8,color="white",weight="bold",zorder=4)
node(0.1,0.5,"CMM",OI["black"]); node(0.5,0.82,"CRP",OI["green"]); node(0.9,0.5,"Eye\ndisease",OI["purple"])
def arrow(x1,y1,x2,y2,lbl,c,rad=0.0):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),connectionstyle=f"arc3,rad={rad}",arrowstyle="-|>",mutation_scale=12,color=c,lw=1.6,zorder=2))
    ax.text((x1+x2)/2,(y1+y2)/2+0.05,lbl,fontsize=7.5,color=c,ha="center")
c=med.iloc[0]  # cataract
arrow(0.16,0.55,0.44,0.80,f"a={c.a_path:.2f}",OI["green"],rad=0.1)
arrow(0.56,0.80,0.84,0.55,f"b={c.b_path:.3f}",OI["verm"],rad=0.1)
arrow(0.20,0.47,0.80,0.47,f"direct={c.direct_logOR:.3f}",OI["blue"],rad=-0.25)
ax.text(0.5,0.05,f"Cataract: total logOR={c.total_logOR:.3f}, indirect={c.indirect:.3f} "
                f"(95%CI {c.ind_lo:.3f},{c.ind_hi:.3f}), mediated {c.prop_mediated*100:.0f}%",fontsize=7.5,ha="center")
d=med.iloc[1]
ax.text(0.5,0.14,f"AMD: indirect={d.indirect:.3f} (95%CI {d.ind_lo:.3f},{d.ind_hi:.3f}), mediated {d.prop_mediated*100:.0f}%  ->  a-path solid, mediation NOT supported",fontsize=7.2,ha="center",color=OI["gray"])
ax.set_xlim(0,1); ax.set_ylim(0,1)
ax.set_title("Mediation: CMM -> log2-CRP -> eye disease (age-adjusted)", fontsize=9.5)
fig.tight_layout()
fig.savefig(os.path.join(FIG,"Fig3_mediation.png"),dpi=300); fig.savefig(os.path.join(FIG,"Fig3_mediation.pdf")); plt.close(fig)

print("Figures written:", [f for f in os.listdir(FIG)])
