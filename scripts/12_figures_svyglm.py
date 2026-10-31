"""Nature-style figures: each panel is a 4 x 4 cm forest plot with OR (95% CI) and P
printed to its right; panels are assembled (Figures 1-2 side by side, Figure 3 2+1).
No in-figure notes beyond panel names/letters.
"""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.transforms as mtrans, pandas as pd, numpy as np, os
from matplotlib.ticker import ScalarFormatter, NullFormatter

OUT = r"H:\眼科论文\眼科2\analysis\output"; FIG = r"H:\眼科论文\眼科2\figures"
os.makedirs(FIG, exist_ok=True)
CM = 1/2.54
plt.rcParams.update({
    "font.family":"sans-serif", "font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
    "font.size":7.2, "axes.labelsize":6.8, "axes.titlesize":7.2, "xtick.labelsize":6.5,
    "ytick.labelsize":6.6, "legend.fontsize":6.2, "axes.linewidth":0.5,
    "xtick.major.width":0.5, "ytick.major.width":0.5, "xtick.major.size":1.8,
    "ytick.major.size":1.8, "xtick.direction":"out", "ytick.direction":"out",
    "axes.spines.top":False, "axes.spines.right":False, "savefig.dpi":600,
    "pdf.fonttype":42, "ps.fonttype":42,
})
C = dict(blue="#2E5A88", orange="#C46A1B", green="#3A7D44", purple="#7E5A8C",
         grey="#8A8A8A", red="#A63A3A")

def axes_cm(fig, l, b, w=4.0, h=4.0):
    W, H = fig.get_size_inches() * 2.54
    return fig.add_axes([l/W, b/H, w/W, h/H])

def alabel(ax, letter):
    ax.text(-0.16, 1.03, letter, transform=ax.transAxes, fontsize=8.5, fontweight="bold",
            va="bottom", ha="left", clip_on=False)

def fp(p):
    if p is None or not np.isfinite(p): return "\u2014"
    if p < 0.001: return "<0.001"
    if p < 0.01: return f"{p:.3f}"
    return f"{p:.2f}"

def forest_panel(fig, x0, y0, rows, title, letter, left_cm=1.5,
                 w_cm=4.0, h_cm=4.0, p_offset_cm=2.1,
                 xlim=(0.45,3.6), xticks=(0.5,1,2,3), xlog=True, null=1.0,
                 xlabel="OR", header="OR (95% CI)", diff=False):
    ax = axes_cm(fig, x0 + left_cm, y0, w_cm, h_cm)
    p_x = 1.0 + p_offset_cm / w_cm
    ys = np.arange(len(rows))[::-1]
    for yi,(lab,OR,lo,hi,col,mk,p) in zip(ys,rows):
        ax.plot([lo,hi],[yi,yi], color=col, lw=0.9, solid_capstyle="butt", zorder=2)
        ax.plot(OR, yi, marker=mk, ms=3.2, color=col, mec="white", mew=0.4, zorder=3)
    ax.axvline(null, color=C["grey"], lw=0.5, ls=(0, (3,2)), zorder=1)
    ax.set_yticks(ys); ax.set_yticklabels([r[0] for r in rows])
    ax.set_ylim(-0.65, len(rows)-0.35)
    if xlog:
        ax.set_xscale("log"); ax.set_xlim(*xlim); ax.set_xticks(xticks)
        ax.xaxis.set_major_formatter(ScalarFormatter()); ax.xaxis.set_minor_formatter(NullFormatter())
        ax.minorticks_off()
    else:
        ax.set_xlim(*xlim); ax.set_xticks(xticks)
        ax.xaxis.set_major_formatter(ScalarFormatter()); ax.minorticks_off()
    ax.set_xlabel(xlabel, labelpad=1.5)
    ax.set_title(title, pad=2)
    alabel(ax, letter)
    tr = mtrans.blended_transform_factory(ax.transAxes, ax.transData)
    ax.text(1.03, len(rows)-0.30, header, transform=tr, fontsize=5.9, fontweight="bold", va="bottom", ha="left")
    ax.text(p_x, len(rows)-0.30, "P", transform=tr, fontsize=5.9, fontweight="bold", va="bottom", ha="left")
    for yi,(lab,OR,lo,hi,col,mk,p) in zip(ys,rows):
        s = f"{OR:+.3f} ({lo:+.3f}, {hi:+.3f})" if diff else f"{OR:.2f} ({lo:.2f}\u2013{hi:.2f})"
        ax.text(1.03, yi, s, transform=tr, fontsize=5.9, va="center", ha="left")
        ax.text(p_x, yi, fp(p), transform=tr, fontsize=5.9, va="center", ha="left")
    return ax

# ---------------------------------------------------------------- data
prim = pd.read_csv(os.path.join(OUT,"sr_nhanes_primary_svyglm.csv"))
dec  = pd.read_csv(os.path.join(OUT,"sr_cmm3_vs_cmm4.csv"))
sens = pd.read_csv(os.path.join(OUT,"sr_sensitivity_svyglm.csv"))
r5   = pd.read_csv(os.path.join(OUT,"sr2_r5_corrected_models.csv"))
coh  = pd.read_csv(os.path.join(OUT,"sr_cohort_cmm3_svyglm.csv"))
meta = pd.read_csv(os.path.join(OUT,"sr2_meta_audit.csv"))
strat= pd.read_csv(os.path.join(OUT,"sr_age_strata.csv"))
boot = pd.read_csv(os.path.join(OUT,"sr2_mediation_boot.csv"))

def pget(df, outcome, key="outcome", **kw):
    d=df
    for k,v in kw.items(): d=d[d[k]==v]
    c="outcome" if "outcome" in d.columns else "term"
    return d[d[c]==outcome].iloc[0]
def crow(cohort,outcome):
    r=coh[(coh.cohort==cohort)&(coh.outcome==outcome)].iloc[0]
    return r.OR,r.lo,r.hi,r.p

# ---------------------------------------------------------------- Figure 1
r = pget(prim,"cataract",model="CMM4_total"); c4c=(r.OR,r.lo,r.hi,r.p)
r = pget(prim,"amd",model="CMM4_total");      c4a=(r.OR,r.lo,r.hi,r.p)
c3c = crow("NHANES","cataract"); c3a = crow("NHANES","amd")
dc = dec[dec.outcome=="cataract"].iloc[0]; da = dec[dec.outcome=="amd"].iloc[0]
r = sens[(sens.analysis=="no_obesity_CMM2")&(sens.outcome=="cataract")].iloc[0]; dm=(r.OR,r.lo,r.hi,r.p)
r = sens[(sens.analysis=="access_adjusted")&(sens.outcome=="cataract")].iloc[0]; acc_c=(r.OR,r.lo,r.hi,r.p)
r = sens[(sens.analysis=="access_adjusted")&(sens.outcome=="amd")].iloc[0]; acc_a=(r.OR,r.lo,r.hi,r.p)
o = r5[r5.analysis=="amd_obj_any~CMM3_meas"].iloc[0]; obj=(o.OR,o.lo,o.hi,o.p)

rows_c = [("CMM4",*c4c[:3],C["blue"],"o",c4c[3]),
          ("CMM3",*c3c[:3],C["green"],"o",c3c[3]),
          ("Dyslipidaemia",dc.dys_OR,dc.dys_lo,dc.dys_hi,C["purple"],"o",dc.dys_p),
          ("Diabetes+HTN",*dm[:3],C["orange"],"o",dm[3]),
          ("CMM3+access",*acc_c[:3],C["grey"],"o",acc_c[3])]
rows_a = [("CMM4",*c4a[:3],C["blue"],"o",c4a[3]),
          ("CMM3",*c3a[:3],C["green"],"o",c3a[3]),
          ("Dyslipidaemia",da.dys_OR,da.dys_lo,da.dys_hi,C["purple"],"o",da.dys_p),
          ("CMM3+access",*acc_a[:3],C["grey"],"o",acc_a[3]),
          ("Objective AMD",*obj[:3],C["red"],"s",obj[3])]
fig = plt.figure(figsize=(17.0*CM, 6.4*CM))
forest_panel(fig, 0.2, 1.2, rows_c, "Cataract surgery", "a", left_cm=1.5)
forest_panel(fig, 8.7, 1.2, rows_a, "AMD", "b", left_cm=1.5)
fig.savefig(os.path.join(FIG,"Fig1_nhanes_forest.png"), dpi=600, bbox_inches="tight", pad_inches=0.02)
fig.savefig(os.path.join(FIG,"Fig1_nhanes_forest.pdf"), bbox_inches="tight", pad_inches=0.02); plt.close(fig); print("Fig1 done")

# ---------------------------------------------------------------- Figure 2
tot = prim[prim.model=="CMM4_total"]; dirn = prim[prim.model=="CMM4_direct_adjCRP"]
def pr(sub,term): x=sub[sub.term==term].iloc[0]; return x.OR,x.lo,x.hi,x.p
ct=pr(tot,"cataract"); ca=pr(tot,"amd"); dt=pr(dirn,"cataract"); da2=pr(dirn,"amd")
rows2a=[("Cataract total",*ct[:3],C["blue"],"o",ct[3]),
        ("Cataract CRP-adj.",*dt[:3],C["orange"],"o",dt[3]),
        ("AMD total",*ca[:3],C["blue"],"o",ca[3]),
        ("AMD CRP-adj.",*da2[:3],C["orange"],"o",da2[3])]
ind_c=boot["tot_cat"]-boot["dir_cat"]; ind_a=boot["tot_amd"]-boot["dir_amd"]
def bp(v):
    v=v[np.isfinite(v)]; return min(1.0,2*min((v<=0).mean(),(v>=0).mean()))
ic=(np.mean(ind_c),np.quantile(ind_c,.025),np.quantile(ind_c,.975),bp(ind_c))
ia=(np.mean(ind_a),np.quantile(ind_a,.025),np.quantile(ind_a,.975),bp(ind_a))
rows2b=[("Cataract",*ic[:3],C["red"],"o",ic[3]),("AMD",*ia[:3],C["red"],"o",ia[3])]
fig = plt.figure(figsize=(17.0*CM, 6.4*CM))
forest_panel(fig, 0.2, 1.2, rows2a, "Total vs CRP-adjusted", "a", left_cm=1.8)
forest_panel(fig, 8.7, 1.2, rows2b, "Difference", "b", left_cm=1.2, p_offset_cm=2.6,
             xlim=(-0.07,0.05), xticks=(-0.06,-0.04,-0.02,0.0,0.02,0.04),
             xlog=False, null=0.0, xlabel="Difference", header="Difference (95% CI)",
             diff=True)
fig.savefig(os.path.join(FIG,"Fig2_mediation.png"), dpi=600, bbox_inches="tight", pad_inches=0.02)
fig.savefig(os.path.join(FIG,"Fig2_mediation.pdf"), bbox_inches="tight", pad_inches=0.02); plt.close(fig); print("Fig2 done")

# ---------------------------------------------------------------- Figure 3
def crows(outcome):
    rows=[]
    for name,col in [("NHANES",C["blue"]),("NHIS",C["green"])]:
        rows.append((name,*crow(name,outcome)[:3],col,"o",crow(name,outcome)[3]))
    m=meta[meta.outcome==outcome].iloc[0]
    rows.append(("Fixed",float(m.OR_fe),float(m.fe_lo),float(m.fe_hi),C["orange"],"D",float(m.p_fe)))
    rows.append(("Random",float(m.OR_re_hk),float(m.re_lo),float(m.re_hi),C["red"],"D",float(m.p_re_hk)))
    return rows
rows3a=crows("cataract"); rows3b=crows("amd")
rows3c=[]
for age,lab in [("40-59","40/45\u201359"),("60-69","60\u201369"),("70+","\u2265 70")]:
    for name,col in [("NHANES",C["blue"]),("NHIS",C["green"])]:
        r=strat[(strat.cohort==name)&(strat.outcome=="cataract")&(strat.agegrp==age)].iloc[0]
        rows3c.append((f"{lab}, {name}",r.OR,r.lo,r.hi,col,"o",r.p))
fig = plt.figure(figsize=(17.0*CM, 12.2*CM))
forest_panel(fig, 0.2, 7.0, rows3a, "Cataract surgery", "a", left_cm=1.4)
forest_panel(fig, 8.7, 7.0, rows3b, "AMD", "b", left_cm=1.4)
forest_panel(fig, 0.2, 0.9, rows3c, "By age stratum", "c", left_cm=2.0, w_cm=8.0)
fig.savefig(os.path.join(FIG,"Fig3_cross_survey.png"), dpi=600, bbox_inches="tight", pad_inches=0.02)
fig.savefig(os.path.join(FIG,"Fig3_cross_survey.pdf"), bbox_inches="tight", pad_inches=0.02); plt.close(fig); print("Fig3 done")
