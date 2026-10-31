"""Nature-style figures, rebuilt from the verified CSVs. Rules:
- height 5 cm; width 6.5 cm (multi-panel) or 5 cm (single column).
- no in-figure title, no explanatory notes; only axes, tick labels, panel letters and legends.
- Helvetica/Arial, 5.5-7 pt, thin lines, left/bottom spines, colour-blind-safe muted palette.
"""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, pandas as pd, numpy as np, os
from matplotlib.patches import Patch
from matplotlib.ticker import ScalarFormatter, NullFormatter

OUT = r"H:\眼科论文\眼科2\analysis\output"; FIG = r"H:\眼科论文\眼科2\figures"
os.makedirs(FIG, exist_ok=True)
CM = 1/2.54
plt.rcParams.update({
    "font.family":"sans-serif", "font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
    "font.size":5.2, "axes.labelsize":5.5, "axes.titlesize":5.8, "xtick.labelsize":5.0,
    "ytick.labelsize":5.0, "legend.fontsize":4.8, "axes.linewidth":0.5,
    "xtick.major.width":0.5, "ytick.major.width":0.5, "xtick.major.size":2,
    "ytick.major.size":2, "xtick.direction":"out", "ytick.direction":"out",
    "axes.spines.top":False, "axes.spines.right":False, "savefig.dpi":600,
    "pdf.fonttype":42, "ps.fonttype":42,
})
C = dict(blue="#2E5A88", orange="#C46A1B", green="#3A7D44", purple="#7E5A8C",
         grey="#8A8A8A", red="#A63A3A")

def axes_cm(fig, l, b, w=5.0, h=5.0):
    W, H = fig.get_size_inches() * 2.54
    return fig.add_axes([l/W, b/H, w/W, h/H])

def alabel(ax, letter):
    ax.text(-0.08, 1.04, letter, transform=ax.transAxes, fontsize=7, fontweight="bold",
            va="bottom", ha="left", clip_on=False)

def panel(ax, letter, dx=-0.100, dy=0.000):
    b = ax.get_position()
    ax.figure.text(b.x0 + dx, b.y1 + dy, letter, fontsize=7, fontweight="bold",
                   va="bottom", ha="left", clip_on=False)

def forest(ax, rows, xlabel="OR (95% CI)", xlim=(0.45, 3.6), xticks=(0.5, 1, 2, 3)):
    ys = np.arange(len(rows))[::-1]
    for yi, (lab, OR, lo, hi, col, mk) in zip(ys, rows):
        ax.plot([lo, hi], [yi, yi], color=col, lw=0.9, solid_capstyle="butt", zorder=2)
        ax.plot(OR, yi, marker=mk, ms=3.6, color=col, mec="white", mew=0.4, zorder=3)
    ax.axvline(1, color=C["grey"], lw=0.5, ls=(0, (3, 2)), zorder=1)
    ax.set_yticks(ys); ax.set_yticklabels([r[0] for r in rows])
    ax.set_ylim(-0.7, len(rows) - 0.3)
    ax.set_xscale("log"); ax.set_xlim(*xlim); ax.set_xticks(xticks)
    ax.xaxis.set_major_formatter(ScalarFormatter()); ax.xaxis.set_minor_formatter(NullFormatter())
    ax.minorticks_off(); ax.set_xlabel(xlabel, labelpad=2)

# ---------------------------------------------------------------- Figure 1 (6.5 x 5 cm)
prim = pd.read_csv(os.path.join(OUT,"sr_nhanes_primary_svyglm.csv"))
dec  = pd.read_csv(os.path.join(OUT,"sr_cmm3_vs_cmm4.csv"))
sens = pd.read_csv(os.path.join(OUT,"sr_sensitivity_svyglm.csv"))
r5   = pd.read_csv(os.path.join(OUT,"sr2_r5_corrected_models.csv"))
rows_c, rows_a = [], []
r = prim[(prim.model=="CMM4_total") & (prim.term=="cataract")].iloc[0]; rows_c.append(("CMM4", r.OR, r.lo, r.hi, C["blue"], "o"))
r = prim[(prim.model=="CMM4_total") & (prim.term=="amd")].iloc[0];     rows_a.append(("CMM4", r.OR, r.lo, r.hi, C["blue"], "o"))
d = dec[dec.outcome=="cataract"].iloc[0]
rows_c += [("CMM3", d.CMM3_OR, d.CMM3_lo, d.CMM3_hi, C["green"], "o"),
           ("Dyslipidaemia", d.dys_OR, d.dys_lo, d.dys_hi, C["purple"], "o")]
d = dec[dec.outcome=="amd"].iloc[0]
rows_a += [("CMM3", d.CMM3_OR, d.CMM3_lo, d.CMM3_hi, C["green"], "o"),
           ("Dyslipidaemia", d.dys_OR, d.dys_lo, d.dys_hi, C["purple"], "o")]
r = sens[(sens.analysis=="no_obesity_CMM2") & (sens.outcome=="cataract")].iloc[0]
rows_c.append(("Diabetes+HTN", r.OR, r.lo, r.hi, C["orange"], "o"))
r = sens[(sens.analysis=="access_adjusted") & (sens.outcome=="cataract")].iloc[0]
rows_c.append(("CMM3+access", r.OR, r.lo, r.hi, C["grey"], "o"))
r = sens[(sens.analysis=="access_adjusted") & (sens.outcome=="amd")].iloc[0]
rows_a.append(("CMM3+access", r.OR, r.lo, r.hi, C["grey"], "o"))
o = r5[r5.analysis=="amd_obj_any~CMM3_meas"].iloc[0]
rows_a.append(("Objective AMD", o.OR, o.lo, o.hi, C["red"], "s"))

# two 5 x 5 cm panels, assembled side by side
fig = plt.figure(figsize=(13.9*CM, 7.0*CM))
ax = axes_cm(fig, 2.4, 1.3); forest(ax, rows_c); ax.set_title("Cataract surgery", pad=2); alabel(ax, "a")
ax = axes_cm(fig, 8.5, 1.3); forest(ax, rows_a); ax.set_title("AMD", pad=2); alabel(ax, "b")
fig.savefig(os.path.join(FIG,"Fig1_nhanes_forest.png"), dpi=600)
fig.savefig(os.path.join(FIG,"Fig1_nhanes_forest.pdf")); plt.close(fig); print("Fig1 done")

# ---------------------------------------------------------------- Figure 2 (6.5 x 5 cm, side by side)
tot = prim[prim.model=="CMM4_total"]; dirn = prim[prim.model=="CMM4_direct_adjCRP"]
boot = pd.read_csv(os.path.join(OUT,"sr2_mediation_boot_summary.csv"), index_col=0)
# two 5 x 5 cm panels, assembled side by side
fig = plt.figure(figsize=(13.4*CM, 7.0*CM))
ax = axes_cm(fig, 1.8, 1.4)
xo = np.arange(2); w = 0.32
vals = {}
for nm, sub in [("total", tot), ("direct", dirn)]:
    vals[nm] = dict(OR=[float(sub[sub.term=="cataract"].OR.iloc[0]), float(sub[sub.term=="amd"].OR.iloc[0])],
                    lo=[float(sub[sub.term=="cataract"].lo.iloc[0]), float(sub[sub.term=="amd"].lo.iloc[0])],
                    hi=[float(sub[sub.term=="cataract"].hi.iloc[0]), float(sub[sub.term=="amd"].hi.iloc[0])])
for x, nm, col, off in [(xo-w/2,"total",C["blue"],-1),(xo+w/2,"direct",C["orange"],1)]:
    v = vals[nm]; err = [np.array(v["OR"])-v["lo"], np.array(v["hi"])-v["OR"]]
    ax.bar(x, v["OR"], width=w, color=col, edgecolor="none", yerr=err, capsize=1.5,
           error_kw=dict(elinewidth=0.6, capthick=0.6, ecolor="#333333"))
ax.axhline(1, color=C["grey"], lw=0.5, ls=(0, (3, 2)))
ax.set_xticks(xo); ax.set_xticklabels(["Cataract","AMD"]); ax.set_ylabel("CMM4 OR")
ax.set_ylim(0.85,1.55); ax.set_yticks([0.9,1.0,1.1,1.2,1.3,1.4])
ax.legend([Patch(color=C["blue"]), Patch(color=C["orange"])], ["Total","CRP-adjusted"], frameon=False,
          loc="upper left", ncol=1, handlelength=0.8, labelspacing=0.3, borderaxespad=0.1)
alabel(ax, "a")
ax = axes_cm(fig, 8.1, 1.4)
ind = [float(boot.loc["indirect_cataract","est"]), float(boot.loc["indirect_amd","est"])]
ilo = [float(boot.loc["indirect_cataract","lo"]), float(boot.loc["indirect_amd","lo"])]
ihi = [float(boot.loc["indirect_cataract","hi"]), float(boot.loc["indirect_amd","hi"])]
ax.errorbar(ind, [1,0], xerr=[np.array(ind)-ilo, np.array(ihi)-ind], fmt="o", ms=3.6,
            color=C["red"], elinewidth=0.9, capsize=1.8, mec="white", mew=0.4)
ax.axvline(0, color=C["grey"], lw=0.5, ls=(0, (3, 2)))
ax.set_yticks([0,1]); ax.set_yticklabels(["AMD","Cataract"]); ax.set_ylim(-0.6,1.6)
ax.set_xlim(-0.07,0.05); ax.set_xlabel("Difference in log-odds", labelpad=2)
ax.xaxis.set_major_formatter(ScalarFormatter()); ax.minorticks_off()
alabel(ax, "b")
fig.savefig(os.path.join(FIG,"Fig2_mediation.png"), dpi=600)
fig.savefig(os.path.join(FIG,"Fig2_mediation.pdf")); plt.close(fig); print("Fig2 done")

# ---------------------------------------------------------------- Figure 3 (6.5 x 5 cm, 2x2 with c spanning)
coh = pd.read_csv(os.path.join(OUT,"sr_cohort_cmm3_svyglm.csv"))
meta = pd.read_csv(os.path.join(OUT,"sr2_meta_audit.csv"))
strat = pd.read_csv(os.path.join(OUT,"sr_age_strata.csv"))
fig = plt.figure(figsize=(6.5*CM, 5*CM))
gs = fig.add_gridspec(2, 2, hspace=0.9, wspace=0.75, left=0.215, right=0.97, top=0.92, bottom=0.19)
for j, outcome in enumerate(["cataract","amd"]):
    ax = fig.add_subplot(gs[0,j])
    s = coh[coh.outcome==outcome]; m = meta[meta.outcome==outcome].iloc[0]
    rows = [("NHANES", float(s[s.cohort=="NHANES"].OR.iloc[0]), float(s[s.cohort=="NHANES"].lo.iloc[0]), float(s[s.cohort=="NHANES"].hi.iloc[0]), C["blue"], "o"),
            ("NHIS", float(s[s.cohort=="NHIS"].OR.iloc[0]), float(s[s.cohort=="NHIS"].lo.iloc[0]), float(s[s.cohort=="NHIS"].hi.iloc[0]), C["green"], "o"),
            ("Fixed", float(m.OR_fe), float(m.fe_lo), float(m.fe_hi), C["orange"], "D"),
            ("Random", float(m.OR_re_hk), float(m.re_lo), float(m.re_hi), C["red"], "D")]
    forest(ax, rows, xlabel="OR per condition", xlim=(0.45,3.6), xticks=(0.5,1,2))
    ax.set_title("Cataract surgery" if outcome=="cataract" else "AMD", pad=3)
    panel(ax, "ab"[j])
ax = fig.add_subplot(gs[1,:])
cats = ["40\u201359 / 45\u201359","60\u201369","\u2265 70"]
yb = np.arange(len(cats))[::-1]
for k,(cohort,col) in enumerate([("NHANES",C["blue"]),("NHIS",C["green"])]):
    ss = strat[(strat.cohort==cohort)&(strat.outcome=="cataract")].set_index("agegrp").reindex(["40-59","60-69","70+"])
    yy = yb + (0.16 if k==0 else -0.16)
    ax.errorbar(ss.OR, yy, xerr=[ss.OR-ss.lo, ss.hi-ss.OR], fmt="o", ms=3.6, color=col,
                elinewidth=0.9, capsize=1.8, mec="white", mew=0.4, label=cohort)
ax.axvline(1, color=C["grey"], lw=0.5, ls=(0,(3,2)))
ax.set_xscale("log"); ax.set_xlim(0.45,3.6); ax.set_xticks([0.5,1,2,3])
ax.xaxis.set_major_formatter(ScalarFormatter()); ax.minorticks_off()
ax.set_yticks(yb); ax.set_yticklabels(cats); ax.set_ylim(-0.75,2.75)
ax.set_xlabel("OR per condition (95% CI)", labelpad=2)
ax.legend(frameon=False, loc="lower right", handlelength=0.9, borderaxespad=0.2)
panel(ax, "c")
fig.savefig(os.path.join(FIG,"Fig3_cross_survey.png"), dpi=600)
fig.savefig(os.path.join(FIG,"Fig3_cross_survey.pdf")); plt.close(fig); print("Fig3 done")
