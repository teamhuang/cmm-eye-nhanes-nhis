"""Fig4: two-cohort forest plot (NHANES + NHIS) with meta summary, per outcome."""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, pandas as pd, numpy as np, os
from matplotlib.ticker import ScalarFormatter, NullFormatter
OUT=r"H:\眼科论文\眼科2\analysis\output"; FIG=r"H:\眼科论文\眼科2\figures"
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","DejaVu Sans"],
                     "svg.fonttype":"none","font.size":9,"axes.spines.top":False,"axes.spines.right":False})
OI=dict(blue="#0072B2",orange="#E69F00",verm="#D55E00",green="#009E73",gray="#888888",black="#000")
coh=pd.read_csv(os.path.join(OUT,"cohort_results.csv")); meta=pd.read_csv(os.path.join(OUT,"meta_nhanes_nhis.csv"))
fig,axes=plt.subplots(1,2,figsize=(7.2,3.0))
for ax,outcome,title in zip(axes,["cataract","amd"],["Cataract surgery","AMD (self-report)"]):
    s=coh[coh.outcome==outcome]; m=meta[meta.outcome==outcome].iloc[0]
    rows=list(s.itertuples())+[("Meta",m.OR_re,m.re_lo,m.re_hi)]
    labels=[r.cohort for r in s.itertuples()]+["Random-effects meta"]
    ys=np.arange(len(rows))[::-1]
    for yi,r in zip(ys,list(s.itertuples())+[[None,m.OR_re,m.re_lo,m.re_hi]]):
        OR,lo,hi=(r.OR,r.lo,r.hi) if hasattr(r,'OR') else (r[1],r[2],r[3])
        ismeta=not hasattr(r,'OR')
        c=OI['verm'] if ismeta else OI['blue']
        ax.plot([lo,hi],[yi,yi],color=c,lw=1.6); ax.plot(OR,yi,marker="D" if ismeta else "o",ms=7 if ismeta else 6,color=c,mec="black",mew=.4)
    ax.axvline(1,ls="--",color=OI['gray'],lw=1)
    ax.set_yticks(ys); ax.set_yticklabels(labels,fontsize=8.5)
    ax.set_xscale("log"); ax.set_xlabel("Adjusted OR per CMM3 condition (95% CI)")
    ax.xaxis.set_major_formatter(ScalarFormatter()); ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xticks([0.9,1.0,1.2,1.5,2.0]); ax.minorticks_off()
    ax.set_title(f"{title}   (I²={m.I2:.0f}%)",fontsize=9.5)
fig.suptitle("Cross-cohort replication: cardiometabolic multimorbidity (CMM3) and eye disease", fontsize=10)
fig.tight_layout(rect=[0,0,1,0.93])
fig.savefig(os.path.join(FIG,"Fig4_meta.png"),dpi=300); fig.savefig(os.path.join(FIG,"Fig4_meta.pdf")); plt.close(fig)
print("wrote Fig4")
