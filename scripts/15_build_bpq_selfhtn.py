"""Build the self-reported-hypertension CMM3 sensitivity variable (Section 3.7).
Reads the local NHANES BPQ_D.xpt / BPQ_E.xpt (CDC public-use, variable BPQ020 = ever told
high blood pressure) and merges it into the analytic file. Produces
nhanes_analysis_selfhtn.csv with CMM3_selfhtn = diabetes + self-reported HTN + obesity.
(Analysis only; no estimates are computed here.)"""
import pyreadstat, pandas as pd, numpy as np, os
RAW = r"H:\眼科论文\眼科2\data\raw\nhanes"
OUT = r"C:\ops2\output"
if not os.path.exists(OUT): OUT = r"H:\眼科论文\眼科2\analysis\output"

parts = []
for f in ["BPQ_D.xpt", "BPQ_E.xpt"]:
    d, _ = pyreadstat.read_xport(os.path.join(RAW, f)); d.columns = [c.upper() for c in d.columns]
    parts.append(d[["SEQN", "BPQ020"]])
b = pd.concat(parts, ignore_index=True).drop_duplicates("SEQN")
b["BPQ020"] = pd.to_numeric(b["BPQ020"], errors="coerce")
b["htn_self"] = np.where(b["BPQ020"] == 1, 1.0, np.where(b["BPQ020"] == 2, 0.0, np.nan))

nh = pd.read_csv(os.path.join(OUT, "nhanes_analytic_age40_complete.csv")); nh = nh[nh.educ.isin([1, 2, 3, 4, 5])]
m = nh.merge(b[["SEQN", "BPQ020", "htn_self"]], on="SEQN", how="left")
m["CMM3_selfhtn"] = m["diabetes"] + m["htn_self"] + m["obesity"]
print("merged", len(m), "self-HTN non-missing", int(m.htn_self.notna().sum()),
      "measured %.3f self %.3f" % (m.hypertension.mean(), m.htn_self.mean()))
m.to_csv(os.path.join(OUT, "nhanes_analysis_selfhtn.csv"), index=False)
print("wrote nhanes_analysis_selfhtn.csv")
