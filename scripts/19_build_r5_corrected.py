"""Round-5 corrections: rebuild objective AMD with the CORRECT worse-eye variables,
add self-reported BMI from the NHANES Weight History (WHQ) files, and add a common-age (>=45)
flag. Reads local XPT files; writes nhanes_analysis_r5.csv. No estimates computed here."""
import pyreadstat, pandas as pd, numpy as np, os
RAW = r"H:\眼科论文\眼科2\data\raw\nhanes"; OUT = r"C:\ops2\output"
if not os.path.exists(os.path.join(OUT, "nhanes_analytic_age40_complete.csv")): OUT = r"H:\眼科论文\眼科2\analysis\output"

def rd(f):
    d, _ = pyreadstat.read_xport(os.path.join(RAW, f)); d.columns = [c.upper() for c in d.columns]; return d

# --- objective AMD: WORSE-EYE variables (OPDUARM 0/1/2; OPDULARM 0/1) ---
opx = pd.concat([rd("OPXRET_D.xpt"), rd("OPXRET_E.xpt")], ignore_index=True)
opx = opx[["SEQN","OPDUARM","OPDULARM","OPDDARM","OPDDLARM"]].drop_duplicates("SEQN")
for v in ["OPDUARM","OPDULARM","OPDDARM","OPDDLARM"]:
    opx[v] = pd.to_numeric(opx[v], errors="coerce")
opx["amd_obj_any"]     = np.where(opx["OPDUARM"].notna(),  (opx["OPDUARM"] >= 1).astype(float), np.nan)   # any ARM, worse eye
opx["amd_obj_late"]    = np.where(opx["OPDULARM"].notna(), (opx["OPDULARM"] == 1).astype(float), np.nan)  # late ARM, worse eye
opx["amd_obj_any_RE"]  = np.where(opx["OPDDARM"].notna(),  (opx["OPDDARM"] >= 1).astype(float), np.nan)    # right-eye any (old approach)

# --- self-reported height/weight (WHQ): WHD010 inches, WHD020 pounds; 9999/7777/777/999 = missing ---
whq = pd.concat([rd("WHQ_D.xpt"), rd("WHQ_E.xpt")], ignore_index=True)
whq = whq[["SEQN","WHD010","WHD020"]].drop_duplicates("SEQN")
for v in ["WHD010","WHD020"]:
    whq[v] = pd.to_numeric(whq[v], errors="coerce")
bad = {9999, 7777, 777, 999, 8888, 888}
whq.loc[whq["WHD010"].isin(bad) | (whq["WHD010"] < 48) | (whq["WHD010"] > 84), "WHD010"] = np.nan
whq.loc[whq["WHD020"].isin(bad) | (whq["WHD020"] < 70) | (whq["WHD020"] > 700), "WHD020"] = np.nan
whq["bmi_self"] = 703.0 * whq["WHD020"] / (whq["WHD010"] ** 2)
whq["obese_self"] = np.where(whq["bmi_self"].notna(), (whq["bmi_self"] >= 30).astype(float), np.nan)

# --- self-reported hypertension (BPQ020) ---
bpq = pd.concat([rd("BPQ_D.xpt"), rd("BPQ_E.xpt")], ignore_index=True)[["SEQN","BPQ020"]].drop_duplicates("SEQN")
bpq["BPQ020"] = pd.to_numeric(bpq["BPQ020"], errors="coerce")
bpq["htn_self"] = np.where(bpq["BPQ020"] == 1, 1.0, np.where(bpq["BPQ020"] == 2, 0.0, np.nan))

nh = pd.read_csv(os.path.join(OUT, "nhanes_analytic_age40_complete.csv")); nh = nh[nh.educ.isin([1,2,3,4,5])]
nh = nh.drop(columns=[c for c in ["amd_obj_any","amd_obj_late","amd_obj_any_RE"] if c in nh.columns])
m = nh.merge(opx, on="SEQN", how="left").merge(whq, on="SEQN", how="left").merge(bpq[["SEQN","htn_self"]], on="SEQN", how="left")
m["CMM3_meas"] = m["diabetes"] + m["hypertension"] + m["obesity"]
m["CMM3_selfHTN"] = m["diabetes"] + m["htn_self"] + m["obesity"]
m["CMM3_selfALL"] = m["diabetes"] + m["htn_self"] + m["obese_self"]   # fully self-reported, matches NHIS
m["age45"] = (m["age"] >= 45).astype(int)
print("n=", len(m))
print("objective AMD (worse eye): any non-missing", int(m.amd_obj_any.notna().sum()), "any cases", int(m.amd_obj_any.sum()), "late cases", int(m.amd_obj_late.sum()))
print("old right-eye approach cases (for audit):", int(m.amd_obj_any_RE.sum()))
print("self-reported BMI non-missing", int(m.bmi_self.notna().sum()), " self-obese %", round(100*m.obese_self.mean(),1), " measured obese %", round(100*m.obesity.mean(),1))
print("CMM3_selfALL non-missing", int(m.CMM3_selfALL.notna().sum()))
m.to_csv(os.path.join(OUT, "nhanes_analysis_r5.csv"), index=False)
print("wrote nhanes_analysis_r5.csv")
