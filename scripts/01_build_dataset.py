"""Build NHANES analytic dataset for Topic 3 (CMM -> cataract / AMD, CRP mediator).

Freeze rules applied (admission_record.md):
- Analysis unit = person.
- Cycles 2005-06 (SDDSRVYR=4) and 2007-08 (SDDSRVYR=5) combined -> WTMEC2YR / 2.
- Analytic set: age >= 40, valid cataract surgery + valid self-reported AMD outcome.
- log2(CRP + 0.01) because CRP is right-skewed intensity (manual: log scale).
- CMM = count of measured cardiometabolic conditions (diabetes / hypertension / dyslipidemia / obesity).
Do NOT change definitions post hoc. Missingness reported honestly.
"""
import pyreadstat, pandas as pd, numpy as np, os, sys, json

RAW = r"H:\眼科论文\眼科2\data\raw\nhanes"
OUT = r"H:\眼科论文\眼科2\analysis\output"
os.makedirs(OUT, exist_ok=True)

def read(f):
    path = os.path.join(RAW, f)
    if not os.path.exists(path):
        print(f"MISSING FILE: {f}"); return None
    df, _ = pyreadstat.read_xport(path)
    df.columns = [c.upper() for c in df.columns]
    return df

def cyc(suffix):
    """load a variable's two cycles."""
    a = read(f"{suffix}_D.xpt"); b = read(f"{suffix}_E.xpt")
    if a is None or b is None: return None
    return pd.concat([a, b], ignore_index=True)

# --- demographics + weights ---
demo = cyc("DEMO")
demo = demo[['SEQN','SDDSRVYR','RIDSTATR','RIAGENDR','RIDAGEYR','RIDRETH1','DMDEDUC2','INDFMPIR','WTMEC2YR','SDMVPSU','SDMVSTRA']]

# --- vision questionnaire: cataract surgery (VIQ071), glaucoma (VIQ090), AMD (VIQ310), self-rated vision (VIQ031), trouble seeing (VIQ041), eye exam (VIQ056) ---
viq  = cyc("VIQ")[['SEQN','VIQ017','VIQ031','VIQ041','VIQ056','VIQ061','VIQ071','VIQ090','VIQ310']]

# --- CRP ---
crp  = cyc("CRP")[['SEQN','LBXCRP']]

# --- diabetes ---
diq  = cyc("DIQ")[['SEQN','DIQ010']]

# --- BP: use readings 2,3,4 mean (NHANES standard) ---
bpx  = cyc("BPX")[['SEQN','BPXSY2','BPXSY3','BPXSY4','BPXDI2','BPXDI3','BPXDI4']]
bpx['SBP_mean'] = bpx[['BPXSY2','BPXSY3','BPXSY4']].mean(axis=1, skipna=True)
bpx['DBP_mean'] = bpx[['BPXDI2','BPXDI3','BPXDI4']].mean(axis=1, skipna=True)
bpx = bpx[['SEQN','SBP_mean','DBP_mean']].drop_duplicates('SEQN')

# --- BMI ---
bmx  = cyc("BMX")[['SEQN','BMXBMI']].drop_duplicates('SEQN')

# --- lipids ---
tc   = cyc("TCHOL")[['SEQN','LBXTC']].drop_duplicates('SEQN')
hdl  = cyc("HDL")[['SEQN','LBDHDD']].drop_duplicates('SEQN')
trg  = cyc("TRIGLY")[['SEQN','LBXTR','LBDLDL']].drop_duplicates('SEQN')

# --- smoking: SMQ020 = ever smoked >=100 cigarettes (asked of ALL adults); SMQ040 only for ever-smokers ---
smq  = cyc("SMQ")[['SEQN','SMQ020']].drop_duplicates('SEQN')

# --- retinal grading (for sensitivity only; merged as extra cols) ---
opx_d = read("OPXRET_D.xpt"); opx_e = read("OPXRET_E.xpt")
opx = pd.concat([opx_d, opx_e], ignore_index=True) if (opx_d is not None and opx_e is not None) else None
if opx is not None:
    keep = ['SEQN','OPDDARM','OPDDLARM','OPDDARMA','OPDDEXU','OPDDAGA']
    keep = [c for c in keep if c in opx.columns]
    opx = opx[keep].drop_duplicates('SEQN')

# --- merge ---
d = demo.copy()
for right in [viq, crp, diq, bpx, bmx, tc, hdl, trg, smq, opx]:
    if right is not None:
        d = d.merge(right, on='SEQN', how='left')

print("merged rows (all ages):", len(d))

# --- recode 7/9 (missing) → NaN on binary questionnaire vars ---
for v in ['VIQ017','VIQ031','VIQ041','VIQ056','VIQ061','VIQ071','VIQ090','VIQ310','DIQ010','SMQ020']:
    if v in d:
        d[v] = pd.to_numeric(d[v], errors='coerce').replace({7.0:np.nan, 9.0:np.nan, 77.0:np.nan, 99.0:np.nan})

# numeric coercion for measured variables
for v in ['LBXCRP','BMXBMI','LBXTC','LBDHDD','LBXTR','LBDLDL','SBP_mean','DBP_mean','RIDAGEYR','WTMEC2YR','SDMVPSU','SDMVSTRA','INDFMPIR','DMDEDUC2','RIDRETH1','RIAGENDR','SDDSRVYR']:
    if v in d:
        d[v] = pd.to_numeric(d[v], errors='coerce')

# --- define indicators ---
d['diabetes']     = (d['DIQ010']==1).astype('float'); d.loc[d['DIQ010'].isna(),'diabetes']=np.nan
d['hypertension'] = (pd.to_numeric(d['SBP_mean'],errors='coerce').ge(140) | pd.to_numeric(d['DBP_mean'],errors='coerce').ge(90)).astype('float')
d.loc[d['SBP_mean'].isna() & d['DBP_mean'].isna(), 'hypertension'] = np.nan

def dyslip(row):
    tc = row.get('LBXTC'); hdl = row.get('LBDHDD'); sex = row.get('RIAGENDR')
    if pd.isna(tc) and pd.isna(hdl): return np.nan
    hi_tc = pd.notna(tc) and tc >= 240
    lo_hdl = pd.notna(hdl) and ((sex==1 and hdl<40) or (sex==2 and hdl<50))
    return float(hi_tc or lo_hdl)
d['dyslipidemia'] = d.apply(dyslip, axis=1)
d['obesity'] = pd.to_numeric(d['BMXBMI'], errors='coerce').ge(30).astype('float')
d.loc[d['BMXBMI'].isna(),'obesity']=np.nan

# CMM count
ind = d[['diabetes','hypertension','dyslipidemia','obesity']]
d['CMM'] = ind.sum(axis=1, skipna=False)  # any NaN -> NaN
d['CMM_valid'] = ind.notna().all(axis=1)

# outcomes (binary)
d['cataract_surg'] = (d['VIQ071']==1).astype('float'); d.loc[d['VIQ071'].isna(),'cataract_surg']=np.nan
d['amd_self']      = (d['VIQ310']==1).astype('float'); d.loc[d['VIQ310'].isna(),'amd_self']=np.nan
d['glaucoma_self'] = (d['VIQ090']==1).astype('float'); d.loc[d['VIQ090'].isna(),'glaucoma_self']=np.nan

# mediator
d['log2crp'] = np.log2(pd.to_numeric(d['LBXCRP'], errors='coerce') + 0.01)

# covariates
d['age'] = d['RIDAGEYR']
d['male'] = (d['RIAGENDR']==1).astype('float')
d['educ'] = pd.Categorical(d['DMDEDUC2']).codes  # 1-5
d['educ'] = pd.to_numeric(d['DMDEDUC2'], errors='coerce')
d['pir'] = d['INDFMPIR']
d['ever_smoke'] = (d['SMQ020']==1).astype('float'); d.loc[d['SMQ020'].isna(),'ever_smoke']=np.nan

# combined weight (2 cycles)
d['wt'] = d['WTMEC2YR']/2.0

# objective AMD (sensitivity): OPDDARM levels — >1 = AMD; OPDDLARM>1 = late ARM
d['amd_obj_any']  = pd.to_numeric(d.get('OPDDARM'),errors='coerce').gt(1).astype('float')
d.loc[pd.to_numeric(d.get('OPDDARM'),errors='coerce').isna(),'amd_obj_any']=np.nan
d['amd_obj_late'] = pd.to_numeric(d.get('OPDDLARM'),errors='coerce').gt(1).astype('float')
d.loc[pd.to_numeric(d.get('OPDDLARM'),errors='coerce').isna(),'amd_obj_late']=np.nan

# --- age >= 40 restriction ---
d40 = d[(d['age']>=40)].copy()
print("rows age>=40:", len(d40))

# --- primary analytic: valid cataract & AMD & CMM & log2crp & weight ---
need = ['cataract_surg','amd_self','CMM','log2crp','wt','age','male','educ','pir','ever_smoke']
mask_all = d40[need].notna().all(axis=1)
dA = d40[mask_all].copy()
print("rows full-complete-case:", len(dA))

# missingness in age>=40
miss = d40[[c for c in d40.columns if c in set(need)|{'SBP_mean','DBP_mean','LBXTC','LBDHDD','BMXBMI','DIQ010','SMQ040'}]].isna().mean().round(3).sort_values(ascending=False)
print("\nMissingness (age>=40):"); print(miss.to_string())

# save
csv_main = os.path.join(OUT, 'nhanes_analytic_age40_complete.csv')
csv_all  = os.path.join(OUT, 'nhanes_full_age40_with_missing.csv')
dA.to_csv(csv_main, index=False)
d40.to_csv(csv_all, index=False)

# summary json
summ = {
  'rows_merged_all_ages': int(len(d)),
  'rows_age40': int(len(d40)),
  'rows_complete_case_primary': int(len(dA)),
  'cycle_counts_40': {int(k):int(v) for k,v in d40['SDDSRVYR'].value_counts().to_dict().items()},
  'cycle_counts_cc': {int(k):int(v) for k,v in dA['SDDSRVYR'].value_counts().to_dict().items()},
  'prevalence_cc': {
     'cataract_surg': float(dA['cataract_surg'].mean()),
     'amd_self': float(dA['amd_self'].mean()),
     'glaucoma_self': float(dA['glaucoma_self'].mean()),
     'diabetes': float(dA['diabetes'].mean()),
     'hypertension': float(dA['hypertension'].mean()),
     'dyslipidemia': float(dA['dyslipidemia'].mean()),
     'obesity': float(dA['obesity'].mean()),
  },
  'CMM_dist_cc': {int(k):int(v) for k,v in dA['CMM'].value_counts().sort_index().to_dict().items()},
  'crp_median_cc': float(dA['LBXCRP'].median()),
  'crp_iqr_cc': [float(dA['LBXCRP'].quantile(.25)), float(dA['LBXCRP'].quantile(.75))],
}
print("\nSUMMARY:"); print(json.dumps(summ, indent=2))
with open(os.path.join(OUT, 'nhanes_build_summary.json'),'w') as f:
    json.dump(summ,f,indent=2)
print("Saved:", csv_main)
print("Saved:", csv_all)
