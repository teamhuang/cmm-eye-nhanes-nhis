"""Build NHIS 2023 harmonized analytic set for cross-cohort replication (Topic 3).
Harmonized definitions (common to NHANES & NHIS):
  exposure CMM3 = diabetes + hypertension + obesity   (dyslipidemia dropped: no NHIS item)
  cataract = cataract SURGERY (VIMCSURG_A==1)  [matches NHANES VIQ071]
  amd      = told macular degeneration (VIMMDEV_A==1)  [matches NHANES VIQ310]
  covars: age, male, 3-level education, ever-smoke; weight WTFA_A; design PPSU/PSTRAT.
Value recode: 1=yes, 2=no, else (7/8/9/97/99) -> NaN. Restrict age>=40, complete on needed vars.
"""
import pandas as pd, numpy as np, os
RAW = r"H:\眼科论文\眼科2\data\raw\nhis"
OUT = r"C:\ops2\output"
os.makedirs(OUT, exist_ok=True)
usecols = ['AGEP_A','SEX_A','EDUCP_A','SMKEV_A','DIBEV_A','HYPEV_A','BMICAT_A',
           'VIMCSURG_A','VIMMDEV_A','VIMCAEV_A','VIMGLEV_A','VIMDREV_A','WTFA_A','PPSU','PSTRAT']
df = pd.read_csv(os.path.join(RAW,'adult23.csv'), usecols=usecols)
def yesno(s):
    s = pd.to_numeric(s, errors='coerce')
    return pd.Series(np.where(s==1,1.0, np.where(s==2,0.0, np.nan)), index=s.index)

d = pd.DataFrame()
d['age']   = pd.to_numeric(df['AGEP_A'], errors='coerce')
d['male']  = (pd.to_numeric(df['SEX_A'],errors='coerce')==1).astype(float)
d['smoke'] = yesno(df['SMKEV_A'])
d['diabetes']= yesno(df['DIBEV_A'])
d['hypertension']= yesno(df['HYPEV_A'])
bmi = pd.to_numeric(df['BMICAT_A'], errors='coerce')
d['obesity'] = np.where(bmi==4,1.0, np.where(bmi.isin([1,2,3]),0.0,np.nan))
d['cataract']= yesno(df['VIMCSURG_A'])
d['amd']     = yesno(df['VIMMDEV_A'])
e = pd.to_numeric(df['EDUCP_A'], errors='coerce')
d['educ'] = np.where(e.isin([1,2,3]),0.0, np.where(e.isin([4,5]),1.0, np.where(e.isin([6,7,8,9,10]),2.0,np.nan)))
d['wt']=pd.to_numeric(df['WTFA_A'],errors='coerce'); d['psu']=df['PPSU']; d['strata']=df['PSTRAT']

d['CMM'] = d[['diabetes','hypertension','obesity']].sum(axis=1, skipna=False)
need=['age','male','smoke','CMM','cataract','amd','educ','wt']
mask=(d['age']>=40) & d[need].notna().all(axis=1)
da=d[mask].copy()
print('NHIS rows>=40 complete:', len(da))
print('prevalences:', {'cataract_surg':round(da.cataract.mean(),4),'amd':round(da.amd.mean(),4),
      'diabetes':round(da.diabetes.mean(),4),'hypertension':round(da.hypertension.mean(),4),
      'obesity':round(da.obesity.mean(),4)})
print('CMM dist:', da.CMM.value_counts().sort_index().to_dict())
print('PSUs:', da.psu.nunique(), 'strata:', da.strata.nunique())
da.to_csv(os.path.join(OUT,'nhis_analytic.csv'), index=False)
print('wrote', os.path.join(OUT,'nhis_analytic.csv'))
