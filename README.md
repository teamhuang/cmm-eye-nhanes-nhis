# Cardiometabolic multimorbidity and age-related eye disease (NHANES 2005–2008; NHIS 2023)

Analysis code and result files for the manuscript:

> **Cardiometabolic multimorbidity and age-related eye disease: cross-sectional analyses of two US national surveys**
> Jianhuang Huang, Jing Wang

This repository contains every script needed to build the analytic files, fit the models and regenerate the tables, figures and audit values reported in the paper. It does **not** contain individual-level survey data; those are freely available from CDC/NCHS (see below).

## Repository contents

```
scripts/   R and Python scripts (run order below)
outputs/   result CSVs (sr_*.csv, sr2_*.csv) used for the manuscript
figures/   Figure 1–3 as 600 dpi PNG and vector PDF
docs/      pre-specified analysis plan, STROBE checklist, supplementary audit
```

## Data (not included)

Download the public-use files from CDC/NCHS and place them under `data/raw/`:

- **NHANES 2005–2006 and 2007–2008** component files (DEMO, VIQ, OPXRET, CRP, DIQ, BPX, BMX, TCHOL, HDL, TRIGLY, SMQ, BPQ, WHQ):
  `https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/<year>/DataFiles/<FILE>.xpt`
  (e.g. `.../2005/DataFiles/BPQ_D.xpt`; the older `/nchs/nhanes/<cycle>/` path may return an HTML block page).
- **NHIS 2023 adult file**:
  `https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Datasets/NHIS/2023/adult23csv.zip`

## Environment

- R 4.3.3 with `survey` 4.4.8 and `splines`.
- Python 3.10 with `pandas`, `numpy`, `pyreadstat`, `matplotlib`.
- Random seed for the mediation bootstrap: **20261029**.
- The analysis was run on Windows. Because the project path contains non-ASCII characters, R wrote to an ASCII staging directory; each script defines `OUT` (output directory) and `RAW` (data directory) at the top. Set these two constants to your local paths before running.

## Run order

```
01_build_dataset.py      -> NHANES analytic file (n = 5,366)
06_build_nhis.py         -> NHIS 2023 harmonised file (n = 17,613)
10_svyglm_all.R          -> primary, cross-survey, age, dose, sensitivity models
13_audit_round2.R        -> design diagnostics, meta audit, Holm, decomposition, IPC
14_boot_mediation_v2.R   -> stratified PSU bootstrap for the descriptive CRP analysis
12_figures_svyglm.py     -> Figures 1–3
15_build_bpq_selfhtn.py  -> self-reported-hypertension variable
16_selfhtn_sensitivity.R -> measurement-harmonisation analysis
17_ipc_likeforlike.R     -> like-for-like inverse-probability-of-completion analysis
18_final_consistency.R   -> weighted Table 1, sample flow, JKn validation
19_build_r5_corrected.py -> corrected worse-eye objective AMD; self-reported BMI
20_r5_corrected_models.R -> corrected objective AMD and self-reported CMM3 models
```

## Script → manuscript item map

A detailed map is in `docs/Supplementary_audit.md` (Section S-A). In brief:

- Table 1: `18_final_consistency.R` (`sr2_table1_weighted.csv`)
- Table 2 and Figure 2: `10_svyglm_all.R`, `14_boot_mediation_v2.R`
- Table 3 and Figure 3: `10_svyglm_all.R`, `13_audit_round2.R`, `12_figures_svyglm.py`
- Table 4: `10_svyglm_all.R`
- Figure 1: `10_svyglm_all.R`, `19/20_r5_corrected_models`
- Design/variance diagnostics: `13_audit_round2.R`, `18_final_consistency.R`

## Notes on key definitions

- All models are survey-weighted (`survey::svyglm`, Taylor linearisation); weights are rescaled to mean 1.
- The primary cross-survey adjustment set is age spline, sex, education and smoking; income is a NHANES-only sensitivity.
- CMM3 = diabetes + hypertension + obesity (harmonised across surveys); CMM4 adds dyslipidaemia in NHANES.
- Objective AMD uses worse-eye grading: any ARM = `OPDUARM >= 1`; late ARM = `OPDULARM = 1`.
- NHIS eye-outcome questions are asked only from age 45; the NHANES sample is aged >= 40.

## Result files

`outputs/` contains the exact CSVs used for the manuscript (including `sr2_meta_audit.csv`, `sr2_jackknife_validation.csv`, `sr2_ipc_likeforlike_t.csv`, `sr2_mediation_boot_summary.csv`, `sr2_r5_corrected_models.csv`). These allow the reported numbers to be checked without re-running the models.

## Citation

If you use this code, please cite the manuscript and this repository:

> Huang J, Wang J. Cardiometabolic multimorbidity and age-related eye disease: cross-sectional analyses of two US national surveys. GitHub repository, version v1.0.3. https://github.com/teamhuang/cmm-eye-nhanes-nhis

A versioned archive DOI (Zenodo) will be added if required by the journal.

## Contact

Jianhuang Huang — teamhuang@sina.com
Jing Wang — 13273530801@163.com
