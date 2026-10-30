# Supplementary Audit — code-level verification of reported statistics

This file documents that every quantity in `眼科2.md` (Tables 1–4, Figures 1–3) is produced by the scripts in `analysis/scripts/` from the individual-level public-use files in `data/raw/`, and records the design/variance diagnostics requested in review. Estimator for all associations: `survey::svyglm` (Taylor linearisation), R 4.3.3, survey 4.4.8. Result files: `analysis/output/sr_*.csv`, `sr2_*.csv`.

## Supplementary Tables S1–S4

**Supplementary Table S1. Survey design and variance diagnostics.**

| Cohort | n | Strata | PSUs | Design df | df used for CI | Taylor SE | JKn SE | Ratio |
|---|---|---|---|---|---|---|---|---|
| NHANES | 5,366 | 31 | 62 | 31 | 22 | 0.0708 | 0.0709 | 1.00 |
| NHIS | 17,613 | 52 | 659 | 607 | 598 | 0.0298 | 0.0299 | 1.00 |

SE is for the CMM3–cataract coefficient. A hand-coded delete-one-PSU variant gave 0.080 (ratio 1.13) and is reported only as a conservative approximation.

**Supplementary Table S2. Meta-analysis inputs and results for CMM3.**

| Outcome | logOR NHANES (SE) | logOR NHIS (SE) | FE weight NHANES/NHIS | FE OR (95% CI) | RE-HK OR (95% CI) | τ² | Q (p) | I² |
|---|---|---|---|---|---|---|---|---|
| Cataract surgery | 0.1875 (0.0708) | 0.2967 (0.0298) | 15% / 85% | 1.32 (1.25–1.40) | 1.30 (0.51–3.27) | 0.0030 | 2.02 (0.16) | 51% |
| AMD | 0.0514 (0.1100) | 0.2081 (0.0480) | 16% / 84% | 1.20 (1.10–1.31) | 1.17 (0.36–3.87) | 0.0051 | 1.70 (0.19) | 41% |

**Supplementary Table S3. Stratified bootstrap diagnostics for the descriptive CRP analysis.**

| Quantity | Estimate | Bootstrap SE | 95% CI | MC SE (lower) | MC SE (upper) | Converged replicates |
|---|---|---|---|---|---|---|
| Indirect difference, cataract surgery | +0.011 | 0.0111 | −0.011 to 0.032 | 0.0010 | 0.0008 | 1000/1000 |
| Indirect difference, AMD | −0.020 | 0.0158 | −0.053 to 0.009 | 0.0018 | 0.0010 | 1000/1000 |
| a-path (CMM4 → log2-CRP) | 0.525 | 0.0297 | 0.468 to 0.585 | 0.0024 | 0.0028 | 1000/1000 |
| Total-effect log-OR (validation) | 0.109 | 0.0565 | 0.004 to 0.221 | 0.0040 | 0.0050 | 1000/1000 |

The total-effect bootstrap SE (0.057) is below the Taylor design SE (0.071), indicating mild anti-conservatism of the cluster bootstrap; inflating the indirect SE by 1.25 still leaves both intervals spanning zero.

**Supplementary Table S4. Missingness, selection and inverse-probability-of-completion diagnostics (NHANES).**

| Step | n |
|---|---|
| Examined, aged ≥ 40 years | 7,081 |
| Valid cataract-surgery data | 7,074 |
| Valid AMD data | 7,005 |
| Complete CMM4 | 5,773 |
| Complete CRP | 5,770 |
| Complete covariates | 5,366 |

| IPC analysis | Complete-case OR (95% CI) | IPC-weighted OR (95% CI) |
|---|---|---|
| Cataract surgery, primary covariates | 1.21 (1.04–1.40) | 1.20 (1.04–1.39) |
| Cataract surgery, income-adjusted | 1.19 (1.02–1.39) | — |
| AMD, primary covariates | 1.05 (0.84–1.32) | 1.07 (0.86–1.34) |

Inverse-probability ratio: median 1.19, 95th percentile 1.52, maximum 2.02; effective sample size 3,309 (unweighted) versus 3,404 (weighted).

## S-A. Script → output map

| Script | Produces |
|---|---|
| `00_discover_nhanes.py` | variable inventory from `.XPT` files |
| `01_build_dataset.py` | `nhanes_analytic_age40_complete.csv`, `nhanes_full_age40_with_missing.csv`, `nhanes_build_summary.json` |
| `06_build_nhis.py` | `nhis_analytic.csv` |
| `10_svyglm_all.R` | `sr_design.csv`, `sr_cohort_cmm3_svyglm.csv`, `sr_cmm3_vs_cmm4.csv`, `sr_age_strata.csv`, `sr_age_interaction.csv`, `sr_progressive.csv`, `sr_dose_response.csv`, `sr_sensitivity_svyglm.csv` |
| `13_audit_round2.R` | `sr2_design_diagnostics.csv`, `sr2_cohort_cmm3_raw.csv`, `sr2_meta_audit.csv`, `sr2_jackknife_vs_taylor.csv`, `sr2_holm_family1_cmm4.csv`, `sr2_cmm3_dys_joint.csv`, `sr2_exposure_shift.csv`, `sr2_stepwise_exclusion.csv`, `sr2_missingness.csv`, `sr2_completers_comparison.csv`, `sr2_ipc_weighting.csv`, `sr2_dose_trend_vs_nonlinear.csv`, `sr2_age_modification_cross_survey.csv` |
| `14_boot_mediation_v2.R` | `sr2_mediation_boot.csv`, `sr2_mediation_boot_summary.csv` |
| `12_figures_svyglm.py` | `figures/Fig1_nhanes_forest.*`, `Fig2_cross_survey.*`, `Fig3_mediation.*` (read the CSVs above) |

## S-B. Survey design verification

- **NHANES**: 2005–2006 + 2007–2008 pooled with `WTMEC2YR / 2` (÷ number of 2-year cycles). Strata = `SDMVSTRA`, PSU = `SDMVPSU`, `nest = TRUE`. Analytic sample: **31 strata, 62 stratum:PSU pairs, every stratum exactly 2 PSUs, 0 singleton strata**. `SDMVPSU` codes (1, 2) repeat across strata → nesting required. All weights > 0.
- **NHIS 2023**: full-sample adult weight `WTFA_A`; strata = `PSTRAT`, PSU = `PPSU` — NCHS's pseudo-strata / pseudo-PSU provided for variance estimation. Analytic sample: **52 strata, 659 stratum:PSU pairs, minimum 3 PSUs/stratum, 0 singleton strata**. `PPSU` codes repeat across strata (80 of 108 raw codes) → `nest = TRUE` required and applied. Weights 1,792–32,759, none ≤ 0.
- **CI degrees of freedom** (design df minus model rank): NHANES 31 − 9 = **22**; NHIS 607 − 9 = **598** (`sr2_cohort_cmm3_raw.csv`, column `df_used_for_CI`). t-multipliers therefore ≈ 2.07 (NHANES) and 1.96 (NHIS).
- **Sub-sample handling**: CRP is measured on the full examined sample → MEC weight applies directly. Retinal photography (≥ 40 examinees) is analysed by `subset()` of the design (not re-weighted as a fresh sample); the resulting 35-case objective-AMD estimate is reported as descriptive/underpowered. 【Confirm at submission that no dedicated imaging-subsample weight is required for NHANES 2005–2008 retinal photography.】

## S-C. Variance-estimator concordance (identical model + sample)

`sr2_jackknife_vs_taylor.csv` — model = cataract ~ CMM3 + ns(age,4) + male + educ3 + smoke; jackknife = delete-one stratum:PSU pair, weights renormalised per replicate:

| Cohort | Taylor SE | Jackknife SE | Ratio | Replicates used |
|---|---|---|---|---|
| NHANES | 0.0708 | 0.0800 | 1.13 | 62/62 |
| NHIS | 0.0298 | 0.0299 | 1.005 | 659/659 |

## S-D. Meta-analysis audit (raw log-OR / SE, not CI back-derivation)

`sr2_meta_audit.csv`:

| Outcome | logOR(NHANES) | SE | logOR(NHIS) | SE | FE weight NHANES/NHIS | FE OR (95% CI) | RE-HK OR (95% CI) | τ² | Q (df=1), p | I² | t used |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Cataract | 0.1875 | 0.0708 | 0.2967 | 0.0298 | 15% / 85% | 1.32 (1.25–1.40) | 1.30 (0.51–3.27) | 0.0030 | 2.02, p=0.16 | 51% | 12.71 |
| AMD | 0.0514 | 0.1100 | 0.2081 | 0.0480 | 16% / 84% | 1.20 (1.10–1.31) | 1.17 (0.36–3.87) | 0.0051 | 1.70, p=0.19 | 41% | 12.71 |

Per-cohort estimates are always presented before the pooled estimate in the text; the fixed-effect summary is flagged as NHIS-dominated.

## S-E. Holm families (regenerated with `p.adjust(method="holm")`)

- **Family 1** (`sr2_holm_family1_cmm4.csv`): NHANES CMM4 → cataract raw p = 0.073, Holm p = 0.121; AMD raw p = 0.061, Holm p = 0.121.
- **Family 2** (`sr2_cohort_cmm3_raw.csv`): CMM3 → NHIS cataract Holm 3.9e-21; NHIS AMD Holm 5.2e-05; NHANES cataract raw 0.0147 → Holm 0.029; NHANES AMD Holm 0.64.
- Everything else (age strata/interaction, dose–response, decomposition, sensitivity, IPC) is exploratory with unadjusted p-values.

## S-F. Mediation bootstrap (stratified) diagnostics

`14_boot_mediation_v2.R`, `sr2_mediation_boot_summary.csv`:
- Resampling: PSUs drawn **within each stratum** with replacement (preserves stratification); each replicate refits the weighted point models (identical coefficients to `svyglm`).
- **Convergence: 1000/1000 replicates converged for all five models (0 failures, 0 separations).**
- Indirect (total − direct, log-odds), percentile 95% CI, with Monte-Carlo SE of each endpoint: cataract **+0.011 (−0.011 to 0.032)**, MC SE ≤ 0.001; AMD **−0.020 (−0.053 to 0.009)**, MC SE ≤ 0.002. a-path 0.525 (0.468–0.585).
- **Validation of the resampling variance** (distinct from the point-estimate identity): bootstrap SE of the total-effect log-OR = 0.057 vs Taylor design SE = 0.071 (ratio 0.80) — the two-PSU-per-stratum cluster bootstrap is mildly anti-conservative. Inflating each indirect SE by 1/0.80 = 1.25 still leaves both intervals spanning zero (cataract −0.016 to 0.039; AMD −0.059 to 0.019), so the "no mediation" conclusion is robust to this.
- The difference-of-coefficients is reported as **descriptive, not an identified causal indirect effect** (non-collapsibility; unmeasured mediator–outcome confounding; single time point).

## S-G. Exposure-definition audit

- `sr2_cmm3_dys_joint.csv` — dyslipidaemia prevalence within CMM3 = 0/1/2/3: 34% / 45% / 53% / 47% (n per cell reported). Dyslipidaemia is a common co-occurrence at every burden level.
- `sr2_exposure_shift.csv` — cataract: CMM3 1.19, CMM4 1.12, dyslipidaemia-own-term 0.94 (0.74–1.19); AMD: CMM3 1.05, CMM4 1.15, dyslipidaemia 1.54 (1.15–2.06). Estimates reported with CIs; the dilution account is framed as plausible, not proven.

## S-H. Missing-data / selection audit

- `sr2_stepwise_exclusion.csv`: 7,081 → valid cataract 7,074 → +valid AMD 7,005 → +complete CMM4 5,773 → +complete CRP 5,770 → +complete covariates 5,366.
- `sr2_missingness.csv`: CMM 17.8%, CRP 9.2%, income 8.1%, outcomes < 1%.
- `sr2_completers_comparison.csv`: completers vs non-completers (age 59.6/63.9; diabetes 15.4/20.4; cataract 13.0/21.0; AMD 3.7/5.7; education and income higher among completers).
- `sr2_ipc_weighting.csv` + S-F completion-model diagnostics (income-adjusted variant): cataract 1.19 → 1.19 (1.03–1.37); like-for-like primary-covariate comparison in S-K (1.21 → 1.20); AMD 1.05 → 1.07; weights max ratio 2.02 (none > 3), ESS 3,309 → 3,404, 99th-pct trimming unchanged. Residual unmeasured selection bias cannot be excluded.

## S-J. Round-4 additions (measurement harmonisation + like-for-like IPC)

Additional scripts:
| Script | Produces |
|---|---|
| `15_build_bpq_selfhtn.py` | merges NHANES `BPQ020` into the analytic file → `nhanes_analysis_selfhtn.csv` (adds `CMM3_selfhtn`) |
| `16_selfhtn_sensitivity.R` | `sr2_selfhtn_sensitivity.csv` |
| `17_ipc_likeforlike.R` | `sr2_ipc_likeforlike.csv` |

**Like-for-like IPC (Section 3.6).** Using the *same covariates as the primary CMM3 model* (no income term): cataract complete-case OR 1.21 (1.05–1.39) → IPC-weighted 1.20 (1.04–1.38), p = 0.017; AMD 1.05 (0.85–1.31) → 1.07 (0.87–1.33). The 1.19 quoted in §3.3/§3.5 is the **income-adjusted** model on the same cohort (complete-case 1.19 → IPC-weighted 1.19); the 1.21-vs-1.19 difference is therefore covariate adjustment only, not a different sample or exposure. `sr2_ipc_likeforlike.csv`.

**Measurement harmonisation (Section 3.7).** NHANES `BPQ020` (ever told high blood pressure) was used to rebuild CMM3 with self-reported hypertension, matching NHIS, while retaining measured BMI. Self-reported hypertension prevalence 45.5% vs measured 25.7% (as expected with treated/controlled patients). CMM3–cataract OR: measured 1.21 (1.04–1.40) vs self-reported 1.24 (1.05–1.46) — similar under the alternative definition, suggesting this measurement difference is unlikely to fully explain the association. Obesity ascertainment (measured BMI vs NHIS self-reported category) could not be harmonised because NHANES 2005–2008 lacks self-reported height/weight; residual heterogeneity remains for that component. `sr2_selfhtn_sensitivity.csv`. The raw `BPQ_D.xpt` / `BPQ_E.xpt` files were retrieved from the CDC public-use server and are stored in `data/raw/nhanes/`.

## S-I. Reproducibility

Fixed bootstrap seed 20261029; BLAS threads pinned to 1 before library load; ASCII staging directory used for R to avoid non-ASCII path issues. A versioned public repository with these scripts, a `README` (seed, package versions, script→number map) and the frozen environment is to be deposited — 【repository URL / version DOI: AUTHOR-SPECIFIC】.

## S-K. Final-consistency verification (round 5)

**Exact sample flow** (`18_final_consistency.R`, `sr2_sample_flow.csv`): 7,081 examined → 7,074 valid cataract surgery (−7) → 7,005 valid AMD (−69) → 5,773 complete CMM4 (−1,232) → **5,770 complete CRP (−3)** → 5,366 complete covariates (−404). The 5,773 → 5,366 difference is therefore −3 −404 = −407, i.e. 5,366 exactly; the apparent 5,773 − 404 arithmetic in the previous draft omitted the three-person CRP step, which is now shown explicitly in Section 3.1.

**Weighted Table 1** (`sr2_table1_weighted.csv`): survey-weighted means/proportions for both cohorts; NHANES age 56.4 (55.6–57.2), NHIS 62.7 (62.5–62.9); cataract surgery 9.8% vs 18.6%; diabetes 11.2% vs 15.9%; hypertension 21.6% vs 49.5%; obesity 36.3% vs 34.9%; dyslipidaemia 42.3% (NHANES only); AMD 3.0% vs 3.3%; education distribution 18.2/25.8/56.0 vs 14.5/38.8/46.7.

**Variance-estimator validation** (`sr2_jackknife_validation.csv`): Taylor SE 0.0708; survey-package **JKn** replicate design 0.0709 (ratio 1.00); hand-coded delete-one-PSU 0.080 (ratio 1.13, conservative approximation). The JKn/Taylor agreement is the valid concordance check; the hand-coded variant is not used as the reference.

**Like-for-like IPC, t-based** (`sr2_ipc_likeforlike_t.csv`): cataract complete-case 1.21 (1.04–1.40) → IPC 1.20 (1.04–1.39); AMD 1.05 (0.84–1.32) → 1.07 (0.86–1.34); income-adjusted 1.19 (1.02–1.39).

## S-L. Round-6 corrections (age eligibility, retinal grading, self-reported BMI)

**NHIS age eligibility.** The 2023 NHIS eye-outcome questions (`VIMCSURG_A`, `VIMMDEV_A`) are administered only from age 45. Verified in the raw file: of 2,201 respondents aged 40–44, **none** has a non-missing eye outcome; the minimum age with an outcome is 45, and the analytic NHIS sample (n = 17,613) has minimum age **45** (0 aged < 45). The models are therefore unaffected, but the manuscript now states NHIS ≥ 45 years, and the age strata are relabelled (NHANES 40–59; NHIS 45–59).

**Objective AMD — grading corrected.** The previous build used `OPDDARM > 1` (right-eye ARM severity), which is not a person-level "any AMD" variable and gave only 35 cases. Rebuilt with the **worse-eye** variables: any ARM = `OPDUARM ≥ 1` (3-level severity) and late ARM = `OPDULARM = 1`. Among the analysed participants, **4,722 had gradable photographs, with 353 any-ARM and 45 late-ARM cases** (`19_build_r5_corrected.py`, `20_r5_corrected_models.R`, `sr2_r5_corrected_models.csv`). CMM3–any-ARM OR **1.02 (0.83–1.26)**, p = 0.85; CMM3–late-ARM OR 0.93 (0.58–1.50), p = 0.77. The objective analysis is thus a better-powered null, not the earlier underpowered 35-case estimate.

**Self-reported BMI available.** NHANES collects self-reported height and weight in the Weight History files (`WHQ_D/E.xpt`; `WHD010` inches, `WHD020` pounds; valid ranges applied, 9999/7777 etc. set missing). 5,223 analysed participants had self-reported BMI; self-reported obesity was 33.0% versus 37.3% measured. A **fully self-reported CMM3** (self-reported diabetes, self-reported hypertension and self-reported BMI ≥ 30), matching NHIS ascertainment, gave a cataract–surgery OR of **1.26 (1.07–1.49)**. The earlier statement that NHANES lacked self-reported height/weight is corrected.

**Common age range.** Restricting NHANES to ≥ 45 years (n = 4,653) gave CMM3–cataract OR **1.20 (1.03–1.39)**, essentially unchanged from the full ≥ 40 sample (1.21, 1.04–1.40).

**Degrees of freedom audit.** `degf(design)` = 31 (NHANES) and 607 (NHIS); `df.residual(svyglm)` = 22 and 598 (design df minus model rank), and `confint`/`summary` use the same residual df. The hand-coded and JKn jackknife concordance checks are in S-C/S-K.

**Reference verification (PMIDs).** All 24 references were checked via PubMed/Crossref. Because references are numbered in order of first appearance (Nature style), the list was renumbered after the author's own papers were added. PMIDs by current number: 1 29032195; 2 28242111; 3 15708105; 4 21144031; 5 34333624; 6 25060855; 7 24587127; 8 40646243; 9 25474403; 10 19933187; 11 14871913; 12 25207946; 13 10192648; 14 18064739; 15 22788053; 16 24775908; 17 1951829; 18 12939780; 19 24548571; 20 26653405; 21 42039168; 22 41383630; 23 42039201; 24 41327731.

