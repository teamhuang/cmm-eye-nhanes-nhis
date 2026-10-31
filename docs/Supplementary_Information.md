# Supplementary Information

**Cardiometabolic multimorbidity and age-related eye disease: cross-sectional analyses of two US national surveys**

## Supplementary Tables

**Supplementary Table S1. Survey design and variance diagnostics.**

| Cohort | n | Strata | PSUs | Design df | df used for CI | Taylor SE | JKn SE | Ratio |
|---|---|---|---|---|---|---|---|---|
| NHANES | 5,366 | 31 | 62 | 31 | 22 | 0.0708 | 0.0709 | 1.00 |
| NHIS | 17,613 | 52 | 659 | 607 | 598 | 0.0298 | 0.0299 | 1.00 |

SE, standard error for the CMM3–cataract coefficient; JKn, delete-one-PSU jackknife. A hand-coded delete-one-PSU variant gave 0.080 (ratio 1.13) and is reported only as a conservative approximation.

**Supplementary Table S2. Meta-analysis inputs and results for CMM3.**

| Outcome | logOR NHANES (SE) | logOR NHIS (SE) | Fixed-effect weight, NHANES/NHIS | Fixed-effect OR (95% CI) | Random-effects OR (95% CI) | τ² | Q (p) | I² |
|---|---|---|---|---|---|---|---|---|
| Cataract surgery | 0.1875 (0.0708) | 0.2967 (0.0298) | 15% / 85% | 1.32 (1.25–1.40) | 1.30 (0.51–3.27) | 0.0030 | 2.02 (0.16) | 51% |
| AMD | 0.0514 (0.1100) | 0.2081 (0.0480) | 16% / 84% | 1.20 (1.10–1.31) | 1.17 (0.36–3.87) | 0.0051 | 1.70 (0.19) | 41% |

logOR, log odds ratio; SE, standard error; OR, odds ratio; CI, confidence interval; τ², between-study variance; I², heterogeneity.

**Supplementary Table S3. Stratified bootstrap diagnostics for the descriptive CRP analysis.**

| Quantity | Estimate | Bootstrap SE | 95% CI | Monte Carlo SE, lower | Monte Carlo SE, upper | Converged replicates |
|---|---|---|---|---|---|---|
| Difference, cataract surgery | +0.011 | 0.0111 | −0.011 to 0.032 | 0.0010 | 0.0008 | 1000/1000 |
| Difference, AMD | −0.020 | 0.0158 | −0.053 to 0.009 | 0.0018 | 0.0010 | 1000/1000 |
| CMM4 → log2-CRP (a-path) | 0.525 | 0.0297 | 0.468 to 0.585 | 0.0024 | 0.0028 | 1000/1000 |
| Total-effect log-OR (validation) | 0.109 | 0.0565 | 0.004 to 0.221 | 0.0040 | 0.0050 | 1000/1000 |

CRP, C-reactive protein; SE, standard error; CI, confidence interval. The bootstrap standard error of the total-effect log-OR (0.057) was below the Taylor design standard error (0.071); inflating the difference standard errors by 1.25 still leaves both intervals spanning zero.

**Supplementary Table S4. Missingness, selection and inverse-probability-of-completion diagnostics (NHANES).**

| Analysis step | n |
|---|---|
| Examined, aged ≥ 40 years | 7,081 |
| Valid cataract-surgery data | 7,074 |
| Valid AMD data | 7,005 |
| Complete CMM4 | 5,773 |
| Complete CRP | 5,770 |
| Complete covariates | 5,366 |

| Analysis | Complete-case OR (95% CI) | IPC-weighted OR (95% CI) |
|---|---|---|
| Cataract surgery, primary covariates | 1.21 (1.04–1.40) | 1.20 (1.04–1.39) |
| Cataract surgery, income-adjusted | 1.19 (1.02–1.39) | — |
| AMD, primary covariates | 1.05 (0.84–1.32) | 1.07 (0.86–1.34) |

IPC, inverse probability of completion; OR, odds ratio; CI, confidence interval. Inverse-probability ratio: median 1.19, 95th percentile 1.52, maximum 2.02; effective sample size 3,309 (unweighted) versus 3,404 (weighted).

## Supplementary Figures

**Supplementary Figure S1. Participant flow in NHANES.** Numbers show participants remaining after each complete-case requirement.

![Figure S1](figures/FigS1_flow.png){w=8,h=8.9}

**Supplementary Figure S2. Bootstrap distributions of the descriptive CRP difference.** Histograms show the stratified PSU-bootstrap distributions of the total-minus-CRP-adjusted difference for cataract surgery and AMD; the dashed line marks zero.

![Figure S2](figures/FigS2_bootstrap.png){w=8,h=5.3}

**Supplementary Figure S3. Prevalence of dyslipidaemia within each CMM3 level in NHANES.** Bars show the percentage with dyslipidaemia among participants with 0, 1, 2 or 3 of the other cardiometabolic conditions.

![Figure S3](figures/FigS3_cmm3_dys.png){w=8,h=5.3}
