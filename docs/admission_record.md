# Admission Record — Topic 3 (RUN-FREEZE, write before any model fit)

**Study title (working)**: Cardiometabolic multimorbidity burden and eye disease in older adults — Cataract and Age-Related Macular Degeneration — with C-Reactive Protein as mediator: a NHANES-anchored, internationally replicated secondary analysis (NHANES US; ELSA UK; SHARE Europe).

**Freeze date**: 2026-10-28. This file records the pre-registered claim, metric, direction, expected effect ± tolerance, data, analytic sample, and modeling choices. After running, changes to any of these are reported as **a NEW analysis** with the original retained; do not silently re-tune.

## Primary claim (NHANES derivation cohort)
Higher cardiometabolic multimorbidity (CMM) burden is associated with higher odds of (a) cataract / cataract surgery and (b) age-related macular degeneration (AMD), and this association is **partially mediated** by systemic low-grade inflammation (log-CRP).

## Analysis unit
Person. Survey weights = NHANES MEC 2-year weights combined across 2005-06 & 2007-08 with the standard 1/2 scaling; design SDMVPSU / SDMVSTRA. Nested eye-level measures are aggregated to worst-better-eye (AMD) or person-level (cataract surgery, self-report) so no pseudo-replication.

## Data (already verified as free direct-download from CDC)
NHANES 2005-2008 (cycles _D and _E):
- DEMO_{D,E}: demographics + weights.
- VIQ_{D,E}: vision questionnaire (self-reported eye conditions: cataract / surgery / glaucoma / AMD; vision trouble).
- OPXRET_{D,E}: retinal photofundus grading (objective AMD severity OPDDARM etc.; DR grades available but NOT used as primary given user-excluded DR framing).
- CRP_{D,E}: C-reactive protein (mg/dL).
- BPX_{D,E}: blood pressure.
- TCHOL/TRIGLY/HDL_{D,E} (+ LDL_D for 2005-06): lipids.
- BMX_{D,E}: anthropometrics (BMI).
- DIQ_{D,E}: diabetes dx.
- SMQ_{D,E}: smoking.
- GLU_{D,E}: fasting glucose (secondary check).

Analytic cohort: NHANES participants **aged ≥40** (age at which cataract/AMD incidence becomes clinically meaningful and objective grading is validated). Primary = **complete-case** with valid CRP + valid eye outcome + valid weights.

## Exposure (frozen definition)
**CMM count = number of cardiometabolic conditions present (0-4)**:
1. Diabetes: DIQ010==1 OR (fasting glucose ≥126 OR HbA1c ≥6.5 — GLU file if usable).
2. Hypertension: mean SBP ≥140 OR mean DBP ≥90 (BPX, 3rd/4th/5th/6th reading mean).
3. Dyslipidemia: TC ≥240 mg/dL OR HDL <40 (M)/<50 (W) OR TG ≥200 OR (calculated/measured LDL ≥160).
4. Obesity: BMI ≥30.
Rationale: use **measured** criteria (not self-report) because NHANES suppressed "told-high-cholesterol" questionnaire for 2005-08.

## Mediator
log2(CRP + 0.01), CRP in mg/dL from CRP_D/E (`CRPP`). log2 required by manual (intensity data → log scale; CRP is right-skewed).

## Outcomes (both reported; Holm across them)
- **Y1 = cataract / cataract surgery** (VIQ file): "ever had cataract operation" (primary), "diagnosed cataract" if present (sensitivity).
- **Y2 = AMD**: self-reported "ever told had macular degeneration" (VIQ, primary harmonized outcome for cross-cooling with ELSA/SHARE) + **objective AMD grading** (OPXRET `OPDDARM`, any-ARM vs late-ARM) as **sensitivity** in NHANES only (ELSA/SHARE have no imaging).

## Direction & expected effects (with tolerance for PASS band)
Expected positive:
- CMM count → cataract surgery OR per-unit ≈ 1.15-1.80 (PASS band: point estimate >1, CI excludes 1, OR Holm-adj p<0.05).
- CMM count → AMD OR ≈ 1.20-2.50.
- CMM count → log2 CRP: positive, +0.05 to +0.25 per unit.
- Proportion mediated via CRP: 5-20%.
Honest null is also acceptable and will be reported as such. Do NOT flip endpoint, dataset, or CMM definition to force significance.

## Covariates (adjusted set)
Age (natural spline, 4 knots), sex, race/ethnicity (5-level NHANES standard), education (<HS/HS/Some-col/College), PIR quintile, current smoking, alcohol (never/former/current via ALQ surrogate or SMQ alcohol where applicable), physical-activity binary (PAQ if available else omit), eye-care access for cataract-surgery sensitivity (VIQ050 "eye exam in past 12 mo"). BMI NOT in adjustment for CMM-count models (obesity is inside CMM); BMI included in "excluding obesity from CMM" sensitivity.

## Statistical models
1. **Exposure-outcome (Y|X)**: survey-weighted logistic (`LogisticRegression`/statsmodels GLM Binomial with design-variance via `svyglm`-style approximation, or R `survey::svyglm` — see runtime notes). OR + 95% CI.
2. **Mediator model (M|X)**: survey-weighted OLS, log2 CRP.
3. **Outcome|mediator (Y|X,M)**: survey-weighted logistic with log2 CRP added.
4. **Mediation estimation** (VanderWeele / Imai product-of-coefficients on appropriate scales): A (X→M), B (M→Y adjusted), indirect ≈ A×B (log-OR scale), natural direct effect = total − indirect. 2,000-replicate BCa bootstrap (seed = 20261028) with the survey weights resampled at PSU level.
5. **Multiple comparisons**: Holm across the 2 outcomes for the primary exposure effect and separately across 2 outcomes for the mediated effect.
6. **Sensitivity**: (a) objective AMD grading instead of self-report; (b) excluding statins/anti-hypertensives (to reduce "on-treatment" attenuation); (c) alternate CMM definitions (syndrome z-score; binary "any"); (d) restricting to non-diabetes subsample for AMD (avoid DR/AMD confound); (e) add eye-care access adjustment.
7. **Cross-cohort replication** (ELSA/SHARE, when files supplied by user): refit Y|X (no mediation — CRP not available in some) on harmonized self-reported cataract / AMD outcomes with the same covariates; compare effect direction & magnitude via a fixed-effect and RE meta across the three cohorts (heterogeneity I^2; REML+HK given small k=3). Mediation replication attempted in ELSA (has CRP).

## Reporting
- STROBE for observational part.
- Report raw p to 2 sig-figs alongside Holm-adjusted; never `p=0.000`.
- Effect sizes + CI only; no causal verbs ("causes/drives/leads to"); use "associated with".
- Fig 1: cohort flow. Fig 2: forest plot of CMM→Y across outcomes and cohorts. Fig 3: mediation path diagram with A/B/indirect. Fig 4 (sensitivity): objective AMD grading.

## Non-negotiables (integrity)
- Do NOT fabricate ELSA/SHARE numbers: those cohorts run only after the user supplies downloaded files.
- Do NOT use DR/age/DME framing that the user has excluded.
- Do NOT pool NHANES objective grading with ELSA/SHARE self-report — they are different measurement constructs; report as sensitivity or separately.
- Do NOT overwrite this file with a "results-fit" version. Findings go into `analysis/output/findings.md`.

## Verdict thresholds
- **PASS** if primary CMM→Y1 and CMM→Y2 have point estimates with 95% CI excluding 1 in expected direction, Holm-adj p<0.05, and CRP mediation produces a positive indirect effect with BCa CI not covering 0.
- **PARTIAL PASS** if some combinations produce honest nulls — report each with CI and interpret without spin.
- **FAIL** only if the pipeline is broken or numbers cannot be sourced.
