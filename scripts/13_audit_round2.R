#!/usr/bin/env Rscript
# ROUND 2 STATISTICAL AUDIT. Recomputes every reported quantity from raw svyglm coefficients/SEs
# (no CI back-derivation), verifies the survey design, tests jackknife-vs-Taylor on the IDENTICAL
# model+sample, defines the two Holm families explicitly, and adds the reviewer-requested
# distributions / selection sensitivity. Outputs sr2_*.csv. Estimator = survey::svyglm (Taylor).
Sys.setenv(OMP_NUM_THREADS=1, MKL_NUM_THREADS=1, OPENBLAS_NUM_THREADS=1, NUMEXPR_NUM_THREADS=1)
suppressMessages({library(survey); library(splines)})
OUT <- "C:/ops2/output"

# ---------- data ----------
nh <- read.csv(file.path(OUT,"nhanes_analytic_age40_complete.csv")); nh <- nh[nh$educ %in% 1:5,]
nh$cataract<-nh$cataract_surg; nh$amd<-nh$amd_self; nh$smoke<-nh$ever_smoke
nh$CMM3 <- nh$diabetes + nh$hypertension + nh$obesity; nh$CMM4 <- nh$CMM
nh$educ3 <- factor(ifelse(nh$educ<=2,0, ifelse(nh$educ==3,1,2)))
nh <- nh[!is.na(nh$CMM3) & !is.na(nh$CMM4),]
nh$psu <- nh$SDMVPSU; nh$strata <- nh$SDMVSTRA
nh$agegrp <- factor(cut(nh$age, c(39,59,69,200), labels=c("40-59","60-69","70+")), levels=c("40-59","60-69","70+"))
des_nh <- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~wt, data=nh, nest=TRUE)

ni <- read.csv(file.path(OUT,"nhis_analytic.csv"))
ni$educ3 <- factor(ni$educ); ni$CMM3 <- ni$CMM
ni$agegrp <- factor(cut(ni$age, c(39,59,69,200), labels=c("40-59","60-69","70+")), levels=c("40-59","60-69","70+"))
des_ni <- svydesign(id=~psu, strata=~strata, weights=~wt, data=ni, nest=TRUE)
COV <- "ns(age, df=4) + male + educ3 + smoke"

# ---------- 1. design diagnostics ----------
desdiag <- function(dat, label){
  psuper <- tapply(dat$psu, dat$strata, function(x) length(unique(x)))
  data.frame(cohort=label, n=nrow(dat), n_strata=length(unique(dat$strata)),
             n_psu=length(unique(paste(dat$strata,dat$psu))), design_df=length(unique(paste(dat$strata,dat$psu)))-length(unique(dat$strata)),
             strata_with_1_psu=sum(psuper==1), strata_with_ge2_psu=sum(psuper>=2), min_psu_per_strata=min(psuper))
}
dd_nh <- desdiag(nh, "NHANES"); dd_ni <- desdiag(ni, "NHIS")
write.csv(rbind(dd_nh, dd_ni), file.path(OUT,"sr2_design_diagnostics.csv"), row.names=FALSE)

# ---------- 2. per-cohort CMM3 from RAW coef/SE ----------
getbse <- function(des, y, spec, term){
  m <- svyglm(as.formula(paste(y,"~",spec)), des, family=quasibinomial())
  c(logOR=unname(coef(m)[term]), se=unname(sqrt(diag(vcov(m)))[term]),
    df=df.residual(m), p=unname(summary(m)$coefficients[term,"Pr(>|t|)"]))
}
rows <- list()
for (cohort in c("NHANES","NHIS")){
  des <- if(cohort=="NHANES") des_nh else des_ni
  for (y in c("cataract","amd")){
    r <- getbse(des, y, paste("CMM3 +",COV), "CMM3")
    rows[[length(rows)+1]] <- data.frame(cohort=cohort, outcome=y, logOR=r["logOR"], se=r["se"],
      df_used_for_CI=r["df"], OR=exp(r["logOR"]), p_raw=r["p"]) } }
cb <- do.call(rbind, rows)
cb$p_holm4 <- p.adjust(cb$p_raw, method="holm")
write.csv(cb, file.path(OUT,"sr2_cohort_cmm3_raw.csv"), row.names=FALSE)

# ---------- 3. meta from TRUE logOR/SE with full audit ----------
meta_audit <- function(sub, label){
  b <- sub$logOR; se <- sub$se; v <- se^2; k <- length(b)
  wf <- 1/v; mu_f <- sum(wf*b)/sum(wf); se_f <- sqrt(1/sum(wf))
  Q <- sum(wf*(b-mu_f)^2); df <- k-1
  tau2 <- max(0,(Q-df)/(sum(wf)-sum(wf^2)/sum(wf)))
  wr <- 1/(v+tau2); mu_r <- sum(wr*b)/sum(wr); q2 <- max(1,Q/df); se_hk <- sqrt(q2/sum(wr))
  tcrit <- qt(0.975, max(1,df))
  data.frame(outcome=label, k=k,
    logOR_NHANES=b[1], se_NHANES=se[1], logOR_NHIS=b[2], se_NHIS=se[2],
    weight_fe_NHANES=wf[1]/sum(wf), weight_fe_NHIS=wf[2]/sum(wf),
    weight_re_NHANES=wr[1]/sum(wr), weight_re_NHIS=wr[2]/sum(wr),
    OR_fe=exp(mu_f), fe_lo=exp(mu_f-1.96*se_f), fe_hi=exp(mu_f+1.96*se_f), p_fe=2*pnorm(-abs(mu_f/se_f)),
    OR_re_hk=exp(mu_r), re_lo=exp(mu_r-tcrit*se_hk), re_hi=exp(mu_r+tcrit*se_hk), p_re_hk=2*pt(-abs(mu_r/se_hk),df=max(1,df)),
    tau2=tau2, Q=Q, df_het=df, p_het=pchisq(Q,df,lower.tail=FALSE), I2=max(0,(Q-df)/Q)*100, t_crit_used=tcrit)
}
ma <- rbind(meta_audit(cb[cb$outcome=="cataract",],"cataract"), meta_audit(cb[cb$outcome=="amd",],"amd"))
write.csv(ma, file.path(OUT,"sr2_meta_audit.csv"), row.names=FALSE)

# ---------- 4. jackknife vs Taylor on the IDENTICAL model + sample ----------
jk_vs_taylor <- function(des, dat, y){
  m <- svyglm(as.formula(paste(y,"~ CMM3 +",COV)), des, family=quasibinomial())
  taylor_se <- unname(sqrt(diag(vcov(m)))["CMM3"])
  d <- dat; d$grp <- paste(d$strata, d$psu); d$w1 <- d$wt/mean(d$wt)
  f <- as.formula(paste(y,"~ CMM3 + ns(age,df=4) + male + educ3 + smoke"))
  g <- unique(d$grp)
  j <- sapply(g, function(gr){ dd <- d[d$grp!=gr,]; dd$w1 <- dd$wt/mean(dd$wt)
       tryCatch(unname(coef(glm(f,data=dd,weights=w1,family=binomial()))["CMM3"]), error=function(e) NA) })
  j <- j[is.finite(j)]; R <- length(j)
  jk_se <- sqrt((R-1)/R*sum((j-mean(j))^2))
  data.frame(model="CMM3+COV", taylor_se=taylor_se, jackknife_se=jk_se,
             ratio_jk_over_taylor=jk_se/taylor_se, jackknife_reps_used=R, psu_total=length(g))
}
jt <- rbind(cbind(cohort="NHANES", jk_vs_taylor(des_nh, nh, "cataract")),
            cbind(cohort="NHIS",   jk_vs_taylor(des_ni, ni, "cataract")))
write.csv(jt, file.path(OUT,"sr2_jackknife_vs_taylor.csv"), row.names=FALSE)

# ---------- 5. Holm family 1: NHANES CMM4 primary (2 tests) ----------
r_cat <- getbse(des_nh, "cataract", paste("CMM4 +",COV,"+ pir"), "CMM4")
r_amd <- getbse(des_nh, "amd",      paste("CMM4 +",COV,"+ pir"), "CMM4")
fam1 <- data.frame(outcome=c("cataract","amd"), logOR=c(r_cat["logOR"],r_amd["logOR"]),
  se=c(r_cat["se"],r_amd["se"]), df=c(r_cat["df"],r_amd["df"]), p_raw=c(r_cat["p"],r_amd["p"]))
fam1$OR <- exp(fam1$logOR); fam1$lo <- exp(fam1$logOR - qt(0.975,fam1$df)*fam1$se); fam1$hi <- exp(fam1$logOR + qt(0.975,fam1$df)*fam1$se)
fam1$p_holm2 <- p.adjust(fam1$p_raw, method="holm")
write.csv(fam1[,c("outcome","OR","lo","hi","p_raw","p_holm2")], file.path(OUT,"sr2_holm_family1_cmm4.csv"), row.names=FALSE)

# ---------- 6. CMM3 x dyslipidaemia joint distribution + exposure shift ----------
jd <- as.data.frame.table(table(CMM3=nh$CMM3, dyslipidemia=nh$dyslipidemia), responseName="n")
jd <- jd[jd$n>0,]
jd$pct_within_CMM3 <- with(jd, 100*n/ave(n, CMM3, FUN=sum))
write.csv(jd, file.path(OUT,"sr2_cmm3_dys_joint.csv"), row.names=FALSE)
shift <- list()
for (y in c("cataract","amd")){
  a <- getbse(des_nh, y, paste("CMM3 +",COV,"+ pir"), "CMM3")
  b <- getbse(des_nh, y, paste("CMM4 +",COV,"+ pir"), "CMM4")
  c <- getbse(des_nh, y, paste("CMM3 + dyslipidemia +",COV,"+ pir"), "dyslipidemia")
  shift[[length(shift)+1]] <- data.frame(outcome=y,
    CMM3_OR=exp(a["logOR"]), CMM3_lo=exp(a["logOR"]-1.96*a["se"]), CMM3_hi=exp(a["logOR"]+1.96*a["se"]),
    CMM4_OR=exp(b["logOR"]), CMM4_lo=exp(b["logOR"]-1.96*b["se"]), CMM4_hi=exp(b["logOR"]+1.96*b["se"]),
    dys_OR=exp(c["logOR"]), dys_lo=exp(c["logOR"]-1.96*c["se"]), dys_hi=exp(c["logOR"]+1.96*c["se"])) }
write.csv(do.call(rbind,shift), file.path(OUT,"sr2_exposure_shift.csv"), row.names=FALSE)

# ---------- 7. missingness + stepwise exclusion + completers vs non-completers ----------
full <- read.csv(file.path(OUT,"nhanes_full_age40_with_missing.csv"))
need <- c("cataract_surg","amd_self","CMM","log2crp","wt","age","male","educ","pir","ever_smoke")
write.csv(data.frame(variable=need, pct_missing=round(100*sapply(full[need], function(x) mean(is.na(x))),2)),
          file.path(OUT,"sr2_missingness.csv"), row.names=FALSE)
ok <- function(...) Reduce(`&`, lapply(list(...), function(e) !is.na(e)))
vcat <- !is.na(full$cataract_surg); vamd <- !is.na(full$amd_self); vCMM <- !is.na(full$CMM)
vcrp <- !is.na(full$log2crp); veduc <- full$educ %in% 1:5; vpir <- !is.na(full$pir); vsmk <- !is.na(full$ever_smoke)
step <- data.frame(step=c("age>=40 examined (full file)","valid cataract surgery","+ valid AMD","+ complete CMM4","+ complete CRP","+ complete covariates (educ/pir/smoke)"),
  n=c(nrow(full), sum(vcat), sum(vcat&vamd), sum(vcat&vamd&vCMM), sum(vcat&vamd&vCMM&vcrp),
      sum(vcat&vamd&vCMM&vcrp&veduc&vpir&vsmk)))
write.csv(step, file.path(OUT,"sr2_stepwise_exclusion.csv"), row.names=FALSE)
full$inc <- full$SEQN %in% nh$SEQN
cmp_vars <- c("age","RIAGENDR","RIDRETH1","DMDEDUC2","INDFMPIR","hypertension","obesity","diabetes","dyslipidemia","cataract_surg","amd_self")
getm <- function(v, grp){ x <- suppressWarnings(as.numeric(full[[v]])); if(v=="RIAGENDR") x <- as.numeric(full[[v]]==1); mean(x[full$inc==grp], na.rm=TRUE) }
comp <- do.call(rbind, lapply(cmp_vars, function(v) data.frame(var=v, completers=round(getm(v,TRUE),3), noncompleters=round(getm(v,FALSE),3))))
write.csv(comp, file.path(OUT,"sr2_completers_comparison.csv"), row.names=FALSE)

# ---------- 8. IPC-weighted sensitivity (selection-probability weighting) ----------
full$included <- full$inc
desfull <- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~WTMEC2YR, data=full, nest=TRUE)
ipm <- svyglm(included ~ ns(age,df=4) + RIAGENDR + RIDRETH1 + SDDSRVYR, design=desfull, family=quasibinomial())
full$p_inc <- predict(ipm, type="response")
nhm <- merge(nh, full[,c("SEQN","p_inc")], by="SEQN"); nhm$ipc_w <- nhm$wt/nhm$p_inc
des_ipc <- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~ipc_w, data=nhm, nest=TRUE)
ipc_rows <- list()
for (y in c("cataract","amd")){
  cc <- getbse(des_ipc, y, paste("CMM3 +",COV,"+ pir"), "CMM3")
  cc2 <- getbse(des_nh, y, paste("CMM3 +",COV,"+ pir"), "CMM3")
  ipc_rows[[length(ipc_rows)+1]] <- data.frame(outcome=y, OR_completecase=exp(cc2["logOR"]),
    OR_ipcw=exp(cc["logOR"]), ipcw_lo=exp(cc["logOR"]-1.96*cc["se"]), ipcw_hi=exp(cc["logOR"]+1.96*cc["se"]), p_ipcw=cc["p"]) }
write.csv(do.call(rbind,ipc_rows), file.path(OUT,"sr2_ipc_weighting.csv"), row.names=FALSE)

# ---------- 9. dose: trend vs nonlinearity; cross-cohort age-modification test ----------
dose <- list()
for (cohort in c("NHANES","NHIS")){
  dat <- if(cohort=="NHANES") nh else ni
  dat$CMMsq <- dat$CMM3^2
  des <- svydesign(id=~psu, strata=~strata, weights=~wt, data=dat, nest=TRUE)
  for (y in c("cataract","amd")){
    lin <- getbse(des, y, paste("CMM3 +",COV), "CMM3")
    mq <- svyglm(as.formula(paste(y,"~ CMM3 + CMMsq +",COV)), des, family=quasibinomial())
    nonlin_p <- unname(summary(mq)$coefficients["CMMsq","Pr(>|t|)"])
    dose[[length(dose)+1]] <- data.frame(cohort=cohort, outcome=y,
      linear_trend_OR=exp(lin["logOR"]), linear_trend_p=lin["p"], p_departure_from_linear=nonlin_p) } }
write.csv(do.call(rbind,dose), file.path(OUT,"sr2_dose_trend_vs_nonlinear.csv"), row.names=FALSE)

age_mod <- function(des, y){
  m <- svyglm(as.formula(paste(y,"~ CMM3*agegrp +",COV)), des, family=quasibinomial())
  nm <- grep("CMM3:agegrp70", names(coef(m)), value=TRUE)
  c(delta = -unname(coef(m)[nm]), se = unname(sqrt(diag(vcov(m)))[nm]))  # logOR(40-59)-logOR(70+)
}
am_nh <- age_mod(des_nh,"cataract"); am_ni <- age_mod(des_ni,"cataract")
z <- (am_nh["delta"]-am_ni["delta"])/sqrt(am_nh["se"]^2+am_ni["se"]^2)
amtab <- data.frame(outcome="cataract", delta_NHANES=am_nh["delta"], se_NHANES=am_nh["se"],
  delta_NHIS=am_ni["delta"], se_NHIS=am_ni["se"], z_cross_survey=z, p_cross_survey=2*pnorm(-abs(z)))
write.csv(amtab, file.path(OUT,"sr2_age_modification_cross_survey.csv"), row.names=FALSE)

cat("ROUND2 AUDIT DONE\n")
cat("\n== design diagnostics ==\n"); print(rbind(dd_nh,dd_ni), row.names=FALSE)
cat("\n== jackknife vs Taylor (identical model/sample) ==\n"); print(jt, row.names=FALSE)
cat("\n== meta audit (cataract) ==\n"); print(ma[1,], row.names=FALSE)
cat("\n== cross-cohort age-modification test ==\n"); print(amtab, row.names=FALSE)
cat("\n== IPC weighting ==\n"); print(do.call(rbind,ipc_rows), row.names=FALSE)
