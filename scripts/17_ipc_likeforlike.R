#!/usr/bin/env Rscript
# Section 3.6 like-for-like IPC comparison: complete-case vs IPC-weighted, using the SAME covariates
# as the primary CMM3 model (no income). Also reports the income-adjusted variant to explain why the
# paper quotes both 1.21 (primary) and 1.19 (income-adjusted) for the same analytic cohort.
suppressMessages({library(survey); library(splines)})
OUT <- "C:/ops2/output"
if (!file.exists(file.path(OUT, "nhanes_analytic_age40_complete.csv"))) OUT <- "H:/眼科论文/眼科2/analysis/output"
mkSpline <- function(d){ b <- as.data.frame(ns(d$age, df=4)); names(b) <- paste0("ag", 1:4); cbind(d, b) }
nh <- read.csv(file.path(OUT, "nhanes_analytic_age40_complete.csv")); nh <- nh[nh$educ %in% 1:5,]
nh$cataract <- nh$cataract_surg; nh$amd <- nh$amd_self; nh$smoke <- nh$ever_smoke
nh$CMM3 <- nh$diabetes + nh$hypertension + nh$obesity; nh$CMM4 <- nh$CMM
nh$educ3 <- factor(ifelse(nh$educ <= 2, 0, ifelse(nh$educ == 3, 1, 2)))
nh <- nh[!is.na(nh$CMM3) & !is.na(nh$CMM4),]; nh <- mkSpline(nh)
full <- read.csv(file.path(OUT, "nhanes_full_age40_with_missing.csv")); full <- full[!is.na(full$age),]; full <- mkSpline(full)
full$included <- full$SEQN %in% nh$SEQN
desfull <- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~WTMEC2YR, data=full, nest=TRUE)
ipm <- svyglm(included ~ ag1+ag2+ag3+ag4 + RIAGENDR + RIDRETH1 + SDDSRVYR, design=desfull, family=quasibinomial())
full$p_inc <- predict(ipm, type="response")
nhm <- merge(nh, full[,c("SEQN","p_inc")], by="SEQN"); nhm <- nhm[is.finite(nhm$p_inc) & nhm$p_inc > 0,]; nhm$ipc_w <- nhm$wt / nhm$p_inc
des_cc <- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~wt, data=nhm, nest=TRUE)
des_ipc <- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~ipc_w, data=nhm, nest=TRUE)
row <- function(des, y, rhs, tag){ m <- svyglm(as.formula(paste(y, "~", rhs)), des, family=quasibinomial())
  b <- unname(coef(m)["CMM3"]); s <- unname(sqrt(vcov(m)["CMM3","CMM3"]))
  data.frame(analysis=tag, OR=exp(b), lo=exp(b-1.96*s), hi=exp(b+1.96*s), p=unname(summary(m)$coefficients["CMM3","Pr(>|t|)"])) }
P <- "CMM3 + ag1+ag2+ag3+ag4 + male + educ3 + smoke"; Pinc <- paste(P, "+ pir")
out <- rbind(
  row(des_cc, "cataract", P,    "cataract_complete_case_PRIMARY_no_income"),
  row(des_ipc,"cataract", P,    "cataract_IPC_weighted_PRIMARY_no_income"),
  row(des_cc, "cataract", Pinc, "cataract_complete_case_WITH_income"),
  row(des_ipc,"cataract", Pinc, "cataract_IPC_weighted_WITH_income"),
  row(des_cc, "amd",      P,    "amd_complete_case_PRIMARY_no_income"),
  row(des_ipc,"amd",      P,    "amd_IPC_weighted_PRIMARY_no_income"))
print(out, row.names=FALSE); write.csv(out, file.path(OUT, "sr2_ipc_likeforlike.csv"), row.names=FALSE)
