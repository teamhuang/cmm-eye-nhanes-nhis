#!/usr/bin/env Rscript
# Section 3.7 measurement-harmonisation sensitivity: CMM3 with self-reported physician-diagnosed
# hypertension (BPQ020, matches NHIS) vs measured hypertension. Same covariates as the primary model.
suppressMessages({library(survey); library(splines)})
OUT <- "C:/ops2/output"
if (!file.exists(file.path(OUT, "nhanes_analysis_selfhtn.csv"))) OUT <- "H:/眼科论文/眼科2/analysis/output"
d <- read.csv(file.path(OUT, "nhanes_analysis_selfhtn.csv"))
d <- d[!is.na(d$CMM3_selfhtn) & !is.na(d$cataract_surg),]
d$cataract <- d$cataract_surg; d$smoke <- d$ever_smoke
d$educ3 <- factor(ifelse(d$educ <= 2, 0, ifelse(d$educ == 3, 1, 2)))
d$CMM3_meas <- d$diabetes + d$hypertension + d$obesity
bb <- as.data.frame(ns(d$age, df=4)); names(bb) <- paste0("ag", 1:4); d <- cbind(d, bb)
des <- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~wt, data=d, nest=TRUE)
fit <- function(exp){
  m <- svyglm(as.formula(paste("cataract ~", exp, "+ ag1+ag2+ag3+ag4 + male + educ3 + smoke")), des, family=quasibinomial())
  b <- unname(coef(m)[exp]); s <- unname(sqrt(vcov(m)[exp, exp]))
  c(OR=exp(b), lo=exp(b-qt(0.975,df.residual(m))*s), hi=exp(b+qt(0.975,df.residual(m))*s), p=unname(summary(m)$coefficients[exp,"Pr(>|t|)"])) }
rm_ <- fit("CMM3_meas"); rs <- fit("CMM3_selfhtn")
out <- data.frame(exposure=c("CMM3 measured HTN (primary)", "CMM3 self-reported HTN (harmonised to NHIS)"),
  OR=c(rm_["OR"], rs["OR"]), lo=c(rm_["lo"], rs["lo"]), hi=c(rm_["hi"], rs["hi"]), p=c(rm_["p"], rs["p"]))
print(out, row.names=FALSE)
write.csv(out, file.path(OUT, "sr2_selfhtn_sensitivity.csv"), row.names=FALSE)
