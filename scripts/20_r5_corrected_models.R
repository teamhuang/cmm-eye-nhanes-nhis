#!/usr/bin/env Rscript
# Round-5 corrected analyses (design rebuilt per subset): objective AMD worse-eye (any/late);
# fully self-reported CMM3; self-HTN; primary; common-age >=45 sensitivity; design-df audit.
suppressMessages({library(survey); library(splines)})
OUT <- "C:/ops2/output"
if (!file.exists(file.path(OUT,"nhanes_analysis_r5.csv"))) OUT <- "H:/眼科论文/眼科2/analysis/output"
d <- read.csv(file.path(OUT,"nhanes_analysis_r5.csv"))
d$cataract <- d$cataract_surg; d$smoke <- d$ever_smoke
d$educ3 <- factor(ifelse(d$educ<=2,0, ifelse(d$educ==3,1,2)))
bb <- as.data.frame(ns(d$age,df=4)); names(bb) <- paste0("ag",1:4); d <- cbind(d,bb)
COV <- "ag1+ag2+ag3+ag4 + male + educ3 + smoke"
fitOR <- function(dat, y, expo){
  des <- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~wt, data=dat, nest=TRUE)
  m <- svyglm(as.formula(paste(y,"~",expo,"+",COV)), des, family=quasibinomial())
  b <- unname(coef(m)[expo]); s <- unname(sqrt(vcov(m)[expo,expo])); tc <- qt(0.975, df.residual(m))
  data.frame(analysis=paste0(y,"~",expo), n=nrow(dat), events=sum(dat[[y]],na.rm=TRUE),
    OR=exp(b), lo=exp(b-tc*s), hi=exp(b+tc*s), p=unname(summary(m)$coefficients[expo,"Pr(>|t|)"]),
    df_resid=df.residual(m), df_design=degf(des)) }
o <- rbind(
  fitOR(d[!is.na(d$amd_obj_any),],  "amd_obj_any",  "CMM3_meas"),
  fitOR(d[!is.na(d$amd_obj_late),], "amd_obj_late", "CMM3_meas"),
  fitOR(d[!is.na(d$CMM3_selfALL),], "cataract",     "CMM3_selfALL"),
  fitOR(d[!is.na(d$CMM3_selfHTN),], "cataract",     "CMM3_selfHTN"),
  fitOR(d,                          "cataract",     "CMM3_meas"),
  fitOR(d[d$age45==1,],             "cataract",     "CMM3_meas"))
print(o, row.names=FALSE); write.csv(o, file.path(OUT,"sr2_r5_corrected_models.csv"), row.names=FALSE)
cat("\nobjective AMD any cases:", sum(d$amd_obj_any,na.rm=TRUE), " late:", sum(d$amd_obj_late,na.rm=TRUE), "\n")
cat("self-BMI available:", sum(!is.na(d$bmi_self)), " self-obese %:", round(100*mean(d$obese_self,na.rm=TRUE),1),
    " measured-obese %:", round(100*mean(d$obesity,na.rm=TRUE),1), "\n")
cat("NHANES >=45 n:", sum(d$age45), "\n")
