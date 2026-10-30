#!/usr/bin/env Rscript
# Topic 3 — pre-registered sensitivity analyses (NHANES).
# (a) OBJECTIVE retinal AMD grading (stronger phenotype) instead of self-report.
# (b) CMM dose-response (categorical 0/1/2/3-4).
# (c) cataract surgery adjusted for eye-care access (VIQ056 = eye exam in past 12 mo) — access confounding.
# (d) CMM excluding obesity (obesity is the component closest to CRP -> avoids mediator overlap).
# Same design-weighted point estimates + delete-a-PSU jackknife variance. Report honestly.
suppressMessages({library(splines)})
OUT <- "C:/ops2/output"
prep <- function(dat){
  dat<-dat[dat$educ %in% 1:5,]
  dat$educf<-factor(dat$educ)
  ab<-as.data.frame(ns(dat$age,df=4)); colnames(ab)<-paste0("ag",seq_len(ncol(ab))); dat<-cbind(dat,ab)
  dat$grp<-paste(dat$SDMVSTRA,dat$SDMVPSU); dat$w1<-dat$wt/mean(dat$wt); dat
}
cfb <- function(f,dat,fam,nm) unname(coef(glm(f,data=dat,weights=w1,family=fam,control=list(maxit=200)))[nm])
jk_stat <- function(f,dat,fam,nm){
  pt<-cfb(f,dat,fam,nm); g<-unique(dat$grp)
  j<-sapply(g,function(gr){dd<-dat[dat$grp!=gr,]; dd$w1<-dd$wt/mean(dd$wt); tryCatch(cfb(f,dd,fam,nm),error=function(e)NA)})
  j<-j[is.finite(j)]; se<-sqrt((length(j)-1)/length(j)*sum((j-mean(j))^2))
  c(est=pt,se=se,lo=pt-1.96*se,hi=pt+1.96*se,p=2*(1-pnorm(abs(pt/se))))
}
srows <- list()
addrow <- function(analysis, outcome, term, r, n){
  srows[[length(srows)+1]] <<- data.frame(analysis=analysis, outcome=outcome, term=term,
     OR=exp(r["est"]), OR_lo=exp(r["lo"]), OR_hi=exp(r["hi"]), p=r["p"], n=n, check.names=FALSE)
}

cov<-"ag1+ag2+ag3+ag4 + male + educf + pir + ever_smoke"

# ---- MAIN complete set ----
dm <- prep(read.csv(file.path(OUT,"nhanes_analytic_age40_complete.csv")))

# (b) dose-response, cataract & AMD (factor CMM, ref=0)
dm$cmm_cat<-factor(ifelse(dm$CMM>=3,"3+",as.character(dm$CMM)), levels=c("0","1","2","3+"))
for (y in c("cataract_surg","amd_self")){
  f<-as.formula(paste0(y," ~ cmm_cat + ",cov))
  for (lv in c("1","2","3+")){
    nm<-paste0("cmm_cat",lv); r<-jk_stat(f,dm,binomial(),nm)
    addrow("dose_response", y, paste0("CMM=",lv), r, nrow(dm))
    cat(sprintf("[dose] %-14s CMM=%-2s OR=%.3f (%.3f-%.3f) p=%.3g\n", y, lv, exp(r["est"]),exp(r["lo"]),exp(r["hi"]),r["p"]))
  }
}

# (d) CMM excluding obesity
dm$CMM_noob <- dm$diabetes + dm$hypertension + dm$dyslipidemia
for (y in c("cataract_surg","amd_self")){
  f<-as.formula(paste0(y," ~ CMM_noob + ",cov)); r<-jk_stat(f,dm,binomial(),"CMM_noob")
  addrow("CMM_no_obesity", y, "per-unit", r, nrow(dm))
  cat(sprintf("[no-obesity] %-14s OR=%.3f (%.3f-%.3f) p=%.3g\n", y, exp(r["est"]),exp(r["lo"]),exp(r["hi"]),r["p"]))
}

# (c) cataract with eye-care access adjustment (VIQ056)
if("VIQ056" %in% names(dm)){
  dm$exam<- (dm$VIQ056==1)*1
  dacc<-dm[!is.na(dm$exam),]
  f<-as.formula(paste0("cataract_surg ~ CMM + exam + ",cov)); r<-jk_stat(f,dacc,binomial(),"CMM")
  addrow("access_adjusted", "cataract_surg", "per-unit", r, nrow(dacc))
  cat(sprintf("[cataract+access] OR=%.3f (%.3f-%.3f) p=%.3g (n=%d)\n", exp(r["est"]),exp(r["lo"]),exp(r["hi"]),r["p"],nrow(dacc)))
}

# ---- (a) OBJECTIVE AMD from the with-missing file (retinal subsample) ----
dfo <- read.csv(file.path(OUT,"nhanes_full_age40_with_missing.csv"))
need<-c("amd_obj_any","amd_obj_late","CMM","log2crp","age","male","educ","pir","ever_smoke","wt","SDMVSTRA","SDMVPSU")
dfo<-dfo[complete.cases(dfo[,need]),]; dfo<-dfo[dfo$amd_obj_any %in% 0:1,]
dfo<-prep(dfo)
cat("\n[objective AMD] n=",nrow(dfo)," cases=",sum(dfo$amd_obj_any),"\n")
f<-as.formula(paste0("amd_obj_any ~ CMM + ",cov)); r<-jk_stat(f,dfo,binomial(),"CMM")
addrow("objective_AMD_any", "amd_obj_any", "per-unit", r, nrow(dfo))
cat(sprintf("[objective AMD any] OR=%.3f (%.3f-%.3f) p=%.3g\n", exp(r["est"]),exp(r["lo"]),exp(r["hi"]),r["p"]))
# late AMD
if(any(dfo$amd_obj_late %in% 0:1,na.rm=TRUE) && sum(dfo$amd_obj_late)>0){
  dfol<-dfo[dfo$amd_obj_late %in% 0:1,]; f<-as.formula(paste0("amd_obj_late ~ CMM + ",cov)); rl<-jk_stat(f,dfol,binomial(),"CMM")
  addrow("objective_AMD_late", "amd_obj_late", "per-unit", rl, nrow(dfol))
  cat(sprintf("[objective AMD late] OR=%.3f (%.3f-%.3f) p=%.3g (n=%d cases=%d)\n", exp(rl["est"]),exp(rl["lo"]),exp(rl["hi"]),rl["p"],nrow(dfol),sum(dfol$amd_obj_late)))
} else {
  cat("[objective AMD late] 0 cases -> uninformative, not estimated\n")
}

out<-do.call(rbind, srows)
write.csv(out, file.path(OUT,"sensitivity_results.csv"), row.names=FALSE)
cat("\nwrote sensitivity_results.csv rows=", nrow(out), "\n")
