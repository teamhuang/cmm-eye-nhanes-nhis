#!/usr/bin/env Rscript
# Topic 3 — two-cohort replication & meta (NHANES 2005-08 + NHIS 2023), harmonized CMM3.
# Same frozen model: eye ~ CMM3 + ns(age,4) + male + educ3 + smoke; design-weighted point,
#   delete-a-PSU jackknife variance; Holm across outcomes; 2-cohort FE + REML/Hartung-Knapp RE meta.
# Honest scope: NHIS is an independent US sample (same country, different era), all self-report,
#   no CRP -> this replicates the ASSOCIATION only, not mediation, and not internationally.
suppressMessages({library(splines)})
OUT <- "C:/ops2/output"

load_nhanes <- function(){
  d <- read.csv(file.path(OUT,"nhanes_analytic_age40_complete.csv"))
  d <- d[d$educ %in% 1:5,]
  d$cataract <- d$cataract_surg; d$amd <- d$amd_self; d$smoke <- d$ever_smoke
  d$psu <- d$SDMVPSU; d$strata <- d$SDMVSTRA
  d$CMM <- d$diabetes + d$hypertension + d$obesity          # harmonized CMM3
  d$educ3 <- factor(ifelse(d$educ<=2,0, ifelse(d$educ==3,1,2)))
  d[!is.na(d$CMM), c("age","male","educ3","smoke","CMM","cataract","amd","wt","psu","strata")]
}
load_nhis <- function(){
  d <- read.csv(file.path(OUT,"nhis_analytic.csv"))
  d$educ3 <- factor(d$educ); d$male <- d$male
  d[, c("age","male","educ3","smoke","CMM","cataract","amd","wt","psu","strata")]
}
prep <- function(d){ ab<-as.data.frame(ns(d$age,df=4)); colnames(ab)<-paste0("ag",seq_len(ncol(ab))); d<-cbind(d,ab); d$grp<-paste(d$strata,d$psu); d$w1<-d$wt/mean(d$wt); d }
cf <- function(f,dat,fam,nm) unname(coef(glm(f,data=dat,weights=w1,family=fam,control=list(maxit=200)))[nm])
jkfit <- function(f,dat,fam,nm){
  pt<-cf(f,dat,fam,nm); g<-unique(dat$grp)
  j<-sapply(g,function(gr){dd<-dat[dat$grp!=gr,]; dd$w1<-dd$wt/mean(dd$wt); tryCatch(cf(f,dd,fam,nm),error=function(e)NA)})
  j<-j[is.finite(j)]; se<-sqrt((length(j)-1)/length(j)*sum((j-mean(j))^2)); c(est=pt,se=se)
}
cohorts <- list(NHANES=prep(load_nhanes()), NHIS=prep(load_nhis()))

res <- do.call(rbind, lapply(names(cohorts), function(nm){
  d<-cohorts[[nm]]; out<-list()
  for(y in c("cataract","amd")){
    f<-as.formula(paste0(y," ~ CMM + ag1+ag2+ag3+ag4 + male + educ3 + smoke"))
    r<-jkfit(f,d,binomial(),"CMM")
    out[[length(out)+1]]<-data.frame(cohort=nm,outcome=y,n=nrow(d),logOR=r["est"],se=r["se"],
        OR=exp(r["est"]),lo=exp(r["est"]-1.96*r["se"]),hi=exp(r["est"]+1.96*r["se"]),
        p=2*(1-pnorm(abs(r["est"]/r["se"]))),check.names=FALSE)
  }
  do.call(rbind,out)
}))
# Holm across the 4 cohort x outcome primary tests
res$p_holm <- p.adjust(res$p, method="holm")
cat("==== per-cohort CMM3 associations ====\n"); print(res[,c("cohort","outcome","n","OR","lo","hi","p","p_holm")], row.names=FALSE)

meta_fn <- function(b,v){
  k<-length(b); # fixed effect
  wf<-1/v; mu_f<-sum(wf*b)/sum(wf); se_f<-sqrt(1/sum(wf))
  Q<-sum(wf*(b-mu_f)^2); df<-k-1
  tau2<-max(0,(Q-df)/(sum(wf)-sum(wf^2)/sum(wf)))
  wr<-1/(v+tau2); mu_r<-sum(wr*b)/sum(wr)
  q2<-max(1,Q/df); se_hk<-sqrt(q2/sum(wr)); tcrit<-qt(0.975, df=max(1,df))
  I2<-max(0,(Q-df)/Q)*100
  data.frame(k=k, OR_fe=exp(mu_f), fe_lo=exp(mu_f-1.96*se_f), fe_hi=exp(mu_f+1.96*se_f),
             p_fe=2*pnorm(-abs(mu_f/se_f)), OR_re=exp(mu_r), re_lo=exp(mu_r-tcrit*se_hk),
             re_hi=exp(mu_r+tcrit*se_hk), p_re_hk=2*pt(-abs(mu_r/se_hk),df=max(1,df)),
             tau2=tau2, I2=I2, Q=Q)
}
metas <- do.call(rbind, lapply(c("cataract","amd"), function(y){
  s<-res[res$outcome==y,]; m<-meta_fn(s$logOR, s$se^2); cbind(outcome=y,m)
}))
write.csv(res,"C:/ops2/output/cohort_results.csv",row.names=FALSE)
write.csv(metas,"C:/ops2/output/meta_nhanes_nhis.csv",row.names=FALSE)
cat("\n==== 2-cohort meta (CMM3) ====\n"); print(metas,row.names=FALSE)
cat("\nDONE.\n")
