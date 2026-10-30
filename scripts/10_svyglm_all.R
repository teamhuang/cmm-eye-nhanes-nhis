#!/usr/bin/env Rscript
# Reviewer-response round 1 — ALL analyses re-derived on one estimator: design-weighted
# svyglm (Taylor linearization) via svydesign(PSU,strata,weight). NHANES 2005-08 + NHIS 2023.
# Outputs sr_*.csv in OUT. Deterministic (no seeds needed except the separate bootstrap script 11).
suppressMessages({library(survey); library(splines)})
OUT <- "C:/ops2/output"

# ---------- data ----------
nh <- read.csv(file.path(OUT,"nhanes_analytic_age40_complete.csv")); nh <- nh[nh$educ %in% 1:5,]
nh$cataract<-nh$cataract_surg; nh$amd<-nh$amd_self; nh$smoke<-nh$ever_smoke
nh$CMM3 <- nh$diabetes + nh$hypertension + nh$obesity
nh$CMM4 <- nh$CMM
nh$educ3 <- factor(ifelse(nh$educ<=2,0, ifelse(nh$educ==3,1,2)))
nh <- nh[!is.na(nh$CMM3) & !is.na(nh$CMM4),]
nh$agegrp <- factor(cut(nh$age, c(39,59,69,200), labels=c("40-59","60-69","70+")))
des_nh <- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~wt, data=nh, nest=TRUE)

ni <- read.csv(file.path(OUT,"nhis_analytic.csv"))
ni$educ3 <- factor(ni$educ); ni$CMM3 <- ni$CMM
ni$agegrp <- factor(cut(ni$age, c(39,59,69,200), labels=c("40-59","60-69","70+")))
des_ni <- svydesign(id=~psu, strata=~strata, weights=~wt, data=ni, nest=TRUE)

AG <- "ns(age, df=4)"
COV <- paste(AG, "+ male + educ3 + smoke")          # standard set (NHANES adds pir in M3 only)

# ---------- helpers ----------
or_ci <- function(m, term){
  cf <- coef(m); s <- sqrt(diag(vcov(m))); if(!term %in% names(cf)) return(c(OR=NA,lo=NA,hi=NA,p=NA))
  est <- unname(cf[term]); se <- unname(s[term]); df <- df.residual(m); tc <- qt(0.975, df)
  p <- unname(summary(m)$coefficients[term,"Pr(>|t|)"])
  r <- c(exp(est), exp(est-tc*se), exp(est+tc*se), p); names(r) <- c("OR","lo","hi","p"); r
}
orvec <- function(m, term){ r <- or_ci(m, term); data.frame(OR=r[["OR"]], lo=r[["lo"]], hi=r[["hi"]], p=r[["p"]]) }
rowof <- function(cohort,outcome,extra,term,m){ cbind(data.frame(cohort=cohort,outcome=outcome,extra), orvec(m,term)) }

# ---------- A. design + attrition + selection ----------
write.csv(data.frame(cohort=c("NHANES","NHIS"), n=c(nrow(nh),nrow(ni)),
  strata=c(length(unique(nh$SDMVSTRA)),length(unique(ni$strata))),
  psu=c(length(unique(paste(nh$SDMVSTRA,nh$SDMVPSU))),length(unique(paste(ni$strata,ni$psu)))),
  events_cataract=c(sum(nh$cataract),sum(ni$cataract)), events_amd=c(sum(nh$amd),sum(ni$amd))),
  file.path(OUT,"sr_design.csv"), row.names=FALSE)
full <- read.csv(file.path(OUT,"nhanes_full_age40_with_missing.csv"))
write.csv(data.frame(step=c("age>=40 examined/interviewed subset","complete exposure+mediator+covariates+outcomes (analysis)"),
  n=c(nrow(full), nrow(nh))), file.path(OUT,"sr_attrition_nhanes.csv"), row.names=FALSE)
full$inc <- full$SEQN %in% nh$SEQN
write.csv(data.frame(var=c("age","male","hypertension","obesity","diabetes","cataract_surg","amd_self"),
  included=c(mean(nh$age),mean(nh$male),mean(nh$hypertension,na.rm=TRUE),mean(nh$obesity,na.rm=TRUE),mean(nh$diabetes,na.rm=TRUE),mean(nh$cataract,na.rm=TRUE),mean(nh$amd,na.rm=TRUE)),
  excluded=c(mean(full$age[!full$inc],na.rm=TRUE), mean((full$RIAGENDR[!full$inc]==1),na.rm=TRUE), mean(full$hypertension[!full$inc],na.rm=TRUE), mean(full$obesity[!full$inc],na.rm=TRUE), mean(full$diabetes[!full$inc],na.rm=TRUE), mean(full$cataract_surg[!full$inc],na.rm=TRUE), mean(full$amd_self[!full$inc],na.rm=TRUE))),
  file.path(OUT,"sr_selection_nhanes.csv"), row.names=FALSE)

# ---------- B. NHANES primary CMM4 (total, a-path, direct) ----------
prim <- list()
mtot <- svyglm(as.formula(paste("cataract ~ CMM4 +", COV, "+ pir")), des_nh, family=quasibinomial())
mdat <- svyglm(as.formula(paste("amd ~ CMM4 +", COV, "+ pir")), des_nh, family=quasibinomial())
ma   <- svyglm(as.formula(paste("log2crp ~ CMM4 +", COV, "+ pir")), des_nh, family=gaussian())
mdirc<- svyglm(as.formula(paste("cataract ~ CMM4 +", COV, "+ pir + log2crp")), des_nh, family=quasibinomial())
mdiramd<- svyglm(as.formula(paste("amd ~ CMM4 +", COV, "+ pir + log2crp")), des_nh, family=quasibinomial())
crp_coef <- coef(ma)["CMM4"]; crp_se <- sqrt(diag(vcov(ma)))["CMM4"]
prim <- rbind(
  cbind(data.frame(model="CMM4_total",term="cataract"), orvec(mtot,"CMM4")),
  cbind(data.frame(model="CMM4_total",term="amd"),      orvec(mdat,"CMM4")),
  cbind(data.frame(model="CMM4_direct_adjCRP",term="cataract"), orvec(mdirc,"CMM4")),
  cbind(data.frame(model="CMM4_direct_adjCRP",term="amd"),      orvec(mdiramd,"CMM4")))
write.csv(prim, file.path(OUT,"sr_nhanes_primary_svyglm.csv"), row.names=FALSE)
write.csv(data.frame(path="a_CMM4_to_log2crp", est=as.numeric(crp_coef), se=as.numeric(crp_se),
  lo=as.numeric(crp_coef-1.96*crp_se), hi=as.numeric(crp_coef+1.96*crp_se),
  p=as.numeric(summary(ma)$coefficients["CMM4","Pr(>|t|)"])), file.path(OUT,"sr_apath_svyglm.csv"), row.names=FALSE)

# ---------- C. per-cohort CMM3 + meta ----------
coh <- list()
for (cohort in c("NHANES","NHIS")){
  des <- if(cohort=="NHANES") des_nh else des_ni
  for (y in c("cataract","amd")){
    m <- svyglm(as.formula(paste(y,"~ CMM3 +", COV)), des, family=quasibinomial())
    coh[[length(coh)+1]] <- rowof(cohort,y,"CMM3","CMM3",m) } }
coh <- do.call(rbind,coh); names(coh)[names(coh)=="extra"]<-"term"
write.csv(coh, file.path(OUT,"sr_cohort_cmm3_svyglm.csv"), row.names=FALSE)

meta_fn <- function(b,v){
  k<-length(b); wf<-1/v; mu_f<-sum(wf*b)/sum(wf); se_f<-sqrt(1/sum(wf))
  Q<-sum(wf*(b-mu_f)^2); dfq<-k-1
  tau2<-max(0,(Q-dfq)/(sum(wf)-sum(wf^2)/sum(wf)))
  wr<-1/(v+tau2); mu_r<-sum(wr*b)/sum(wr)
  q2<-max(1,Q/dfq); se_hk<-sqrt(q2/sum(wr)); tc<-qt(0.975,max(1,dfq))
  I2<-max(0,(Q-dfq)/Q)*100
  data.frame(k=k, OR_fe=exp(mu_f), fe_lo=exp(mu_f-1.96*se_f), fe_hi=exp(mu_f+1.96*se_f),
    p_fe=2*pnorm(-abs(mu_f/se_f)), OR_re=exp(mu_r), re_lo=exp(mu_r-tc*se_hk), re_hi=exp(mu_r+tc*se_hk),
    p_re_hk=2*pt(-abs(mu_r/se_hk),df=max(1,dfq)), tau2=tau2, I2=I2, Q=Q)
}
metar<-list()
for (y in c("cataract","amd")){
  sel <- coh[coh$outcome==y,]
  b <- log(sel$OR); v <- ((log(sel$hi)-log(sel$lo))/(2*1.96))^2
  r <- meta_fn(b,v); r$outcome<-y; r$cohorts<-paste(sel$cohort,round(sel$OR,3),collapse=" | ")
  metar[[length(metar)+1]]<-r }
meta <- do.call(rbind,metar); rownames(meta)<-NULL
write.csv(meta[,c("outcome","cohorts","k","OR_fe","fe_lo","fe_hi","p_fe","OR_re","re_lo","re_hi","p_re_hk","tau2","I2","Q")],
  file.path(OUT,"sr_meta_cmm3_svyglm.csv"), row.names=FALSE)

# ---------- D. CMM3 vs CMM4 decomposition ----------
dec <- list()
for (y in c("cataract","amd")){
  m <- svyglm(as.formula(paste(y,"~ CMM3 + dyslipidemia +", COV, "+ pir")), des_nh, family=quasibinomial())
  c3<-or_ci(m,"CMM3"); dy<-or_ci(m,"dyslipidemia")
  dec[[length(dec)+1]]<-data.frame(outcome=y, CMM3_OR=as.numeric(c3["OR"]),CMM3_lo=as.numeric(c3["lo"]),CMM3_hi=as.numeric(c3["hi"]),CMM3_p=as.numeric(c3["p"]),
    dys_OR=as.numeric(dy["OR"]),dys_lo=as.numeric(dy["lo"]),dys_hi=as.numeric(dy["hi"]),dys_p=as.numeric(dy["p"]), dys_prev=mean(nh$dyslipidemia,na.rm=TRUE)) }
write.csv(do.call(rbind,dec), file.path(OUT,"sr_cmm3_vs_cmm4.csv"), row.names=FALSE)

# ---------- E. age strata + CMM x age interaction ----------
str_rows<-list(); int_rows<-list()
for (cohort in c("NHANES","NHIS")){
  des <- if(cohort=="NHANES") des_nh else des_ni; dat <- if(cohort=="NHANES") nh else ni
  for (y in c("cataract","amd")){
    for (lv in levels(dat$agegrp)){
      sub <- subset(des, agegrp==lv)
      if(sum(dat[[y]][dat$agegrp==lv])<5){ next }
      m <- tryCatch(svyglm(as.formula(paste(y,"~ CMM3 +", COV)), sub, family=quasibinomial()), error=function(e) NULL)
      if(is.null(m)) next
      str_rows[[length(str_rows)+1]] <- cbind(data.frame(cohort=cohort,outcome=y,agegrp=lv,n=nrow(dat[dat$agegrp==lv,]),
        events=sum(dat[[y]][dat$agegrp==lv])), orvec(m,"CMM3"))
    }
    mi <- svyglm(as.formula(paste(y,"~ CMM3*agegrp +", COV)), des, family=quasibinomial())
    inms <- grep("CMM3:agegrp", names(coef(mi)), value=TRUE)
    if(length(inms)>=1){ bb <- coef(mi)[inms]; Vv <- vcov(mi)[inms,inms,drop=FALSE]
      chi2 <- as.numeric(t(bb) %*% solve(Vv) %*% bb); dfw <- length(inms)
      int_rows[[length(int_rows)+1]] <- data.frame(cohort=cohort,outcome=y,df=dfw,chi2=chi2,p_int=pchisq(chi2,dfw,lower.tail=FALSE)) }
  }
}
write.csv(do.call(rbind,str_rows), file.path(OUT,"sr_age_strata.csv"), row.names=FALSE)
write.csv(do.call(rbind,int_rows), file.path(OUT,"sr_age_interaction.csv"), row.names=FALSE)

# ---------- F. progressive models ----------
prog<-list()
for (cohort in c("NHANES","NHIS")){
  des <- if(cohort=="NHANES") des_nh else des_ni; dat <- if(cohort=="NHANES") nh else ni
  for (y in c("cataract","amd")){
    specs <- list(M1=paste("CMM3 +",AG), M2=paste0("CMM3 + ",COV))
    if(cohort=="NHANES") specs$M3 <- paste0("CMM3 + ",COV," + pir")
    for (mn in names(specs)){
      m <- svyglm(as.formula(paste(y,"~",specs[[mn]])), des, family=quasibinomial())
      prog[[length(prog)+1]] <- cbind(data.frame(cohort=cohort,outcome=y,model=mn), orvec(m,"CMM3")) } } }
write.csv(do.call(rbind,prog), file.path(OUT,"sr_progressive.csv"), row.names=FALSE)

# ---------- G. dose-response + linear trend + nonlinearity ----------
dose<-list(); nr<-list()
for (cohort in c("NHANES","NHIS")){
  des <- if(cohort=="NHANES") des_nh else des_ni; dat <- if(cohort=="NHANES") nh else ni
  dat$CMMc <- factor(dat$CMM3, levels=0:3)
  des <- svydesign(id=if(cohort=="NHANES") ~SDMVPSU else ~psu, strata=if(cohort=="NHANES") ~SDMVSTRA else ~strata,
                   weights=~wt, data=dat, nest=TRUE)
  for (y in c("cataract","amd")){
    mcat <- svyglm(as.formula(paste(y,"~ CMMc +", COV)), des, family=quasibinomial())
    for (lv in c("1","2","3")){ nm<-paste0("CMMc",lv); if(nm %in% names(coef(mcat)))
      dose[[length(dose)+1]] <- cbind(data.frame(cohort=cohort,outcome=y,cat=paste0("CMM3=",lv," vs 0")), orvec(mcat,nm)) }
    mlin <- svyglm(as.formula(paste(y,"~ CMM3 +", COV)), des, family=quasibinomial())
    dose[[length(dose)+1]] <- cbind(data.frame(cohort=cohort,outcome=y,cat="linear_trend_per_condition"), orvec(mlin,"CMM3"))
    # nonlinearity: add quadratic term to the linear model
    dat$CMMsq <- dat$CMM3^2; desq <- svydesign(id=if(cohort=="NHANES") ~SDMVPSU else ~psu, strata=if(cohort=="NHANES") ~SDMVSTRA else ~strata, weights=~wt, data=dat, nest=TRUE)
    mq <- svyglm(as.formula(paste(y,"~ CMM3 + CMMsq +", COV)), desq, family=quasibinomial())
    nq <- or_ci(mq,"CMMsq"); nr[[length(nr)+1]] <- data.frame(cohort=cohort,outcome=y,p_nonlinearity=nq["p"])
  }
}
write.csv(do.call(rbind,dose), file.path(OUT,"sr_dose_response.csv"), row.names=FALSE)
write.csv(do.call(rbind,nr), file.path(OUT,"sr_nonlinearity.csv"), row.names=FALSE)

# ---------- H. NHANES sensitivity (objective AMD, no-obesity, access-adjusted) ----------
sens<-list()
mobj <- subset(des_nh, !is.na(amd_obj_any))
tryCatch({ mo <- svyglm(amd_obj_any ~ CMM3 + ns(age,df=4) + male + educ3 + smoke + pir, mobj, family=quasibinomial())
  sens[[length(sens)+1]] <- cbind(data.frame(analysis="objective_AMD_retinal", outcome="amd_obj_any", note=sprintf("cases=%d", sum(nh$amd_obj_any,na.rm=TRUE))), orvec(mo,"CMM3")) }, error=function(e) NULL)
# no-obesity CMM = diabetes+hypertension only
nh$CMM2 <- nh$diabetes + nh$hypertension
des2 <- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~wt, data=nh, nest=TRUE)
for (y in c("cataract","amd")){ m <- svyglm(as.formula(paste(y,"~ CMM2 +", COV, "+ pir")), des2, family=quasibinomial())
  sens[[length(sens)+1]] <- cbind(data.frame(analysis="no_obesity_CMM2", outcome=y, note=""), orvec(m,"CMM2")) }
# access-adjusted (recent eye exam VIQ017/VIQ056 as proxy)
nh$access <- ifelse(is.na(nh$VIQ056), 1, ifelse(nh$VIQ056==1,1,0))
desa <- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~wt, data=nh, nest=TRUE)
for (y in c("cataract","amd")){ m <- svyglm(as.formula(paste(y,"~ CMM3 +", COV, "+ pir + access")), desa, family=quasibinomial())
  sens[[length(sens)+1]] <- cbind(data.frame(analysis="access_adjusted", outcome=y, note=""), orvec(m,"CMM3")) }
write.csv(do.call(rbind,sens), file.path(OUT,"sr_sensitivity_svyglm.csv"), row.names=FALSE)

cat("ALL svyglm ANALYSES DONE\n")
