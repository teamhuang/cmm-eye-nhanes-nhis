#!/usr/bin/env Rscript
# Topic 3 — EXTERNAL REPLICATION (ELSA UK + SHARE Europe) + 3-cohort meta-analysis.
# RUNS ONLY on user-supplied harmonised CSVs (see ELSA_SHARE_data_request.md).
# Frozen a priori: primary exposure = full CMM 0-4; outcomes cataract & AMD; Holm across outcomes;
# design-weighted logistic with ns(age,4)+male+educ+smoke; REML + Hartung-Knapp RE meta (small k=3).
# If a cohort file is absent, the script reports it as NOT RUN (never fabricates).
suppressMessages({library(splines)})
RAW <- "H:/眼科论文/眼科2/data/raw"
OUT <- "C:/ops2/output"

fit_cohort <- function(dat, name){
  dat <- dat[dat$age>=40 & !is.na(dat$CMM) & !is.na(dat$cataract) & !is.na(dat$amd), ]
  if(nrow(dat)<200){ cat(sprintf("[%s] too few rows (n=%d) -> skipped\n", name, nrow(dat))); return(NULL) }
  dat$educf <- factor(dat$educ)
  ab <- as.data.frame(ns(dat$age, df=4)); colnames(ab) <- paste0("ag", seq_len(ncol(ab)))
  dat <- cbind(dat, ab); dat$w1 <- dat$wt/mean(dat$wt)
  if(!("psu" %in% names(dat)) || all(dat$psu==1 & dat$strata==1)){
    dat$grp <- seq_len(nrow(dat))  # unweighted / no design -> row jackknife (conservative-ish)
  } else dat$grp <- paste(dat$strata, dat$psu)
  cf <- function(f,dd,fam,nm) unname(coef(glm(f,data=dd,weights=dd$w1,family=fam,control=list(maxit=200)))[nm])
  res <- list(cohort=name, n=nrow(dat))
  for (y in c("cataract","amd")){
    f <- as.formula(paste0(y," ~ CMM + ag1+ag2+ag3+ag4 + male + educf + smoke"))
    pt <- cf(f,dat,binomial(),"CMM")
    g <- unique(dat$grp); j <- sapply(g, function(gr){ dd<-dat[dat$grp!=gr,]; dd$w1<-dd$wt/mean(dd$wt); tryCatch(cf(f,dd,binomial(),"CMM"),error=function(e)NA)})
    j <- j[is.finite(j)]; se <- sqrt((length(j)-1)/length(j)*sum((j-mean(j))^2))
    res[[y]] <- c(logOR=pt, se=se)
    cat(sprintf("[%s] %-9s OR=%.3f (%.3f-%.3f) p=%.3g  n=%d\n", name, y, exp(pt), exp(pt-1.96*se), exp(pt+1.96*se), 2*(1-pnorm(abs(pt/se))), nrow(dat)))
  }
  # optional CRP mediation (ELSA)
  if(all(c("crp") %in% names(dat)) && sum(!is.na(dat$crp))>200){
    dd <- dat[!is.na(dat$crp),]; dd$w1<-dd$wt/mean(dd$wt); dd$log2crp<-log2(dd$crp+0.01)
    fa <- cf(as.formula("log2crp ~ CMM + ag1+ag2+ag3+ag4 + male + educf + smoke"),dd,gaussian(),"CMM")
    cat(sprintf("[%s] CMM->log2CRP a-path=%.3f\n", name, fa))
  }
  res
}

read_if <- function(p){ if(file.exists(p)) tryCatch(read.csv(p), error=function(e) NULL) else NULL }
nhanes <- NULL
# NHANES point estimates come from the completed main analysis (logOR,se) if present
if(file.exists(file.path(OUT,"jackknife_draws.csv"))){
  nj <- read.csv(file.path(OUT,"jackknife_draws.csv")); 
  # recompute NHANES logOR/se for cataract & AMD from saved main table
  mt <- read.csv(file.path(OUT,"main_total_effects.csv"))
  nhanes <- list(cohort="NHANES",
                 cataract=c(logOR=log(mt$OR[mt$outcome=="cataract_surg"]), se=mt$se[mt$outcome=="cataract_surg"]),
                 amd=c(logOR=log(mt$OR[mt$outcome=="amd_self"]), se=mt$se[mt$outcome=="amd_self"]))
}
elsa  <- read_if(file.path(RAW,"elsa","elsa_eye.csv"));  elsa  <- if(!is.null(elsa)) fit_cohort(elsa,"ELSA")
share <- read_if(file.path(RAW,"share","share_eye.csv")); share <- if(!is.null(share)) fit_cohort(share,"SHARE")

coh <- Filter(Negate(is.null), list(nhanes, elsa, share))
if(length(coh)<2){
  cat("\nNOT RUN: fewer than 2 cohorts available. Provide elsa_eye.csv and/or share_eye.csv under data/raw/. This script never imputes cohort results.\n")
  quit(save="no")
}
re_meta <- function(y){
  es <- do.call(rbind, lapply(coh, function(c) c(y)))
  b <- es[,1]; v <- es[,2]^2; k <- length(b)
  # REML tau^2
  tau2 <- 0; for(it in 1:200){ w<-1/(v+tau2); mu<-sum(w*b)/sum(w); 
    new<-max(0,(sum(w*(b-mu)^2)-(k-1-sum(w)+sum(w^2)/sum(w)))/(sum(w)-sum(w^2)/sum(w)))
    if(abs(new-tau2)<1e-8){tau2<-new;break}; tau2<-new }
  w<-1/(v+tau2); mu<-sum(w*b)/sum(w); se<-sqrt(1/sum(w))
  # Hartung-Knapp
  Q<-sum(w*(b-mu)^2); q2<-max(1, Q/(k-1)); se_hk<-sqrt(q2/sum(w)); 
  I2<-max(0,(Q-(k-1))/Q)*100
  data.frame(outcome=y, k=k, OR=exp(mu), lo=exp(mu-1.96*se_hk), hi=exp(mu+1.96*se_hk),
             p_hk=2*pt(-abs(mu/se_hk), df=max(1,k-1)), tau2=tau2, I2=I2)
}
meta <- rbind(re_meta("cataract"), re_meta("amd"))
write.csv(meta, file.path(OUT,"meta_results.csv"), row.names=FALSE)
cat("\n==== 3-cohort random-effects meta (REML + Hartung-Knapp) ====\n"); print(meta, row.names=FALSE)
cat("\nMeta complete. Copy meta_results.csv back to project output.\n")
