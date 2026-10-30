#!/usr/bin/env Rscript
# Reviewer-response analyses (round 1) — NHANES + NHIS. Weighted GLM + delete-a-PSU jackknife.
# All association models adjust for age with a natural cubic spline (df=4) unless a coarse age
# stratum is the stratifier. Outputs sr_*.csv in OUT.
suppressMessages({library(splines)})
OUT <- "C:/ops2/output"

# ---------- data ----------
nh <- read.csv(file.path(OUT,"nhanes_analytic_age40_complete.csv")); nh <- nh[nh$educ %in% 1:5,]
nh$cataract<-nh$cataract_surg; nh$amd<-nh$amd_self; nh$smoke<-nh$ever_smoke
nh$psu<-nh$SDMVPSU; nh$strata<-nh$SDMVSTRA
nh$CMM3 <- nh$diabetes + nh$hypertension + nh$obesity
nh$CMM4 <- nh$CMM
nh$educ3 <- factor(ifelse(nh$educ<=2,0, ifelse(nh$educ==3,1,2)))
nh <- nh[!is.na(nh$CMM3) & !is.na(nh$CMM4),]
nh$agegrp <- factor(cut(nh$age, c(39,59,69,200), labels=c("40-59","60-69","70+")))

ni <- read.csv(file.path(OUT,"nhis_analytic.csv"))
ni$educ3 <- factor(ni$educ); ni$CMM3 <- ni$CMM
ni$agegrp <- factor(cut(ni$age, c(39,59,69,200), labels=c("40-59","60-69","70+")))

# age spline on each cohort
addspline <- function(d){ ab<-as.data.frame(ns(d$age,df=4)); colnames(ab)<-paste0("ag",seq_len(ncol(ab))); cbind(d,ab) }
nh <- addspline(nh); ni <- addspline(ni)
AG <- "ag1 + ag2 + ag3 + ag4"
COV <- paste(AG, "+ male + educ3 + smoke")   # standard adjusted covariate set (no pir)

# ---------- helpers ----------
mm <- function(spec, dat){ X <- model.matrix(as.formula(paste("~",spec)), dat); X <- X[, colnames(X)!="(Intercept)", drop=FALSE]; X <- X[, colSums(is.na(X))==0, drop=FALSE]; colnames(X) <- make.names(colnames(X)); X }
fitX <- function(y, X, w){ df<-data.frame(y=y, X, w1=w); f<-as.formula(paste("y ~", paste(colnames(X),collapse=" + ")));
  glm(f, data=df, weights=w1, family=binomial(), control=list(maxit=500)) }
jack <- function(y, X, dat){
  ncoef <- ncol(X) + 1
  safe <- function(mdl){ if(is.null(mdl)||inherits(mdl,"try-error")||!isTRUE(mdl$converged)) return(rep(NA,ncoef))
    cc <- coef(mdl); if(any(!is.finite(cc)) || any(abs(cc) > 50)) return(rep(NA,ncoef)); cc }
  pt <- safe(tryCatch(fitX(y, X, dat$wt/mean(dat$wt)), error=function(e) NULL))
  g <- unique(paste(dat$strata, dat$psu))
  one <- function(gr){ m <- paste(dat$strata, dat$psu)!=gr; dd<-dat[m,]
    safe(tryCatch(fitX(y[m], X[m,,drop=FALSE], dd$wt/mean(dd$wt)), error=function(e) NULL)) }
  reps <- t(sapply(g, one))
  keep <- apply(reps, 1, function(r) all(is.finite(r))); reps <- reps[keep, , drop=FALSE]; R <- nrow(reps)
  if(R < length(g)) cat(sprintf("jack: dropped %d/%d replicates (non-converged/separating)\n", length(g)-R, length(g)))
  mu <- colMeans(reps); V <- ((R-1)/R)*crossprod(reps - mu)
  list(pt=pt, mu=mu, V=V, names=names(pt), R=R, Rtot=length(g))
}
cirow <- function(j,nm){ i<-match(nm,j$names); est<-j$pt[i]; s<-sqrt(j$V[i,i]); c(OR=exp(est), lo=exp(est-1.96*s), hi=exp(est+1.96*s), p=2*(1-pnorm(abs(est/s)))) }
wald <- function(j, nms){ idx<-match(nms, j$names); idx<-idx[!is.na(idx)]; b<-j$pt[idx]; Vs<-j$V[idx,idx,drop=FALSE]
  W<-as.numeric(t(b)%*%solve(Vs)%*%b); df<-length(idx); c(chi2=W, df=df, p=pchisq(W,df,lower.tail=FALSE)) }

# ---------- 1. design + attrition + selection ----------
full <- read.csv(file.path(OUT,"nhanes_full_age40_with_missing.csv"))
write.csv(data.frame(cohort=c("NHANES","NHIS"), n=c(nrow(nh),nrow(ni)),
  strata=c(length(unique(nh$strata)),length(unique(ni$strata))),
  psu=c(length(unique(paste(nh$strata,nh$psu))),length(unique(paste(ni$strata,ni$psu)))),
  events_cataract=c(sum(nh$cataract),sum(ni$cataract)), events_amd=c(sum(nh$amd),sum(ni$amd))),
  file.path(OUT,"sr_design.csv"), row.names=FALSE)
write.csv(data.frame(step=c("age>=40 examined/interviewed subset","complete exposure+mediator+covariates+outcomes (analysis)"),
  n=c(nrow(full), nrow(nh))), file.path(OUT,"sr_attrition_nhanes.csv"), row.names=FALSE)
full$inc <- full$SEQN %in% nh$SEQN
write.csv(data.frame(var=c("age","male","hypertension","obesity","diabetes","cataract_surg","amd_self"),
  included=c(mean(nh$age),mean(nh$male),mean(nh$hypertension,na.rm=TRUE),mean(nh$obesity,na.rm=TRUE),mean(nh$diabetes,na.rm=TRUE),mean(nh$cataract,na.rm=TRUE),mean(nh$amd,na.rm=TRUE)),
  excluded=c(mean(full$age[!full$inc],na.rm=TRUE), mean((full$RIAGENDR[!full$inc]==1),na.rm=TRUE), mean(full$hypertension[!full$inc],na.rm=TRUE), mean(full$obesity[!full$inc],na.rm=TRUE), mean(full$diabetes[!full$inc],na.rm=TRUE), mean(full$cataract_surg[!full$inc],na.rm=TRUE), mean(full$amd_self[!full$inc],na.rm=TRUE))),
  file.path(OUT,"sr_selection_nhanes.csv"), row.names=FALSE)

# ---------- 2. CMM3 vs CMM4 decomposition (NHANES; dyslipidemia as its own term) ----------
dec <- list()
for (y in c("cataract","amd")){
  j <- jack(nh[[y]], mm(paste0("CMM3 + dyslipidemia + ", COV, " + pir"), nh), nh)
  c3<-cirow(j,"CMM3"); dy<-cirow(j,"dyslipidemia")
  dec[[y]] <- data.frame(outcome=y, CMM3_OR=as.numeric(c3["OR"]),CMM3_lo=as.numeric(c3["lo"]),CMM3_hi=as.numeric(c3["hi"]),CMM3_p=as.numeric(c3["p"]),
     dys_OR=as.numeric(dy["OR"]),dys_lo=as.numeric(dy["lo"]),dys_hi=as.numeric(dy["hi"]),dys_p=as.numeric(dy["p"]), dys_prev=mean(nh$dyslipidemia,na.rm=TRUE))
}
write.csv(do.call(rbind,dec), file.path(OUT,"sr_cmm3_vs_cmm4.csv"), row.names=FALSE)

# ---------- 3. age stratification + CMM x age interaction ----------
str_rows<-list(); int_rows<-list()
for (cohort in c("NHANES","NHIS")){
  dat <- if(cohort=="NHANES") nh else ni
  for (y in c("cataract","amd")){
    for (lv in levels(dat$agegrp)){
      sub <- dat[dat$agegrp==lv,]
      if(sum(sub[[y]])<5 || nrow(sub)<50){ next }
      j <- jack(sub[[y]], mm(paste0("CMM3 + ", COV), sub), sub)
      str_rows[[length(str_rows)+1]] <- cbind(data.frame(cohort=cohort,outcome=y,agegrp=lv,n=nrow(sub),events=sum(sub[[y]])), as.data.frame(t(cirow(j,"CMM3"))))
    }
    Xi <- mm(paste0("CMM3*agegrp + ", AG, " + male + educ3 + smoke"), dat)
    inms <- grep("CMM3\\.agegrp", colnames(Xi), value=TRUE)
    ji <- jack(dat[[y]], Xi, dat); w <- wald(ji, inms)
    int_rows[[length(int_rows)+1]] <- data.frame(cohort=cohort,outcome=y,terms=length(inms),chi2=w["chi2"],df=w["df"],p_int=w["p"])
  }
}
write.csv(do.call(rbind,str_rows), file.path(OUT,"sr_age_strata.csv"), row.names=FALSE)
write.csv(do.call(rbind,int_rows), file.path(OUT,"sr_age_interaction.csv"), row.names=FALSE)

# ---------- 4. progressive models ----------
prog<-list()
for (cohort in c("NHANES","NHIS")){
  dat <- if(cohort=="NHANES") nh else ni
  for (y in c("cataract","amd")){
    specs <- list(M1=paste("CMM3 +",AG), M2=paste0("CMM3 + ",COV))
    if(cohort=="NHANES") specs$M3 <- paste0("CMM3 + ",COV," + pir")
    for (mn in names(specs)){
      j <- jack(dat[[y]], mm(specs[[mn]], dat), dat)
      prog[[length(prog)+1]] <- cbind(data.frame(cohort=cohort,outcome=y,model=mn,n=nrow(dat)), as.data.frame(t(cirow(j,"CMM3"))))
    }
  }
}
write.csv(do.call(rbind,prog), file.path(OUT,"sr_progressive.csv"), row.names=FALSE)

# ---------- 5. dose-response (CMM3 categories) + linear trend + nonlinearity ----------
dose<-list(); nr<-list()
for (cohort in c("NHANES","NHIS")){
  dat <- if(cohort=="NHANES") nh else ni; dat$CMMc <- factor(dat$CMM3, levels=0:3)
  for (y in c("cataract","amd")){
    jcat <- jack(dat[[y]], mm(paste0("CMMc + ",COV), dat), dat)
    for (lv in c("1","2","3")){ nm<-paste0("CMMc",lv); if(nm %in% jcat$names) dose[[length(dose)+1]] <- cbind(data.frame(cohort=cohort,outcome=y,cat=paste0("CMM3=",lv)), as.data.frame(t(cirow(jcat,nm)))) }
    jlin <- jack(dat[[y]], mm(paste0("CMM3 + ",COV), dat), dat)
    dose[[length(dose)+1]] <- cbind(data.frame(cohort=cohort,outcome=y,cat="linear_trend"), as.data.frame(t(cirow(jlin,"CMM3"))))
    jn <- jack(dat[[y]], mm(paste0("factor(CMM3) + ",COV), dat), dat)
    cnms <- grep("^factor\\.CMM3", colnames(mm(paste0("factor(CMM3) + ",COV), dat)), value=TRUE)
    if(length(cnms)>=2){ w <- wald(jn, cnms); nr[[length(nr)+1]] <- data.frame(cohort=cohort,outcome=y,chi2=w["chi2"],df=w["df"],p_nonlin=w["p"]) }
  }
}
write.csv(do.call(rbind,dose), file.path(OUT,"sr_dose_response.csv"), row.names=FALSE)
write.csv(do.call(rbind,nr), file.path(OUT,"sr_nonlinearity.csv"), row.names=FALSE)
cat("REVIEWER ANALYSES DONE\n")
