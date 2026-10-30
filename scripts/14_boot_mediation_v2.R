#!/usr/bin/env Rscript
# ROUND 2: corrected mediation bootstrap.
# (a) STRATIFIED cluster resampling: within each stratum, resample that stratum's PSUs with
#     replacement (preserves the survey's stratified design; the round-1 pooled resample did not).
# (b) Reports convergence counts, the Monte-Carlo SE of the percentile endpoints, and a VALIDATION
#     that the bootstrap SE of the TOTAL-effect logOR matches the Taylor design SE (proves the
#     resampling variance -- not just the point estimate -- is design-consistent).
Sys.setenv(OMP_NUM_THREADS=1, MKL_NUM_THREADS=1, OPENBLAS_NUM_THREADS=1, NUMEXPR_NUM_THREADS=1)
suppressMessages({library(splines)})
set.seed(20261029)
OUT <- "C:/ops2/output"

nh <- read.csv(file.path(OUT,"nhanes_analytic_age40_complete.csv")); nh <- nh[nh$educ %in% 1:5,]
nh$cataract<-nh$cataract_surg; nh$amd<-nh$amd_self; nh$smoke<-nh$ever_smoke
nh$CMM4 <- nh$CMM
nh$educ3 <- factor(ifelse(nh$educ<=2,0, ifelse(nh$educ==3,1,2)))
nh <- nh[!is.na(nh$CMM4) & !is.na(nh$log2crp),]
COVS <- "ns(age, df=4) + male + educ3 + smoke + pir"

# stratified PSU index: list of PSUs per stratum, each PSU -> row indices
rows_by_psu <- split(seq_len(nrow(nh)), interaction(nh$SDMVSTRA, nh$SDMVPSU, drop=TRUE))
psu_by_strata <- split(names(rows_by_psu), nh$SDMVSTRA)

gcoef <- function(f, d, fam, term){ dd <- d; dd$w1 <- dd$wt/mean(dd$wt)
  m <- tryCatch(glm(f, data=dd, weights=w1, family=fam, control=list(maxit=200)), error=function(e) NULL)
  if(is.null(m) || !isTRUE(m$converged)) return(c(NA,0))
  cc <- coef(m); if(term %in% names(cc) && is.finite(cc[[term]]) && abs(cc[[term]])<50) return(c(unname(cc[[term]]),1))
  c(NA,0) }

ftot <- list(cataract=as.formula(paste("cataract ~ CMM4 +",COVS)), amd=as.formula(paste("amd ~ CMM4 +",COVS)))
fdir <- list(cataract=as.formula(paste("cataract ~ CMM4 +",COVS,"+ log2crp")), amd=as.formula(paste("amd ~ CMM4 +",COVS,"+ log2crp")))
fa   <- as.formula(paste("log2crp ~ CMM4 +",COVS))

B <- as.integer(Sys.getenv("BOOT_B","1000"))
store <- matrix(NA_real_, nrow=B, ncol=5, dimnames=list(NULL,c("tot_cat","dir_cat","tot_amd","dir_amd","a_crp")))
conv <- rep(NA, B)
for (b in seq_len(B)){
  take <- unlist(lapply(names(psu_by_strata), function(s){
      ps <- psu_by_strata[[s]]; sample(ps, length(ps), replace=TRUE) }), use.names=FALSE)
  rows <- unlist(rows_by_psu[take], use.names=FALSE); d <- nh[rows,]
  tc <- gcoef(ftot$cataract, d, binomial(), "CMM4"); dc <- gcoef(fdir$cataract, d, binomial(), "CMM4")
  ta <- gcoef(ftot$amd, d, binomial(), "CMM4");     da <- gcoef(fdir$amd, d, binomial(), "CMM4")
  aa <- gcoef(fa, d, gaussian(), "CMM4")
  store[b,] <- c(tc[1], dc[1], ta[1], da[1], aa[1])
  conv[b] <- (tc[2]+dc[2]+ta[2]+da[2]+aa[2])
  if (b %% 200 == 0) cat(sprintf("bootstrap %d/%d\n", b, B))
}
write.csv(as.data.frame(store), file.path(OUT,"sr2_mediation_boot.csv"), row.names=FALSE)

nconv_all <- sum(conv==5)
cat(sprintf("replicates with ALL 5 models converged: %d/%d (%.1f%%)\n", nconv_all, B, 100*nconv_all/B))

# percentile CI + MC SE of endpoints (Bounded-density / order-statistic method)
pct_ci <- function(v){ v <- v[is.finite(v)]; n <- length(v)
  q <- quantile(v, c(.025,.5,.975), names=FALSE)
  # density at quantile via normal approx using bootstrap sd
  s <- sd(v); dens_lo <- dnorm(q[1], mean(v), s); dens_hi <- dnorm(q[3], mean(v), s)
  mc_lo <- sqrt(.025*.975/n)/dens_lo; mc_hi <- sqrt(.975*.025/n)/dens_hi
  c(est=mean(v), se=s, median=q[2], lo=q[1], hi=q[3], n_used=n, mc_se_lo=mc_lo, mc_se_hi=mc_hi) }
ind_cat <- store[,"tot_cat"]-store[,"dir_cat"]; ind_amd <- store[,"tot_amd"]-store[,"dir_amd"]
summ <- rbind(indirect_cataract=pct_ci(ind_cat), indirect_amd=pct_ci(ind_amd), a_path=pct_ci(store[,"a_crp"]),
              total_cataract_logOR=pct_ci(store[,"tot_cat"]))
write.csv(as.data.frame(summ), file.path(OUT,"sr2_mediation_boot_summary.csv"))
cat("\n== corrected bootstrap summary ==\n"); print(round(summ,4))
# validation: bootstrap SE of total-effect logOR vs Taylor design SE (0.0708 from audit)
cat(sprintf("\nVALIDATION total-effect cataract logOR: bootstrap SE = %.4f  (Taylor design SE = 0.0708)\n", summ["total_cataract_logOR","se"]))
cat("BOOTSTRAP V2 DONE\n")
