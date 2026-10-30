#!/usr/bin/env Rscript
# Mediation indirect-effect CI via PSU cluster bootstrap.
# The indirect effect (log-odds, Cline difference) is a function only of POINT estimates, which
# for design-weighted GLM equal the weighted-MLE coefficients (the survey design affects the SE,
# not the coefficient). So the bootstrap resamples PSUs and refits the weighted glm point model --
# identical estimator to the svyglm primary, but robust to duplicated-PSU designs.
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

cl <- paste(nh$SDMVSTRA, nh$SDMVPSU)
uidx <- split(seq_len(nrow(nh)), cl)
ukeys <- names(uidx)
B <- as.integer(Sys.getenv("BOOT_B", "500"))

gfit <- function(f, d, fam){ dd <- d; dd$w1 <- dd$wt/mean(dd$wt)
  coef(glm(f, data=dd, weights=w1, family=fam, control=list(maxit=200))) }

res <- vector("list", B)
for (b in seq_len(B)){
  take <- sample(ukeys, length(ukeys), replace=TRUE)
  rows <- unlist(uidx[take], use.names=FALSE); d <- nh[rows,]
  out <- c(iter=b)
  for (y in c("cataract","amd")){
    ft <- as.formula(paste(y,"~ CMM4 +",COVS)); fd <- as.formula(paste(y,"~ CMM4 +",COVS,"+ log2crp"))
    mt <- tryCatch(gfit(ft, d, binomial())["CMM4"], error=function(e)NA)
    md <- tryCatch(gfit(fd, d, binomial())["CMM4"], error=function(e)NA)
    if(is.na(mt)||!is.finite(mt)) mt <- NA; if(is.na(md)||!is.finite(md)) md <- NA
    out[[paste0("tot_",y)]] <- unname(mt); out[[paste0("dir_",y)]] <- unname(md)
    out[[paste0("ind_",y)]] <- if(is.finite(mt)&&is.finite(md)) mt-md else NA_real_
  }
  fa <- as.formula(paste("log2crp ~ CMM4 +",COVS))
  ma <- tryCatch(gfit(fa, d, gaussian())["CMM4"], error=function(e)NA); out[["a_crp"]] <- unname(ma)
  res[[b]] <- as.data.frame(as.list(out))
  if (b %% 100 == 0){ cat(sprintf("bootstrap %d/%d done\n", b, B)) }
}
res <- do.call(rbind, res); write.csv(res, file.path(OUT,"sr_mediation_boot.csv"), row.names=FALSE)

summ <- function(v){ v<-v[is.finite(v)]; c(est=mean(v), se=sd(v), lo=quantile(v,0.025), hi=quantile(v,0.975), n=length(v)) }
tab <- data.frame(quantity=c("indirect_cataract","indirect_amd","a_path_log2crp"),
  rbind(summ(res$ind_cataract), summ(res$ind_amd), summ(res$a_crp)))
write.csv(tab, file.path(OUT,"sr_mediation_boot_summary.csv"), row.names=FALSE)
cat("BOOTSTRAP MEDIATION DONE\n"); print(tab)
