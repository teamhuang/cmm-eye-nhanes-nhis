#!/usr/bin/env Rscript
# Topic 3 main analysis — NHANES derivation cohort (v4, design-correct).
# CMM -> cataract surgery (Y1) & self-reported AMD (Y2); log2-CRP mediator.
# Point estimates: design-weighted GLM (survey weights rescaled to mean 1; MLE scale-invariant).
# Variance/CI/p: delete-a-PSU (cluster) JACKKNIFE at the stratum x PSU level — the design-correct
#   variance for NHANES' 2-PSU-per-stratum structure (validated: total-effect SE ~= Taylor 0.0576).
# Holm across the 2 outcomes. No causal verbs. Seed not required (jackknife is deterministic).
suppressMessages({library(splines)})
OUT <- "C:/ops2/output"

d <- read.csv(file.path(OUT, "nhanes_analytic_age40_complete.csv"), stringsAsFactors=FALSE)
d <- d[d$educ %in% 1:5, ]                       # drop invalid education code (9 = refused/DK; 1 obs)
d$educf <- factor(d$educ)
ageB <- as.data.frame(ns(d$age, df=4)); colnames(ageB) <- paste0("ag", seq_len(ncol(ageB)))
d <- cbind(d, ageB); AGP <- paste(paste0("ag", seq_len(ncol(ageB))), collapse=" + ")
d$grp <- paste(d$SDMVSTRA, d$SDMVPSU)
d$w1  <- d$wt/mean(d$wt)

covars <- paste(AGP, "+ male + educf + pir + ever_smoke")
F1 <- as.formula(paste("cataract_surg ~ CMM +", covars))
F2 <- as.formula(paste("amd_self ~ CMM +", covars))
FA <- as.formula(paste("log2crp ~ CMM +", covars))
D1 <- as.formula(paste("cataract_surg ~ CMM + log2crp +", covars))
D2 <- as.formula(paste("amd_self ~ CMM + log2crp +", covars))

cf <- function(f, dat, fam) unname(coef(glm(f, data=dat, weights=w1, family=fam, control=list(maxit=200)))[c("CMM")] )
cfb<- function(f, dat, fam, nm) unname(coef(glm(f, data=dat, weights=w1, family=fam, control=list(maxit=200)))[nm])

statvec <- function(dat){
  t1<-cf(F1,dat,binomial()); t2<-cf(F2,dat,binomial()); a<-cf(FA,dat,gaussian())
  d1<-cfb(D1,dat,binomial(),"CMM"); d2<-cfb(D2,dat,binomial(),"CMM")
  b1<-cfb(D1,dat,binomial(),"log2crp"); b2<-cfb(D2,dat,binomial(),"log2crp")
  c(tot1=t1, tot2=t2, dir1=d1, dir2=d2, a=a, b1=b1, b2=b2,
    ind1=t1-d1, ind2=t2-d2, ab1=a*b1, ab2=a*b2)
}

pt <- statvec(d)

groups <- unique(d$grp); Rl <- length(groups)
jkv <- sapply(groups, function(gr){ dat <- d[d$grp!=gr, ]; dat$w1 <- dat$wt/mean(dat$wt); statvec(dat) })
jkv <- jkv[, complete.cases(t(jkv)), drop=FALSE]   # drop any replicate that failed to fit
Rl  <- ncol(jkv)
jk_mean <- rowMeans(jkv)
sejk <- sqrt((Rl-1)/Rl * rowSums((jkv - jk_mean)^2))   # delete-a-cluster jackknife SE (per statistic)
names(sejk) <- rownames(jkv)
p <- 2*(1 - pnorm(abs(pt/sejk)))
cat("N =", nrow(d), " jackknife replicates kept =", Rl, "/", length(groups), "\n")

tot_tab <- data.frame(
  outcome=c("cataract_surg","amd_self"),
  logOR=c(pt["tot1"], pt["tot2"]), OR=exp(c(pt["tot1"], pt["tot2"])),
  se=c(sejk["tot1"], sejk["tot2"]),
  OR_lo=exp(c(pt["tot1"],pt["tot2"]) - 1.96*c(sejk["tot1"],sejk["tot2"])),
  OR_hi=exp(c(pt["tot1"],pt["tot2"]) + 1.96*c(sejk["tot1"],sejk["tot2"])),
  p=c(p["tot1"], p["tot2"]), check.names=FALSE)
tot_tab$p_holm <- p.adjust(tot_tab$p, method="holm")

med <- data.frame(
  outcome=c("cataract_surg","amd_self"),
  a_path=pt["a"], b_path=c(pt["b1"],pt["b2"]),
  total_logOR=c(pt["tot1"],pt["tot2"]), direct_logOR=c(pt["dir1"],pt["dir2"]),
  indirect=c(pt["ind1"],pt["ind2"]),
  ind_se=c(sejk["ind1"],sejk["ind2"]),
  ind_lo=c(pt["ind1"],pt["ind2"]) - 1.96*c(sejk["ind1"],sejk["ind2"]),
  ind_hi=c(pt["ind1"],pt["ind2"]) + 1.96*c(sejk["ind1"],sejk["ind2"]),
  prop_mediated=c(pt["ind1"]/pt["tot1"], pt["ind2"]/pt["tot2"]),
  p_indirect=c(p["ind1"], p["ind2"]), check.names=FALSE)

write.csv(tot_tab, file.path(OUT,"main_total_effects.csv"), row.names=FALSE)
write.csv(med,     file.path(OUT,"mediation.csv"), row.names=FALSE)
write.csv(as.data.frame(jkv), file.path(OUT,"jackknife_draws.csv"), row.names=FALSE)
saveRDS(list(pt=pt, sejk=sejk, tot_tab=tot_tab, med=med, jkv=jkv), file.path(OUT,"model_objs.rds"))

cat("\n==== TOTAL EFFECTS (OR per CMM unit; jackknife 95%CI) ====\n"); print(tot_tab, row.names=FALSE)
cat("\n==== A-PATH (log2CRP per CMM)=", round(pt["a"],4), " SE=", round(sejk["a"],4), "\n")
cat("\n==== MEDIATION ====\n"); print(med, row.names=FALSE)
cat("\nDONE.\n")
