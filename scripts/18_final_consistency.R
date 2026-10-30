#!/usr/bin/env Rscript
# FINAL consistency: (1) weighted Table 1 for both cohorts; (2) like-for-like IPC with t-based CI
# (matches the primary svyglm interval); (3) package JKn replicate-design SE to validate jackknife;
# (4) exact NHANES sample-flow chain.
suppressMessages({library(survey); library(splines)})
OUT <- "C:/ops2/output"
if (!file.exists(file.path(OUT,"nhanes_analytic_age40_complete.csv"))) OUT <- "H:/眼科论文/眼科2/analysis/output"

nh <- read.csv(file.path(OUT,"nhanes_analytic_age40_complete.csv")); nh <- nh[nh$educ %in% 1:5,]
nh$cataract<-nh$cataract_surg; nh$amd<-nh$amd_self; nh$smoke<-nh$ever_smoke
nh$CMM3 <- nh$diabetes + nh$hypertension + nh$obesity; nh$CMM4 <- nh$CMM
nh$educ3 <- factor(ifelse(nh$educ<=2,0, ifelse(nh$educ==3,1,2)))
nh <- nh[!is.na(nh$CMM3) & !is.na(nh$CMM4),]
nh$psu<-nh$SDMVPSU; nh$strata<-nh$SDMVSTRA
des_nh <- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~wt, data=nh, nest=TRUE)

ni <- read.csv(file.path(OUT,"nhis_analytic.csv")); ni$educ3<-factor(ni$educ); ni$CMM3<-ni$CMM
des_ni <- svydesign(id=~psu, strata=~strata, weights=~wt, data=ni, nest=TRUE)

# ---------- Table 1 (survey-weighted) ----------
pct <- function(des, form){ m <- svymean(form, des, na.rm=TRUE); sprintf("%.1f", 100*unname(coef(m)[1])) }
ageci <- function(des){ m <- svymean(~age, des, na.rm=TRUE); ci <- confint(m); sprintf("%.1f (%.1f-%.1f)", unname(coef(m)[1]), ci[1], ci[2]) }
edu3 <- function(des){ m<-svymean(~educ3, des, na.rm=TRUE); p<-100*coef(m); sprintf("%.1f / %.1f / %.1f", p[1], p[2], p[3]) }
t1 <- data.frame(
  var=c("Age, y (weighted mean, 95% CI)","Male, % (weighted)","Diabetes, %","Hypertension, %","Obesity, %","Dyslipidaemia, %",
        "Cataract surgery, %","AMD (self-report), %","Education (<HS / HS / >HS), %"),
  NHANES=c(ageci(des_nh), pct(des_nh,~I(male==1)), pct(des_nh,~diabetes), pct(des_nh,~hypertension),
           pct(des_nh,~obesity), pct(des_nh,~dyslipidemia), pct(des_nh,~cataract), pct(des_nh,~amd), edu3(des_nh)),
  NHIS=c(ageci(des_ni), pct(des_ni,~I(male==1)), pct(des_ni,~diabetes), pct(des_ni,~hypertension),
         pct(des_ni,~obesity), "not collected", pct(des_ni,~cataract), pct(des_ni,~amd), edu3(des_ni)),
  stringsAsFactors=FALSE)
print(t1, row.names=FALSE); write.csv(t1, file.path(OUT,"sr2_table1_weighted.csv"), row.names=FALSE)

# ---------- IPC like-for-like with t-based CI ----------
mkSpline <- function(d){ b<-as.data.frame(ns(d$age,df=4)); names(b)<-paste0("ag",1:4); cbind(d,b) }
nhs <- mkSpline(nh)
full <- read.csv(file.path(OUT,"nhanes_full_age40_with_missing.csv")); full<-full[!is.na(full$age),]; full<-mkSpline(full)
full$included <- full$SEQN %in% nh$SEQN
desfull <- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~WTMEC2YR, data=full, nest=TRUE)
ipm <- svyglm(included ~ ag1+ag2+ag3+ag4 + RIAGENDR + RIDRETH1 + SDDSRVYR, design=desfull, family=quasibinomial())
full$p_inc <- predict(ipm, type="response")
nhm <- merge(nhs, full[,c("SEQN","p_inc")], by="SEQN"); nhm<-nhm[is.finite(nhm$p_inc)&nhm$p_inc>0,]; nhm$ipc_w<-nhm$wt/nhm$p_inc
des_cc <- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~wt, data=nhm, nest=TRUE)
des_ipc<- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~ipc_w, data=nhm, nest=TRUE)
rowT <- function(des,y,rhs,tag){ m<-svyglm(as.formula(paste(y,"~",rhs)),des,family=quasibinomial())
  b<-unname(coef(m)["CMM3"]); s<-unname(sqrt(vcov(m)["CMM3","CMM3"])); tc<-qt(0.975,df.residual(m))
  data.frame(analysis=tag, OR=exp(b), lo=exp(b-tc*s), hi=exp(b+tc*s), p=unname(summary(m)$coefficients["CMM3","Pr(>|t|)"]), tcrit=tc) }
P <- "CMM3 + ag1+ag2+ag3+ag4 + male + educ3 + smoke"; Pinc <- paste(P,"+ pir")
ipc <- rbind(rowT(des_cc,"cataract",P,"cat_complete_case_PRIMARY"), rowT(des_ipc,"cataract",P,"cat_IPC_PRIMARY"),
             rowT(des_cc,"cataract",Pinc,"cat_complete_case_INCOME"), rowT(des_ipc,"cataract",Pinc,"cat_IPC_INCOME"),
             rowT(des_cc,"amd",P,"amd_complete_case_PRIMARY"), rowT(des_ipc,"amd",P,"amd_IPC_PRIMARY"))
print(ipc, row.names=FALSE); write.csv(ipc, file.path(OUT,"sr2_ipc_likeforlike_t.csv"), row.names=FALSE)

# ---------- JKn replicate design vs Taylor vs manual jackknife ----------
des_nhs <- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~wt, data=nhs, nest=TRUE)
m <- svyglm(cataract ~ CMM3 + ag1+ag2+ag3+ag4 + male + educ3 + smoke, des_nhs, family=quasibinomial())
taylor <- unname(sqrt(diag(vcov(m)))["CMM3"])
des_jkn <- as.svrepdesign(des_nhs, type="JKn")
m_jkn <- svyglm(cataract ~ CMM3 + ag1+ag2+ag3+ag4 + male + educ3 + smoke, des_jkn, family=quasibinomial())
jkn <- unname(sqrt(diag(vcov(m_jkn)))["CMM3"])
# manual delete-one-PSU jackknife
d<-nhs; d$grp<-paste(d$strata,d$psu); d$w1<-d$wt/mean(d$wt); f<-cataract ~ CMM3 + ag1+ag2+ag3+ag4 + male + educ3 + smoke
g<-unique(d$grp)
j<-sapply(g,function(gr){dd<-d[d$grp!=gr,]; dd$w1<-dd$wt/mean(dd$wt); tryCatch(unname(coef(glm(f,data=dd,weights=w1,family=binomial()))["CMM3"]),error=function(e)NA)})
j<-j[is.finite(j)]; R<-length(j); manual<-sqrt((R-1)/R*sum((j-mean(j))^2))
jk <- data.frame(taylor=taylor, svy_JKn=jkn, manual_delete_a_PSU=manual,
  ratio_JKn_over_taylor=jkn/taylor, ratio_manual_over_taylor=manual/taylor, manual_reps=R)
print(jk, row.names=FALSE); write.csv(jk, file.path(OUT,"sr2_jackknife_validation.csv"), row.names=FALSE)

# ---------- sample flow ----------
full0 <- read.csv(file.path(OUT,"nhanes_full_age40_with_missing.csv"))
flow <- data.frame(step=c("age >= 40 examined (full file)","valid cataract surgery","+ valid AMD",
  "+ complete CMM4 (4 conditions)","+ complete CRP","+ complete covariates (educ 1-5, PIR, smoking)"),
  n=c(nrow(full0), sum(!is.na(full0$cataract_surg)), sum(!is.na(full0$cataract_surg)&!is.na(full0$amd_self)),
      sum(!is.na(full0$cataract_surg)&!is.na(full0$amd_self)&!is.na(full0$CMM)),
      sum(!is.na(full0$cataract_surg)&!is.na(full0$amd_self)&!is.na(full0$CMM)&!is.na(full0$log2crp)),
      nrow(nh)))
flow$removed <- c(NA, -diff(flow$n))
print(flow, row.names=FALSE); write.csv(flow, file.path(OUT,"sr2_sample_flow.csv"), row.names=FALSE)
