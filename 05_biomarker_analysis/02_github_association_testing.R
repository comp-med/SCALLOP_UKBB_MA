#!/usr/bin/env Rscript

## script to run association testing between proteins and phecodes
## Maik Pietzner 24/05/2023
rm(list=ls())

## get the arguments from the command line
args <- commandArgs(trailingOnly=T)

## little options
options(stringsAsFactors = F)
## avoid conversion of numbers
options(scipen = 1)

## correct directory
setwd("<path to file>")

## packages needed
require(data.table)
require(doMC)
require(rms)
require(survival)

## import phecode to test across proteins
phe    <- args[1]
cohort <- args[2]

cat("run regression with", phe, "in", cohort, "\n")

#----------------------------#
##-- import relevant data --##
#----------------------------#

## covariates
ukb.cov  <- fread("input/UKB.covariates.20230523.txt")
## proteins
ukb.prot <- fread("input/UKB.proteins.20230523.txt")
## relevant phecode
ukb.phe  <- fread("input/UKB.phecodes.update.20230524.txt", select = c("eid", paste0("phecode_", phe, "_bin"), paste0("phecode_", phe, "_time")))
## combine everything
ukb.dat  <- merge(ukb.cov, ukb.phe, by.x="f.eid", by.y="eid")
ukb.dat  <- merge(ukb.dat, ukb.prot, by.x="f.eid", by.y="eid")
## delete what is no longer needed
rm(ukb.cov); rm(ukb.prot); rm(ukb.phe); gc(reset=T)
## drop prevalent cases
ukb.dat  <- ukb.dat[ !is.na(eval(as.name(paste0("phecode_", phe, "_bin")))) ]

## import protein label
lab.prot <- fread("input/label.Olink.proteins.20230523.txt")

#----------------------------#
##--    run associations  --##
#----------------------------#

## compute follow-up time
summary(ukb.dat[, paste0("phecode_", phe, "_time"), with=F])
## drop cases within first six months (given in years)
ukb.dat  <- ukb.dat[ !(eval(as.name(paste0("phecode_", phe, "_bin"))) == 1 & eval(as.name(paste0("phecode_", phe, "_time"))) <= .5) ]

## decide whether or not the outcome is sex-specific
if(cohort != "Both"){
  ukb.dat <- ukb.dat[ sex == cohort ]
  ## formula for adjustment
  adj      <- "age + I(age^2) + log_cysc_cleaned + bmi + smoking + alcohol + pc1 + pc2 + pc3 + rcs(month_blood, c(1,6,12)) + rcs(time_blood.num, c(9,14,20)) + rcs(fast.0, c(.5,3,20)) + rcs(sample_age, c(11,12.7,15))"
}else{
  ## formula for adjustment
  adj      <- "age * sex + I(age^2) + log_cysc_cleaned + bmi + smoking + alcohol + pc1 + pc2 + pc3 + rcs(month_blood, c(1,6,12)) + rcs(time_blood.num, c(9,14,20)) + rcs(fast.0, c(.5,3,20)) + rcs(sample_age, c(11,12.7,15))"
}

## use model matrix to ease coding
tmp           <- model.matrix.lm(as.formula(paste0("~ ", adj)), ukb.dat, na.action = "na.pass")
## edit colnames
colnames(tmp) <- gsub(" |,|\\(|\\)|\\^|'", "_", colnames(tmp))
## define new adjustment set
adj           <- paste(colnames(tmp)[-1], collapse = " + ")
## create new data set
ukb.dat       <- data.table(ukb.dat[, c("f.eid", paste0("phecode_", phe, "_bin"), paste0("phecode_", phe, "_time"), lab.prot$id), with=F], tmp[,-1])
rm(tmp); gc(reset=T)

## intermediate output
cat("run regression with ", nrow(ukb.dat), " samples in total\n")

## do in parallel
registerDoMC(12)

## run testing
res.surv <- mclapply(lab.prot$id, function(x){
  
  print(x)
  
  ## run model: implement cox-prop test and possible step function for time intervals
  ff   <- coxph(as.formula(paste0("Surv(phecode_", phe, "_time,","phecode_", phe,"_bin) ~ ", x, " + ", adj)), ukb.dat, ties = "breslow")
  ## cox-prop test
  ff.p <- cox.zph(ff)
  ## time-varying effect
  ff.s <- survSplit(as.formula(paste0("Surv(phecode_", phe, "_time,","phecode_", phe,"_bin) ~ ", x, " + ", adj)), ukb.dat, cut=c(2,5,10), episode = "tgroup")
  ## fit the model (need to redefine names of variables)
  ff.s <- coxph(as.formula(paste0("Surv(tstart, phecode_", phe, "_time,","phecode_", phe,"_bin) ~ ", x, ":strata(tgroup) + ", adj)), ff.s, ties = "breslow")
  ## compute summary for storage
  ff   <- summary(ff)
  ff.s <- summary(ff.s)
  ## store what is needed
  return(data.frame(id=x, phecode=phe, beta=ff$coefficients[1,1], se=ff$coefficients[1,3], pval=ff$coefficients[1,5], nevent=ff$nevent, nall=sum(!is.na(ukb.dat[, ..x])),
                    p.resid.prot=ff.p$table[x, 3], p.resid.overall=ff.p$table["GLOBAL", 3],
                    beta.2yr=ff.s$coefficients[paste0(x,":strata(tgroup)tgroup=1"),1], se.2yr=ff.s$coefficients[paste0(x,":strata(tgroup)tgroup=1"),3], pval.2yr=ff.s$coefficients[paste0(x,":strata(tgroup)tgroup=1"),5],
                    beta.5yr=ff.s$coefficients[paste0(x,":strata(tgroup)tgroup=2"),1], se.5yr=ff.s$coefficients[paste0(x,":strata(tgroup)tgroup=2"),3], pval.5yr=ff.s$coefficients[paste0(x,":strata(tgroup)tgroup=2"),5],
                    beta.10yr=ff.s$coefficients[paste0(x,":strata(tgroup)tgroup=3"),1], se.10yr=ff.s$coefficients[paste0(x,":strata(tgroup)tgroup=3"),3], pval.10yr=ff.s$coefficients[paste0(x,":strata(tgroup)tgroup=3"),5],
                    beta.18yr=ff.s$coefficients[paste0(x,":strata(tgroup)tgroup=4"),1], se.18yr=ff.s$coefficients[paste0(x,":strata(tgroup)tgroup=4"),3], pval.18yr=ff.s$coefficients[paste0(x,":strata(tgroup)tgroup=4"),5]))
  
  
}, mc.cores = 12) 
## combine into one
res.surv <- do.call(rbind, res.surv)

## store the results
fwrite(res.surv, paste("output/res.surv.v2", phe, "txt", sep="."), sep="\t", row.names = F, na = NA)
