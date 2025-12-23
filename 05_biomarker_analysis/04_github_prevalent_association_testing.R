#!/usr/bin/env Rscript

## script to run association testing between proteins and phecodes (existing disease)
## Maik Pietzner 02/07/2025
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
ukb.cov  <- fread("input/UKB.covariates.20240624.txt")
## proteins
ukb.prot <- fread("input/UKB.proteins.20240624.txt")
## relevant phecode
ukb.phe  <- fread("input/UKB.phecodes.prevalence.20250702.txt", select = c("f.eid", paste0("prev_", phe)))
## combine everything
ukb.dat  <- merge(ukb.cov, ukb.phe)
ukb.dat  <- merge(ukb.dat, ukb.prot, by.x="f.eid", by.y="eid")
## delete what is no longer needed
rm(ukb.cov); rm(ukb.prot); rm(ukb.phe); gc(reset=T)

## import protein label
lab.prot <- fread("input/label.Olink.proteins.20240624.txt")

#----------------------------#
##--    run associations  --##
#----------------------------#

## decide whether or not the outcome is sex-specific
if(cohort != "Both"){
  ukb.dat <- ukb.dat[ sex == cohort ]
  ## formula for adjustment
  adj      <- "age + I(age^2) + rcs(month_blood, c(1,6,12)) + rcs(time_blood.num, c(9,14,20)) + rcs(fast.0, c(.5,3,20)) + rcs(sample_age, c(11,12.7,15))"
}else{
  ## formula for adjustment
  adj      <- "age * sex + I(age^2) + rcs(month_blood, c(1,6,12)) + rcs(time_blood.num, c(9,14,20)) + rcs(fast.0, c(.5,3,20)) + rcs(sample_age, c(11,12.7,15))"
}

## use model matrix to ease coding
tmp           <- model.matrix.lm(as.formula(paste0("~ ", adj)), ukb.dat, na.action = "na.pass")
## edit colnames
colnames(tmp) <- gsub("[[:punct:]]| ", "_", colnames(tmp))
## define new adjustment set
adj           <- paste(colnames(tmp)[-1], collapse = " + ")
## create new data set
ukb.dat       <- data.table(ukb.dat[, c("f.eid", paste0("prev_", phe), lab.prot$id), with=F], tmp[,-1])
rm(tmp); gc(reset=T)

## intermediate output
cat("run regression with ", nrow(ukb.dat), " samples in total\n")

## do in parallel
registerDoMC(10)

## run testing
res.prev <- mclapply(lab.prot$id, function(x){
  
  print(x)
  
  ## run model: implement cox-prop test and possible step function for time intervals
  ff   <- summary(glm(paste0("prev_", phe, " ~ ", x, " + ", adj), ukb.dat, family = binomial(link = "logit")))
  ## store what is needed
  return(data.frame(id=x, phecode=phe, 
                    beta=ff$coefficients[2,1], 
                    se=ff$coefficients[2,2], 
                    pval=ff$coefficients[2,4], 
                    nevent=sum(ukb.dat[, paste0("prev_", phe), with = F] == 1), 
                    nall=sum(!is.na(ukb.dat[, ..x]))))
  
  
  
  
}, mc.cores = 10) 
## combine into one
res.prev <- rbindlist(res.prev)

## store the results
fwrite(res.prev, paste("output/res.prevalent.crude", phe, "txt", sep="."), sep="\t", row.names = F, na = NA)

